"""
Equation numbers are from
Ogden, B. F. L. (n.d.). Parameter Estimation for a Conceptual Functional Equivalen (CFE) Formulation of the National Water Model.
"""

import time
import numpy as np
import pandas as pd
import sys
import json
import matplotlib.pyplot as plt

plt.rcParams.update({"font.size": 20})
plt.rcParams["figure.figsize"] = [16, 4]
from cfe import CFE
import spotpy
from tqdm import tqdm


class BMI_CFE:
    def __init__(self, cfg_file=None, verbose=False):
        """Create a Bmi CFE model that is ready for initialization."""
        super(BMI_CFE, self).__init__()
        self._values = {}
        self._var_loc = "node"
        self._var_grid_id = 0
        self._start_time = 0.0
        self._end_time = np.finfo("d").max

        # these need to be initialized here as scale_output() called in update()
        self.streamflow_cms = 0.0
        self.streamflow_fms = 0.0
        self.surface_runoff_m = 0.0

        # ----------------------------------------------
        # Required, static attributes of the model
        # ----------------------------------------------
        self._att_map = {
            "model_name": "Conceptual Functional Equivalent (CFE)",
            "version": "1.0",
            "author_name": "Jonathan Martin Frame",
            "grid_type": "scalar",
            "time_step_size": 1,
            "time_units": "1 hour",
        }

        # ---------------------------------------------
        # Input variable names (CSDMS standard names)
        # ---------------------------------------------
        self._input_var_names = [
            "atmosphere_water__time_integral_of_precipitation_mass_flux",
            "water_potential_evaporation_flux",
        ]

        # ---------------------------------------------
        # Output variable names (CSDMS standard names)
        # ---------------------------------------------
        self._output_var_names = [
            "land_surface_water__runoff_depth",
            "land_surface_water__runoff_volume_flux",
            "DIRECT_RUNOFF",
            "GIUH_RUNOFF",
            "NASH_LATERAL_RUNOFF",
            "DEEP_GW_TO_CHANNEL_FLUX",
            "SOIL_CONCEPTUAL_STORAGE",
        ]

        # ------------------------------------------------------
        # Create a Python dictionary that maps CSDMS Standard
        # Names to the model's internal variable names.
        # This is going to get long,
        #     since the input variable names could come from any forcing...
        # ------------------------------------------------------
        self._var_name_units_map = {
            "land_surface_water__runoff_volume_flux": ["streamflow_cfs", "ft3 s-1"],
            "land_surface_water__runoff_depth": ["total_discharge", "m"],
            # --------------   Dynamic inputs --------------------------------
            "atmosphere_water__time_integral_of_precipitation_mass_flux": [
                "timestep_rainfall_input_m",
                "kg m-2",
            ],
            "water_potential_evaporation_flux": ["potential_et_m_per_s", "m s-1"],
            "DIRECT_RUNOFF": ["surface_runoff_depth_m", "m"],
            "GIUH_RUNOFF": ["flux_giuh_runoff_m", "m"],
            "NASH_LATERAL_RUNOFF": ["flux_nash_lateral_runoff_m", "m"],
            "DEEP_GW_TO_CHANNEL_FLUX": ["flux_from_deep_gw_to_chan_m", "m"],
            "SOIL_CONCEPTUAL_STORAGE": ["soil_reservoir['storage_m']", "m"],
        }

        # ------------------------------------------------------------
        # this is the bmi configuration file
        self.cfg_file = cfg_file
        self.verbose = verbose

    # __________________________________________________________________
    # __________________________________________________________________
    # BMI: Model Control Function
    def initialize(self, current_time_step=0):
        self.current_time_step = current_time_step

        # ----- Create some lookup tabels from the long variable names --------#
        self._var_name_map_long_first = {
            long_name: self._var_name_units_map[long_name][0]
            for long_name in self._var_name_units_map.keys()
        }
        self._var_name_map_short_first = {
            self._var_name_units_map[long_name][0]: long_name
            for long_name in self._var_name_units_map.keys()
        }
        self._var_units_map = {
            long_name: self._var_name_units_map[long_name][1]
            for long_name in self._var_name_units_map.keys()
        }

        # -------------- Initalize all the variables --------------------------#
        # -------------- so that they'll be picked up with the get functions --#
        for long_var_name in list(self._var_name_units_map.keys()):
            # ---------- All the variables are single values ------------------#
            # ---------- so just set to zero for now.        ------------------#
            self._values[long_var_name] = 0
            setattr(self, self.get_var_name(long_var_name), 0)

        ############################################################
        # ________________________________________________________ #
        # GET VALUES FROM CONFIGURATION FILE.                      #
        self.config_from_json()  #

        # ________________________________________________
        # The configuration should let the BMI know what mode to run in (framework vs standalone)
        # If it is stand alone, then load in the forcing and read the time from the forcig file
        if self.stand_alone == 1:
            self.load_forcing_file()
            self.current_time = pd.Timestamp(
                self.forcing_data["time"][self.current_time_step]
            )

        # ________________________________________________
        # In order to check mass conservation at any time
        self.reset_volume_tracking()

        # ________________________________________________
        # initialize simulation constants
        atm_press_Pa = 101325.0
        unit_weight_water_N_per_m3 = 9810.0

        # ________________________________________________
        # Time control
        self.time_step_size = 3600
        self.timestep_h = self.time_step_size / 3600.0
        self.timestep_d = self.timestep_h / 24.0
        self.current_time_step = 0
        self.current_time = pd.Timestamp(year=1970, month=1, day=1, hour=0)

        # ________________________________________________
        # Inputs
        self.timestep_rainfall_input_m = 0
        self.potential_et_m_per_s = 0

        # ________________________________________________
        # calculated flux variables
        self.flux_overland_m = (
            0  # surface runoff that goes through the GIUH convolution process
        )
        self.flux_perc_m = 0  # flux from soil to deeper groundwater reservoir
        self.flux_lat_m = 0  # lateral flux in the subsurface to the Nash cascade
        self.flux_from_deep_gw_to_chan_m = (
            0  # flux from the deep reservoir into the channels
        )
        self.gw_reservoir_storage_deficit_m = (
            0  # the available space in the conceptual groundwater reservoir
        )
        self.primary_flux = 0  # temporary vars.
        self.secondary_flux = 0  # temporary vars.
        self.total_discharge = 0
        self.diff_infilt = 0
        self.diff_perc = 0

        # ________________________________________________
        # Evapotranspiration
        self.potential_et_m_per_timestep = 0
        self.actual_et_m_per_timestep = 0
        self.reduced_potential_et_m_per_timestep = 0
        self.actual_et_from_rain_m_per_timestep = 0
        self.actual_et_from_soil_m_per_timestep = 0

        # ________________________________________________________
        # Set these values now that we have the information from the configuration file.
        self.runoff_queue_m_per_timestep = np.zeros(len(self.giuh_ordinates))
        self.num_giuh_ordinates = len(self.runoff_queue_m_per_timestep)
        self.num_lateral_flow_nash_reservoirs = self.nash_storage.shape[0]

        # ________________________________________________
        # Local values to be used in setting up soil reservoir
        # This seems to be the key component in the model

        # trigger_z_m = 0.5                           # Added to calibration parameter

        # field_capacity_atm_press_fraction = 0.33     # Added to calibration parameter

        self.D_noahMP = 2  # Soil depth used for threshold calculation (NoahMP grid cell vertical grid size)
        self.trigger_z_m = self.trigger_z_fact * self.D_noahMP

        H_water_table_m = (
            self.field_capacity_atm_press_fraction
            * atm_press_Pa
            / unit_weight_water_N_per_m3
        )  # Lateral flow function 1&2, Equation 3 (Ogden's document)

        soil_water_content_at_field_capacity = self.soil_params["smcmax"] * np.power(
            H_water_table_m / self.soil_params["satpsi"], (1.0 / self.soil_params["bb"])
        )  # Lateral flow function 1&2, Equation 3 (Ogden's document), but this equaton is not used in later process

        Omega = (
            H_water_table_m - self.trigger_z_m
        )  # Lateral flow function 3 (Ogden's document). Soil water storage in the soil column, when the lowest soil discretization is equal to field capacity.

        if Omega < 0:
            integral_lower_lim_m = 0
            integral_upper_lim_m = self.D_noahMP - self.trigger_z_m + H_water_table_m
            additional_term = self.trigger_z_m - H_water_table_m
        else:
            integral_lower_lim_m = Omega
            integral_upper_lim_m = Omega + self.D_noahMP
            additional_term = 0

        lower_lim = np.power(
            integral_lower_lim_m, (1.0 - 1.0 / self.soil_params["bb"])
        ) / (
            1.0 - 1.0 / self.soil_params["bb"]
        )  # Lower limit of the integral in equation 4 (Ogden's document).

        upper_lim = np.power(
            integral_upper_lim_m, (1.0 - 1.0 / self.soil_params["bb"])
        ) / (
            1.0 - 1.0 / self.soil_params["bb"]
        )  # Upper limit of the integral in equation 4 (Ogden's document).

        storage_thresh_pow_term = np.power(
            1.0 / self.soil_params["satpsi"], (-1.0 / self.soil_params["bb"])
        )  # This seems to be redundant? Or for needed for descritization?

        lim_diff = (
            upper_lim - lower_lim
        )  # Integral term of the equation 4 (Ogden's document).

        field_capacity_storage_threshold_m = (
            self.soil_params["smcmax"] * storage_thresh_pow_term * lim_diff
        ) + self.soil_params[
            "smcmax"
        ] * additional_term  # Equation 4 (and probably 5?) (Ogden's document).

        # ________________________________________________
        # lateral flow function parameters
        assumed_near_channel_water_table_slope = 0.01  # [L/L]
        lateral_flow_threshold_storage_m = field_capacity_storage_threshold_m  # Equation 4 (and probably 5?) (Ogden's document).
        #         lateral_flow_linear_reservoir_constant = 2.0 * assumed_near_channel_water_table_slope * \     # Not used      # Equation 9 (Ogden's document).
        #                                                  self.soil_params['mult'] * NWM_soil_params.satdk * \ # Not used
        #                                                  self.soil_params['D'] * drainage_density_km_per_km2  # Not used
        #         lateral_flow_linear_reservoir_constant *= 3600.0                                              # Not used
        self.soil_reservoir_storage_deficit_m = 0

        # ________________________________________________
        # Subsurface reservoirs
        self.gw_reservoir = {
            "is_exponential": True,
            "storage_max_m": self.max_gw_storage,
            # GW primary reservoir params
            "coeff_primary": self.Cgw,
            "exponent_primary": self.expon,
            "storage_threshold_primary_m": 0.0,
            # Currently one storage. The secoundary storage is turned off
            "coeff_secondary": 0.0,
            "exponent_secondary": 1.0,
            "storage_threshold_secondary_m": 0.0,
        }
        self.gw_reservoir["storage_m"] = 1.0e-06  # 1.0e-06  # $self.gw_reservoir[
        # "storage_max_m"
        # ]  * 0.5 # Start from the half groundwater reservoir
        self.volstart += self.gw_reservoir["storage_m"]
        self.vol_in_gw_start = self.gw_reservoir["storage_m"]

        self.soil_reservoir = {
            "is_exponential": False,
            "storage_max_m": self.soil_params["smcmax"]
            * self.soil_params[
                "D"
            ],  # Schaake Partitioning function 1 (Ogden's document)
            "coeff_primary": self.soil_params["satdk"]
            * self.soil_params["slop"]
            * self.time_step_size,  # Controls percolation to GW, Equation 11
            "exponent_primary": 1,  # Controls percolation to GW, FIXED to 1 based on Equation 11
            "storage_threshold_primary_m": field_capacity_storage_threshold_m,  # Equation 4 (and probably 5?) (Ogden's document).
            "coeff_secondary": self.K_lf,  # Controls lateral flow
            "exponent_secondary": 1,  # Controls lateral flow, FIXED to 1 based on schematics in the Ogden's document
            "storage_threshold_secondary_m": lateral_flow_threshold_storage_m,
        }  # Equation 4 (and probably 5?) (Ogden's document).
        self.soil_reservoir["storage_m"] = self.soil_reservoir[
            "storage_max_m"
        ]  # Start from the maximum soil reservoir
        self.volstart += self.soil_reservoir["storage_m"]
        self.vol_soil_start = self.soil_reservoir["storage_m"]

        # print(self.soil_reservoir['coeff_primary'])
        # print(self.soil_reservoir['coeff_secondary'])
        # ________________________________________________
        # Schaake
        # self.refkdt = 3.0       # Added to calibration parameters
        Ks_ref = 2.0e-06  # [m/s]
        self.Schaake_adjusted_magic_constant_by_soil_type = (
            self.refkdt * self.soil_params["satdk"] / Ks_ref
        )  # Schaake partitioning function 2  (Ogden's document)
        self.Schaake_output_runoff_m = 0
        self.infiltration_depth_m = 0

        # ________________________________________________
        # Nash cascade
        self.K_nash = self.K_nash

        # ----------- The output is area normalized, this is needed to un-normalize it
        #                         mm->m                             km2 -> m2          hour->s
        self.output_factor_cms = (
            (1 / 1000)
            * (self.catchment_area_km2 * 1000 * 1000)
            * (1 / self.time_step_size)
        )

        ####################################################################
        # ________________________________________________________________ #
        # ________________________________________________________________ #
        # CREATE AN INSTANCE OF THE CONCEPTUAL FUNCTIONAL EQUIVALENT MODEL #
        self.cfe_model = CFE()
        # ________________________________________________________________ #
        # ________________________________________________________________ #
        ####################################################################

    # __________________________________________________________________________________________________________
    # __________________________________________________________________________________________________________
    # BMI: Model Control Function
    def update(self):
        self.cfe_model.run_cfe(self)
        self.scale_output()

    # __________________________________________________________________________________________________________
    # __________________________________________________________________________________________________________
    # BMI: Model Control Function
    def update_until(self, until, verbose=True):
        for i in range(self.current_time_step, until):
            self.cfe_model.run_cfe(self)
            self.scale_output()
            if verbose:
                print("total discharge: {}".format(self.total_discharge))
                print("at time: {}".format(self.current_time))

    # __________________________________________________________________________________________________________
    # __________________________________________________________________________________________________________
    # BMI: Model Control Function
    def finalize(self, print_mass_balance=False):
        self.finalize_mass_balance(verbose=print_mass_balance)
        self.reset_volume_tracking()

        """Finalize model."""
        self.cfe_model = None
        self.cfe_state = None

    # ________________________________________________
    # Mass balance tracking
    def reset_volume_tracking(self):
        """Reset the mass balance completely after a run. Prepare the model for a new run"""

        # Watershed mass balance
        self.volstart = 0
        self.volin = 0
        self.volout = 0
        self.volend = 0

        # ET
        self.vol_PET = 0
        self.vol_et_to_atm = 0
        self.vol_et_from_soil = 0
        self.vol_et_from_rain = 0

        # Partitioning
        self.vol_sch_runoff = 0
        self.vol_sch_runoff_IOF = 0
        self.vol_sch_runoff_SOF = 0
        self.vol_sch_infilt = 0

        # GIUH
        self.vol_out_giuh = 0
        self.vol_end_giuh = 0

        # Soil
        self.vol_soil_start = 0
        self.vol_to_soil = 0
        self.vol_soil_to_lat_flow = 0
        self.vol_soil_to_gw = 0
        self.vol_soil_end = 0

        # GW
        self.vol_to_gw = 0
        self.vol_from_gw = 0
        self.vol_to_gw_start = 0
        self.vol_to_gw_end = 0

        # Nash
        self.vol_in_nash = 0
        self.vol_in_nash_end = 0
        self.vol_out_nash = 0

        return

    # ________________________________________________________
    def config_from_json(self):
        with open(self.cfg_file) as data_file:
            data_loaded = json.load(data_file)

        # ___________________________________________________
        # MANDATORY CONFIGURATIONS
        self.forcing_file = data_loaded["forcing_file"]
        self.catchment_area_km2 = data_loaded["catchment_area_km2"]

        # SOIL PARAMETERS
        self.trigger_z_fact = data_loaded["trigger_z_fact"]
        self.field_capacity_atm_press_fraction = data_loaded["alpha_fc"]
        self.allow_percolation_below_threshold = data_loaded[
            "allow_percolation_below_threshold"
        ]
        self.soil_params = {}
        self.soil_params["bb"] = data_loaded["soil_params"]["bb"]
        self.soil_params["D"] = data_loaded["soil_params"]["D"]
        self.soil_params["satdk"] = data_loaded["soil_params"]["satdk"]
        self.soil_params["satpsi"] = data_loaded["soil_params"]["satpsi"]
        self.soil_params["slop"] = data_loaded["soil_params"]["slop"]
        self.soil_params["smcmax"] = data_loaded["soil_params"]["smcmax"]
        self.soil_params["wltsmc"] = data_loaded["soil_params"]["wltsmc"]
        # GROUNDWATER PARAMETERS
        self.max_gw_storage = data_loaded["max_gw_storage"]
        self.Cgw = data_loaded["Cgw"]
        self.expon = data_loaded["expon"]

        # Schaake parameters
        self.refkdt = data_loaded["refkdt"]
        # self.lksatfac = data_loaded["lksatfac"]

        #
        self.K_lf = data_loaded["K_lf"]
        # self.K_lf = (
        #     0.02
        #     * self.lksatfac
        #     * data_loaded["soil_params"]["satdk"]
        #     * data_loaded["soil_params"]["D"]
        #     * data_loaded["dd"]
        # )  # Equation 11 in the Ogden's document
        self.K_nash = data_loaded["K_nash"]
        self.nash_storage = np.zeros(int(data_loaded["num_nash_storage"]))
        self.giuh_ordinates = np.array(data_loaded["giuh_ordinates"])

        # ___________________________________________________
        # OPTIONAL CONFIGURATIONS
        if "stand_alone" in data_loaded.keys():
            self.stand_alone = data_loaded["stand_alone"]
        if "forcing_file" in data_loaded.keys():
            self.reads_own_forcing = True
            self.forcing_file = data_loaded["forcing_file"]
        if "unit_test" in data_loaded.keys():
            self.unit_test = data_loaded["unit_test"]
            self.compare_results_file = data_loaded["compare_results_file"]

        return

    # ________________________________________________________
    def reset_volume_tracking_after_warmup(self):
        """Reset the mass balance for the warm-up period"""

        # Watershed mass balance
        self.volstart = (
            self.soil_reservoir["storage_m"] + self.gw_reservoir["storage_m"]
        )
        self.volin = 0
        self.volout = 0

        # ET
        self.vol_PET = 0
        self.vol_et_to_atm = 0
        self.vol_et_from_soil = 0
        self.vol_et_from_rain = 0

        # Partitioning
        self.vol_sch_runoff = 0
        self.vol_sch_runoff_IOF = 0
        self.vol_sch_runoff_SOF = 0
        self.vol_sch_infilt = 0

        # GIUH
        self.vol_out_giuh = 0

        # Soil
        self.vol_soil_start = self.soil_reservoir["storage_m"]
        self.vol_to_soil = 0
        self.vol_soil_to_lat_flow = 0
        self.vol_soil_to_gw = 0

        # GW
        self.vol_in_gw_start = self.gw_reservoir["storage_m"]
        self.vol_to_gw = 0
        self.vol_from_gw = 0
        self.vol_to_gw_start = 0
        self.vol_to_gw_end = 0

        # Nash
        self.vol_in_nash = 0
        self.vol_out_nash = 0

    # ________________________________________________________
    def finalize_mass_balance(self, verbose=True):
        self.volend = self.soil_reservoir["storage_m"] + self.gw_reservoir["storage_m"]
        self.vol_in_gw_end = self.gw_reservoir["storage_m"]

        # the GIUH queue might have water in it at the end of the simulation, so sum it up.
        self.vol_end_giuh = np.sum(self.runoff_queue_m_per_timestep)
        self.vol_in_nash_end = np.sum(self.nash_storage)

        self.vol_soil_end = self.soil_reservoir["storage_m"]

        self.global_residual = (
            self.volstart + self.volin - self.volout - self.volend - self.vol_end_giuh
        )  # (?) - vol_in_nash_end? --- This is already excluded from the system.
        self.schaake_residual = (
            self.volin
            - self.vol_et_from_rain
            - self.vol_sch_runoff
            - self.vol_sch_infilt
        )
        self.giuh_residual = self.vol_sch_runoff - self.vol_out_giuh - self.vol_end_giuh
        self.soil_residual = (
            self.vol_soil_start
            - self.vol_soil_end
            + self.vol_sch_infilt
            - self.vol_soil_to_lat_flow
            - self.vol_soil_to_gw
            - self.vol_et_from_soil
        )
        self.nash_residual = self.vol_in_nash - self.vol_out_nash - self.vol_in_nash_end
        self.gw_residual = (
            self.vol_in_gw_start
            + self.vol_to_gw
            - self.vol_from_gw
            - self.vol_in_gw_end
        )

        # Discharge
        self.cumQ_out = self.vol_out_giuh + self.vol_from_gw + self.vol_out_nash
        self.runoff_ratio = self.cumQ_out / self.volin

        if verbose:
            print("\nGLOBAL MASS BALANCE")
            print("      initial volume: {:8.4f}".format(self.volstart))
            print("        volume input: {:8.4f}".format(self.volin))
            print("       volume output: {:8.4f}".format(self.volout))
            print("        final volume: {:8.4f}".format(self.volend))
            print("            residual: {:6.4e}".format(self.global_residual))

            print("\n Discharge details")
            print("Cumulative discharge: {:1.2f}".format(self.cumQ_out))
            print("        Runoff ratio: {:1.2f}".format(self.runoff_ratio))
            print(
                "      Surface runoff: {:1.2f}".format(
                    self.vol_out_giuh / self.cumQ_out
                )
            )
            print(
                "  Deep GW to channel: {:1.2f}".format(self.vol_from_gw / self.cumQ_out)
            )
            print(
                "        Lateral flow: {:1.2f}".format(
                    self.vol_out_nash / self.cumQ_out
                )
            )

            print("\n AET & PET")
            print("      volume PET: {:8.4f}".format(self.vol_PET))
            print("      volume AET: {:8.4f}".format(self.vol_et_to_atm))
            print("ET from rainfall: {:8.4f}".format(self.vol_et_from_rain))
            print("    ET from soil: {:8.4f}".format(self.vol_et_from_soil))

            print("\nSCHAAKE MASS BALANCE")
            print("    volume input: {:8.4f}".format(self.volin))
            print("ET from rainfall: {:8.4f}".format(self.vol_et_from_rain))
            print("  surface runoff: {:8.4f}".format(self.vol_sch_runoff))
            print("             IOF: {:8.4f}".format(self.vol_sch_runoff_IOF))
            print("             SOF: {:8.4f}".format(self.vol_sch_runoff_SOF))
            print("    infiltration: {:8.4f}".format(self.vol_sch_infilt))
            print("schaake residual: {:6.4e}".format(self.schaake_residual))

            print("\nGIUH MASS BALANCE")
            print("  vol. into giuh: {:8.4f}".format(self.vol_sch_runoff))
            print("   vol. out giuh: {:8.4f}".format(self.vol_out_giuh))
            print(" vol. end giuh q: {:8.4f}".format(self.vol_end_giuh))
            print("   giuh residual: {:6.4e}".format(self.giuh_residual))

            print("\nSOIL WATER CONCEPTUAL RESERVOIR MASS BALANCE")
            print("   init soil vol: {:8.4f}".format(self.vol_soil_start))
            print("  vol. into soil: {:8.4f}".format(self.vol_sch_infilt))
            print("vol.soil2latflow: {:8.4f}".format(self.vol_soil_to_lat_flow))
            print(" vol. soil to gw: {:8.4f}".format(self.vol_soil_to_gw))
            print(" vol. soil to ET: {:8.4f}".format(self.vol_et_from_soil))
            print(" final vol. soil: {:8.4f}".format(self.vol_soil_end))
            print("vol. soil resid.: {:6.4e}".format(self.soil_residual))

            print("\nNASH CASCADE CONCEPTUAL RESERVOIR MASS BALANCE")
            print("    vol. to nash: {:8.4f}".format(self.vol_in_nash))
            print("  vol. from nash: {:8.4f}".format(self.vol_out_nash))
            print(" final vol. nash: {:8.4f}".format(self.vol_in_nash_end))
            print("nash casc resid.: {:6.4e}".format(self.nash_residual))

            print("\nGROUNDWATER CONCEPTUAL RESERVOIR MASS BALANCE")
            print("init gw. storage: {:8.4f}".format(self.vol_in_gw_start))
            print("       vol to gw: {:8.4f}".format(self.vol_to_gw))
            print("     vol from gw: {:8.4f}".format(self.vol_from_gw))
            print("final gw.storage: {:8.4f}".format(self.vol_in_gw_end))
            print("    gw. residual: {:6.4e}".format(self.gw_residual))

        return

    # ________________________________________________________
    def load_forcing_file(self):
        self.forcing_data = pd.read_csv(self.forcing_file)

    # ________________________________________________________
    def load_unit_test_data(self):
        self.unit_test_data = pd.read_csv(self.compare_results_file)
        self.unit_test_data["Soil Moisture Content"] = self.unit_test_data[
            "Soil Moisture Content"
        ]
        self.cfe_output_data = pd.DataFrame().reindex_like(self.unit_test_data)
        return self.unit_test_data

    # ------------------------------------------------------------
    def scale_output(self):
        self.surface_runoff_m = self.total_discharge
        self._values["land_surface_water__runoff_depth"] = self.surface_runoff_m / 1000
        self.streamflow_cms = (
            self._values["land_surface_water__runoff_depth"] * self.output_factor_cms
        )

        self._values["land_surface_water__runoff_volume_flux"] = self.streamflow_cms * (
            1 / 35.314
        )

        self._values["DIRECT_RUNOFF"] = self.surface_runoff_depth_m
        self._values["GIUH_RUNOFF"] = self.flux_giuh_runoff_m
        self._values["NASH_LATERAL_RUNOFF"] = self.flux_nash_lateral_runoff_m
        self._values["DEEP_GW_TO_CHANNEL_FLUX"] = self.flux_from_deep_gw_to_chan_m
        self._values["SOIL_CONCEPTUAL_STORAGE"] = self.soil_reservoir["storage_m"]

    # ----------------------------------------------------------------------------
    def initialize_forcings(self):
        for forcing_name in self.cfg_train["dynamic_inputs"]:
            setattr(self, self._var_name_map_short_first[forcing_name], 0)

    # -------------------------------------------------------------------
    # -------------------------------------------------------------------
    # BMI: Model Information Functions
    # -------------------------------------------------------------------
    # -------------------------------------------------------------------

    def get_attribute(self, att_name):
        try:
            return self._att_map[att_name.lower()]
        except:
            print(" ERROR: Could not find attribute: " + att_name)

    # --------------------------------------------------------
    # Note: These are currently variables needed from other
    #       components vs. those read from files or GUI.
    # --------------------------------------------------------
    def get_input_var_names(self):
        return self._input_var_names

    def get_output_var_names(self):
        return self._output_var_names

    # ------------------------------------------------------------
    def get_component_name(self):
        """Name of the component."""
        return self.get_attribute("model_name")  # JG Edit

    # ------------------------------------------------------------
    def get_input_item_count(self):
        """Get names of input variables."""
        return len(self._input_var_names)

    # ------------------------------------------------------------
    def get_output_item_count(self):
        """Get names of output variables."""
        return len(self._output_var_names)

    # ------------------------------------------------------------
    def get_value(self, var_name):
        """Copy of values.
        Parameters
        ----------
        var_name : str
            Name of variable as CSDMS Standard Name.
        dest : ndarray
            A numpy array into which to place the values.
        Returns
        -------
        array_like
            Copy of values.
        """
        return self.get_value_ptr(var_name)

    # -------------------------------------------------------------------
    def get_value_ptr(self, var_name):
        """Reference to values.
        Parameters
        ----------
        var_name : str
            Name of variable as CSDMS Standard Name.
        Returns
        -------
        array_like
            Value array.
        """
        return self._values[var_name]

    # -------------------------------------------------------------------
    # -------------------------------------------------------------------
    # BMI: Variable Information Functions
    # -------------------------------------------------------------------
    # -------------------------------------------------------------------
    def get_var_name(self, long_var_name):
        return self._var_name_map_long_first[long_var_name]

    # -------------------------------------------------------------------
    def get_var_units(self, long_var_name):
        return self._var_units_map[long_var_name]

    # -------------------------------------------------------------------
    def get_var_type(self, long_var_name):
        """Data type of variable.

        Parameters
        ----------
        var_name : str
            Name of variable as CSDMS Standard Name.

        Returns
        -------
        str
            Data type.
        """
        # JG Edit
        return self.get_value_ptr(long_var_name)  # .dtype

    # ------------------------------------------------------------
    def get_var_grid(self, name):
        # JG Edit
        # all vars have grid 0 but check if its in names list first
        if name in (self._output_var_names + self._input_var_names):
            return self._var_grid_id

    # ------------------------------------------------------------
    def get_var_itemsize(self, name):
        #        return np.dtype(self.get_var_type(name)).itemsize
        return np.array(self.get_value(name)).itemsize

    # ------------------------------------------------------------
    def get_var_location(self, name):
        # JG Edit
        # all vars have location node but check if its in names list first
        if name in (self._output_var_names + self._input_var_names):
            return self._var_loc

    # -------------------------------------------------------------------
    # JG Note: what is this used for?
    def get_var_rank(self, long_var_name):
        return np.int16(0)

    # -------------------------------------------------------------------
    def get_start_time(self):
        return self._start_time  # JG Edit

    # -------------------------------------------------------------------
    def get_end_time(self):
        return self._end_time  # JG Edit

    # -------------------------------------------------------------------
    def get_current_time(self):
        return self.current_time

    # -------------------------------------------------------------------
    def get_time_step(self):
        return self.get_attribute("time_step_size")  # JG: Edit

    # -------------------------------------------------------------------
    def get_time_units(self):
        return self.get_attribute("time_units")

    # -------------------------------------------------------------------
    def set_value(self, var_name, value):
        """Set model values.

        Parameters
        ----------
        var_name : str
            Name of variable as CSDMS Standard Name.
        src : array_like
              Array of new values.
        """
        setattr(self, self.get_var_name(var_name), value)
        self._values[var_name] = value

    # ------------------------------------------------------------
    def set_value_at_indices(self, name, inds, src):
        """Set model values at particular indices.
        Parameters
        ----------
        var_name : str
            Name of variable as CSDMS Standard Name.
        src : array_like
            Array of new values.
        indices : array_like
            Array of indices.
        """
        # JG Note: TODO confirm this is correct. Get/set values
        #        val = self.get_value_ptr(name)
        #        val.flat[inds] = src

        # JMFrame: chances are that the index will be zero, so let's include that logic
        if np.array(self.get_value(name)).flatten().shape[0] == 1:
            self.set_value(name, src)
        else:
            # JMFrame: Need to set the value with the updated array with new index value
            val = self.get_value_ptr(name)
            for i in inds.shape:
                val.flatten()[inds[i]] = src[i]
            self.set_value(name, val)

    # ------------------------------------------------------------
    def get_var_nbytes(self, long_var_name):
        """Get units of variable.
        Parameters
        ----------
        var_name : str
            Name of variable as CSDMS Standard Name.
        Returns
        -------
        int
            Size of data array in bytes.
        """
        # JMFrame NOTE: Had to import sys for this function
        return sys.getsizeof(self.get_value_ptr(long_var_name))

    # ------------------------------------------------------------
    def get_value_at_indices(self, var_name, dest, indices):
        """Get values at particular indices.
        Parameters
        ----------
        var_name : str
            Name of variable as CSDMS Standard Name.
        dest : ndarray
            A numpy array into which to place the values.
        indices : array_like
            Array of indices.
        Returns
        -------
        array_like
            Values at indices.
        """
        # JMFrame: chances are that the index will be zero, so let's include that logic
        if np.array(self.get_value(var_name)).flatten().shape[0] == 1:
            return self.get_value(var_name)
        else:
            val_array = self.get_value(var_name).flatten()
            return np.array([val_array[i] for i in indices])

    # JG Note: remaining grid funcs do not apply for type 'scalar'
    #   Yet all functions in the BMI must be implemented
    #   See https://bmi.readthedocs.io/en/latest/bmi.best_practices.html
    # ------------------------------------------------------------
    def get_grid_edge_count(self, grid):
        raise NotImplementedError("get_grid_edge_count")

    # ------------------------------------------------------------
    def get_grid_edge_nodes(self, grid, edge_nodes):
        raise NotImplementedError("get_grid_edge_nodes")

    # ------------------------------------------------------------
    def get_grid_face_count(self, grid):
        raise NotImplementedError("get_grid_face_count")

    # ------------------------------------------------------------
    def get_grid_face_edges(self, grid, face_edges):
        raise NotImplementedError("get_grid_face_edges")

    # ------------------------------------------------------------
    def get_grid_face_nodes(self, grid, face_nodes):
        raise NotImplementedError("get_grid_face_nodes")

    # ------------------------------------------------------------
    def get_grid_node_count(self, grid):
        raise NotImplementedError("get_grid_node_count")

    # ------------------------------------------------------------
    def get_grid_nodes_per_face(self, grid, nodes_per_face):
        raise NotImplementedError("get_grid_nodes_per_face")

    # ------------------------------------------------------------
    def get_grid_origin(self, grid_id, origin):
        raise NotImplementedError("get_grid_origin")

    # ------------------------------------------------------------
    def get_grid_rank(self, grid_id):
        # JG Edit
        # 0 is the only id we have
        if grid_id == 0:
            return 1

    # ------------------------------------------------------------
    def get_grid_shape(self, grid_id, shape):
        raise NotImplementedError("get_grid_shape")

    # ------------------------------------------------------------
    def get_grid_size(self, grid_id):
        # JG Edit
        # 0 is the only id we have
        if grid_id == 0:
            return 1

    # ------------------------------------------------------------
    def get_grid_spacing(self, grid_id, spacing):
        raise NotImplementedError("get_grid_spacing")

    # ------------------------------------------------------------
    def get_grid_type(self, grid_id=0):
        # JG Edit
        # 0 is the only id we have
        if grid_id == 0:
            return "scalar"

    # ------------------------------------------------------------
    def get_grid_x(self):
        raise NotImplementedError("get_grid_x")

    # ------------------------------------------------------------
    def get_grid_y(self):
        raise NotImplementedError("get_grid_y")

    # ------------------------------------------------------------
    def get_grid_z(self):
        raise NotImplementedError("get_grid_z")

    # ________________________________________________________
    def run_unit_test(
        self,
        plot_lims=list(range(1, 52600)),  # 31062 for MH # 52600 for LW
        plot=False,
        print_fluxes=False,
        warm_up=True,
        warmup_offset=365 * 24 * 2,
        warmup_iteration=1,  # 2 years by default
    ):
        self.load_forcing_file()
        self.load_unit_test_data()

        #################### Warmup ####################
        if warm_up:
            # Read forcing data for warmup
            self.forcing_data_warmup = self.forcing_data.iloc[0:warmup_offset].copy()

            previous_gw_storage = 1

            # Repeat warm-up period for 'warmup_iteration' times
            for i in range(0, warmup_iteration):
                # Run CFE by a timestep
                for t, precipitation_input in enumerate(
                    self.forcing_data_warmup["precip_rate"]
                ):
                    self.timestep_rainfall_input_m = precipitation_input
                    if "PET" in self.forcing_data_warmup.columns:
                        self.potential_et_m_per_timestep = self.forcing_data_warmup.loc[
                            t, "PET"
                        ]
                        self.potential_et_m_per_s = (
                            self.forcing_data_warmup.loc[t, "PET"] / self.time_step_size
                        )
                    self.cfe_model.run_cfe(self)

                diff = (
                    abs(previous_gw_storage - self.gw_reservoir["storage_m"])
                    / previous_gw_storage
                )
                if diff < 1.0e-02:
                    print(f"GW converged <1% after warm-up iteration {i}")
                    break

                previous_gw_storage = self.gw_reservoir["storage_m"]

            # Reset the volume tracking after warm-up. Leave resevoirs (the GW and soil reservoir etc.) as the current state
            self.reset_volume_tracking_after_warmup()

        ### Initilaize run ###
        output_time = [0] * len(self.unit_test_data)
        output_ts = [0] * len(self.unit_test_data)
        output_rainfall = [0] * len(self.unit_test_data)
        output_directrunoff = [0] * len(self.unit_test_data)
        output_GIUHrunoff = [0] * len(self.unit_test_data)
        output_lateralflow = [0] * len(self.unit_test_data)
        output_baseflow = [0] * len(self.unit_test_data)
        output_totdischarge = [0] * len(self.unit_test_data)
        output_flow = [0] * len(self.unit_test_data)
        output_SM = [0] * len(self.unit_test_data)

        output_gwstorage = [0] * len(self.unit_test_data)
        output_gwstorage_in = [0] * len(self.unit_test_data)
        output_gwstorage_out = [0] * len(self.unit_test_data)
        output_giuhstorage = [0] * len(self.unit_test_data)
        output_giuhstorage_in = [0] * len(self.unit_test_data)
        output_giuhstorage_out = [0] * len(self.unit_test_data)
        output_smstorage_in = [0] * len(self.unit_test_data)
        output_smstorage_out = [0] * len(self.unit_test_data)
        output_nashstorage = [0] * len(self.unit_test_data)
        output_nashstorage_in = [0] * len(self.unit_test_data)
        output_nashstorage_out = [0] * len(self.unit_test_data)

        self.current_time = pd.Timestamp(self.forcing_data["time"][0])

        for t, precipitation_input in enumerate(self.forcing_data["precip_rate"]):
            ### Set forcing ###
            self.timestep_rainfall_input_m = precipitation_input

            if "PET" in self.forcing_data.columns:
                self.potential_et_m_per_timestep = self.forcing_data.loc[t, "PET"]
                self.potential_et_m_per_s = (
                    self.forcing_data.loc[t, "PET"] / self.time_step_size
                )

            ### Record forcing ###
            self.current_time = pd.Timestamp(self.forcing_data["time"][t])
            output_time[t] = self.current_time
            output_ts[t] = self.current_time_step
            output_rainfall[t] = self.timestep_rainfall_input_m

            ### CFE run ###
            self.cfe_model.run_cfe(self)

            ### Record output ###
            output_directrunoff[t] = self.surface_runoff_depth_m
            output_GIUHrunoff[t] = self.flux_giuh_runoff_m
            output_lateralflow[t] = self.flux_nash_lateral_runoff_m
            output_baseflow[t] = self.flux_from_deep_gw_to_chan_m
            output_flow[t] = self.flux_Qout_m  #'flow' is in [m/timestep]
            output_totdischarge[t] = self.total_discharge  # 'totdischarge' is in [cms]
            output_SM[t] = self.soil_reservoir["storage_m"]

            output_gwstorage[t] = self.gw_reservoir["storage_m"]
            output_gwstorage_in[
                t
            ] = self.flux_perc_m  # this is percolation minus runoff
            output_gwstorage_out[t] = self.flux_from_deep_gw_to_chan_m + self.diff_perc

            output_giuhstorage[t] = np.sum(self.runoff_queue_m_per_timestep)
            output_giuhstorage_in[t] = (
                self.surface_runoff_depth_m + self.diff_perc + self.diff_infilt
            )
            output_giuhstorage_out[t] = self.flux_giuh_runoff_m

            output_smstorage_in[t] = self.infiltration_depth_m
            output_smstorage_out[t] = (
                self.flux_lat_m
                + self.flux_perc_m
                + self.actual_et_from_soil_m_per_timestep
                + self.diff_infilt
            )

            output_nashstorage[t] = np.sum(self.nash_storage)
            output_nashstorage_in[t] = self.flux_lat_m
            output_nashstorage_out[t] = self.flux_nash_lateral_runoff_m

        self.cfe_output_data["Time"] = output_time
        self.cfe_output_data["Time Step"] = output_ts
        self.cfe_output_data["Rainfall"] = output_rainfall
        self.cfe_output_data["Direct Runoff"] = output_directrunoff
        self.cfe_output_data["GIUH Runoff"] = output_GIUHrunoff
        self.cfe_output_data["Lateral Flow"] = output_lateralflow
        self.cfe_output_data["Base Flow"] = output_baseflow
        self.cfe_output_data["Total Discharge"] = output_totdischarge
        self.cfe_output_data["Flow"] = output_flow

        self.cfe_output_data["GW storage"] = output_gwstorage
        self.cfe_output_data["Influx to GW storage"] = output_gwstorage_in
        self.cfe_output_data["Outflux from GW storage"] = output_gwstorage_out

        self.cfe_output_data["GIUH storage"] = output_giuhstorage
        self.cfe_output_data["Influx to GIUH storage"] = output_giuhstorage_in
        self.cfe_output_data["Outflux from GIUH storage"] = output_giuhstorage_out

        self.cfe_output_data["SM storage"] = output_SM
        self.cfe_output_data["Influx to SM storage"] = output_smstorage_in
        self.cfe_output_data["Outflux from SM storage"] = output_smstorage_out
        self.cfe_output_data["Soil Moisture Content"] = [
            number / self.D_noahMP for number in output_SM
        ]

        self.cfe_output_data["Nash storage"] = output_nashstorage
        self.cfe_output_data["Influx to Nash storage"] = output_nashstorage_in
        self.cfe_output_data["Outflux from Nash storage"] = output_nashstorage_out

        # Plottings
        if plot:
            self.sync_data()
            self.plot_precip(plot_lims)
            self.plot_flow(plot_lims)
            self.plot_SM(plot_lims)

            for output_type in [
                "GW storage",
                "Nash storage",
                "GIUH storage",
            ]:
                self.plot_other_fluxes(plot_lims, output_type)

        return self.cfe_output_data

    def plot_precip(self, plot_lims):
        plt.plot(
            self.cfe_output_data["Rainfall"][plot_lims],
            label="Precipitation",
            c="gray",
            lw=0.3,
        )
        ax = plt.gca()
        ax.set_title("Precipitation", loc="left", fontweight="bold")
        plt.xlabel("Time [hr]")
        plt.ylabel("[m/hr]")
        plt.legend(loc="upper right")
        plt.show()

    def sync_data(self):
        sim0 = self.cfe_output_data
        obs0 = self.unit_test_data

        sim_synced = pd.DataFrame()
        obs_synced = pd.DataFrame()

        # Get the results
        # Get the simulated data
        sim = sim0[["Time", "Flow", "Soil Moisture Content"]].copy()
        sim["Time"] = pd.to_datetime(
            sim["Time"], format="%Y-%m-%d %H:%M:%S"
        )  # Works specifically for CFE

        # Get the comparison data
        obs = obs0[["Time", "Flow", "Soil Moisture Content"]].copy()
        try:
            obs["Time"] = pd.to_datetime(obs["Time"], format="%m/%d/%Y %H:%M")
        except:  # Works specifically for Mahurangi data
            try:
                obs["Time"] = pd.to_datetime(obs["Time"], format="%d-%m-%Y %H:%M:%S")
            except:
                obs["Time"] = pd.to_datetime(obs["Time"], format="%Y-%m-%d %H:%M:%S")

        # Merge observed and simulated timeseries
        df = pd.merge_asof(sim, obs, on="Time")
        self.df_timeaxis = df["Time"]

        sim_synced["Flow"] = df["Flow" + "_x"].copy()
        obs_synced["Flow"] = df["Flow" + "_y"].copy()
        sim_synced["Soil Moisture Content"] = df["Soil Moisture Content" + "_x"].copy()
        obs_synced["Soil Moisture Content"] = df["Soil Moisture Content" + "_y"].copy()

        self.obs_synced = obs_synced
        self.sim_synced = sim_synced

    def plot_evaluate(self, output_type):
        _obs_synced = self.obs_synced[output_type]
        _sim_synced = self.sim_synced[output_type]

        KGE = spotpy.objectivefunctions.kge(
            _obs_synced[~np.isnan(_obs_synced)],
            _sim_synced[~np.isnan(_obs_synced)],
        )
        NSE = spotpy.objectivefunctions.nashsutcliffe(
            _obs_synced[~np.isnan(_obs_synced)],
            _sim_synced[~np.isnan(_obs_synced)],
        )
        return KGE, NSE

    def plot_flow(self, plot_lims):
        output_type = "Flow"
        KGE, NSE = self.plot_evaluate(output_type)

        for i in range(0, 2):
            plt.plot(
                self.cfe_output_data[output_type][plot_lims],
                label="Simulated "
                + "\n(KGE = "
                + "{:1.4f}".format(KGE)
                + ", NSE= "
                + "{:1.4f}".format(NSE)
                + ")",
            )
            plt.plot(
                self.unit_test_data[output_type][plot_lims], "--", label="Observed"
            )
            ax = plt.gca()
            if (i == 1) or (i == 2):
                plt.yscale("log")
                ax.set_title(
                    output_type + " (log scale)", loc="left", fontweight="bold"
                )
            else:
                ax.set_title(output_type, loc="left", fontweight="bold")
            plt.xlabel("Time [hr]")
            plt.ylabel("[m/hr]")
            plt.legend(loc="upper right")
            plt.show()

    def plot_SM(self, plot_lims):
        output_type = "Soil Moisture Content"
        KGE, NSE = self.plot_evaluate(output_type)

        plt.plot(
            self.cfe_output_data[output_type][plot_lims],
            label="Simulated "
            + "\n(KGE = "
            + "{:1.4f}".format(KGE)
            + ", NSE= "
            + "{:1.4f}".format(NSE)
            + ")",
        )
        plt.plot(self.unit_test_data[output_type][plot_lims], "--", label="Observed")
        ax = plt.gca()
        ax.set_title(output_type, loc="left", fontweight="bold")
        plt.xlabel("Time [hr]")
        ax.set_ylabel("Storage [m3/m3]")
        plt.legend(loc="upper right")
        plt.show()

    def plot_other_fluxes(self, plot_lims, output_type):
        plt.plot(
            self.cfe_output_data[output_type][plot_lims],
            color="dodgerblue",
            linestyle="-",
            label="Storage",
        )

        ax = plt.gca()
        ax2 = ax.twinx()
        ax2.plot(
            self.cfe_output_data["Influx to " + output_type][plot_lims].cumsum(),
            color="seagreen",
            linestyle="--",
            label="Influx",
        )
        ax2.plot(
            self.cfe_output_data["Outflux from " + output_type][plot_lims].cumsum(),
            color="crimson",
            linestyle="--",
            label="Outflux",
        )
        ax.set_title(output_type, loc="left", fontweight="bold")
        plt.xlabel("Time [hr]")
        ax.set_ylabel("Storage [m]")
        ax2.set_ylabel("Cumulative flux [m]")
        plt.legend(loc="upper right")
        plt.show()

    def record_output_fluxes(self, t, output_data):
        output_data[t, "Direct Runoff"] = self.surface_runoff_depth_m
        output_data[t, "GIUH Runoff"] = self.flux_giuh_runoff_m
        output_data[t, "Lateral Flow"] = self.flux_nash_lateral_runoff_m
        output_data[t, "Base Flow"] = self.flux_from_deep_gw_to_chan_m
        output_data[t, "Flow"] = self.flux_Qout_m
        output_data[t, "Total Discharge"] = self.total_discharge
        output_data[t, "Soil Moisture Content"] = (
            self.soil_reservoir["storage_m"] / self.D_noahMP
        )
        output_data[t, "GW storage"] = self.gw_reservoir["storage_m"]
        output_data[t, "Influx to GW storage"] = self.flux_perc_m
        output_data[t, "Outflux from GW storage"] = (
            self.flux_from_deep_gw_to_chan_m + self.diff_perc
        )
        output_data[t, "GIUH storage"] = np.sum(self.runoff_queue_m_per_timestep)
        output_data[t, "Influx to GIUH storage"] = (
            self.surface_runoff_depth_m + self.diff_perc + self.diff_infilt
        )
        output_data[t, "Outflux from GIUH storage"] = self.flux_giuh_runoff_m
        output_data[t, "Influx to Soil Moisture Content"] = self.infiltration_depth_m
        output_data[t, "Outflux from Soil Moisture Content"] = (
            self.flux_lat_m
            + self.flux_perc_m
            + self.actual_et_from_soil_m_per_timestep
            + self.diff_infilt
        )
        output_data[t, "Nash storage"] = np.sum(self.nash_storage)
        output_data[t, "Influx to Nash storage"] = self.flux_lat_m
        output_data[t, "Outflux from Nash storage"] = self.flux_nash_lateral_runoff_m

        return output_data

    def record_forcing(self, t, output_data):
        output_data[t, "Time"] = self.current_time
        output_data[t, "Time Step"] = self.current_time_step
        output_data[t, "Rainfall"] = self.timestep_rainfall_input_m
        return output_data
