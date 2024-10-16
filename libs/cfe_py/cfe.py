"""
Equation # are from
Ogden, B. F. L. (n.d.). Parameter Estimation for a Conceptual Functional Equivalen (CFE) Formulation of the National Water Model.
"""

import time
import numpy as np
import pandas as pd
from numba import int32, float32  # import the types
from scipy.integrate import solve_ivp, odeint
import sys
from numba import jit

# This function needs to be located outside of the class (global func)
# Soil reservoir ODE for Zone 1
import math
import warnings
import matplotlib.pyplot as plt


@jit(nopython=True)
def conceptual_reservoir_flux_calc(
    t,
    S,
    storage_threshold_primary_m,
    storage_max_m,
    coeff_primary,
    coeff_secondary,
    PET,
    infilt,
    wltsmc_m,
):
    """
    Case 1: S (Soil moisture storage ) > storage_threshold_primary_m
        Interpretation: When the soil moisture is plenty, AET(=PET), percolation, and lateral flow are all active. Thus,
        Equation: dS/dt = Infiltration - PET - (Klf+Kperc) * (S - storage_threshold_primary_m)/(storage_max_m - storage_threshold_primary_m)
    Case 2: storage_threshold_primary_m > S (Soil moisture storage) > storage_threshold_primary_m - wltsmc
        Interpretation: When the soil moisture is in the middle range, AET is active and proportional to the soil moisture storage ratio
        Equation: dS/dt = Infiltration - PET * (S - wltsmc)/(storage_threshold_primary_m - wltsmc)
    Case 3: wltsmc > S (Soil moisture storage)
        Interpretation: When the soil moisture is depleted, no outflux is active
        Equation: dS/dt = Infitlation
    :param t:
    :param S:
    :param storage_threshold_primary_m:
    :param storage_max_m: maximum soil moisture storage, i.e., porosity
    :param coeff_primary: K_perc, percolation coefficient
    :param coeff_secondary: K_lf, lateral flow coefficient
    :param PET: potential evapotranspiration
    :param infilt: infiltration
    :param wltsmc_m: wilting point (in meter)
    :return: dS
    """
    storage_above_threshold_m = S - storage_threshold_primary_m
    storage_diff = storage_max_m - storage_threshold_primary_m
    storage_ratio = np.minimum(storage_above_threshold_m / storage_diff, 1)

    perc_lat_switch = np.multiply(S - storage_threshold_primary_m > 0, 1)
    ## Water-limited but energy-non-limited
    # ET_switch = np.multiply(S > 0, 1)
    ## Energy-limited but water-non-limited
    ET_switch = np.multiply(S - wltsmc_m > 0, 1)

    storage_above_threshold_m_paw = S - wltsmc_m
    storage_diff_paw = storage_threshold_primary_m - wltsmc_m
    ## Water-limited but energy-non-limited
    # storage_ratio_paw = storage_above_threshold_m_paw / storage_diff_paw
    ## Energy-limited but water-non-limited
    storage_ratio_paw = np.minimum(
        storage_above_threshold_m_paw / storage_diff_paw, 1
    )  # Equation 11 (Ogden's document)
    dS = (
        infilt
        - 1 * perc_lat_switch * (coeff_primary + coeff_secondary) * storage_ratio
        - ET_switch * PET * storage_ratio_paw
    )
    return dS


@jit(nopython=True)
def jac(
    t,
    S,
    storage_threshold_primary_m,
    storage_max_m,
    coeff_primary,
    coeff_secondary,
    PET,
    infilt,
    wltsmc_m,
):
    # The Jacobian matrix of the equation conceptual_reservoir_flux_calc. Calculated as (dS/dt)/dS.
    storage_diff = storage_max_m - storage_threshold_primary_m

    perc_lat_switch = np.multiply(S - storage_threshold_primary_m > 0, 1)
    ET_switch = np.multiply(
        (S - wltsmc_m > 0) and (S - storage_threshold_primary_m < 0), 1
    )

    storage_diff_paw = storage_threshold_primary_m - wltsmc_m

    dfdS = (
        -1 * perc_lat_switch * (coeff_primary + coeff_secondary) * 1 / storage_diff
        - ET_switch * PET * 1 / storage_diff_paw
    )
    return [dfdS]


class CFE:
    def __init__(self):
        super(CFE, self).__init__()

    # __________________________________________________________________________________________________________
    # MAIN MODEL FUNCTION
    def run_cfe(self, cfe_state):
        # ________________________________________________
        # Calculate the input rainfall and PET
        cfe_state.volin += cfe_state.timestep_rainfall_input_m
        cfe_state.potential_et_m_per_timestep = (
            cfe_state.potential_et_m_per_s * cfe_state.time_step_size
        )
        cfe_state.reduced_potential_et_m_per_timestep = (
            cfe_state.potential_et_m_per_s * cfe_state.time_step_size
        )
        cfe_state.vol_PET += cfe_state.potential_et_m_per_timestep

        # ________________________________________________
        # Calculate evaporation from rainfall
        cfe_state.actual_et_from_rain_m_per_timestep = 0
        self.et_from_rainfall(cfe_state)

        cfe_state.vol_et_from_rain += cfe_state.actual_et_from_rain_m_per_timestep
        cfe_state.vol_et_to_atm += cfe_state.actual_et_from_rain_m_per_timestep
        cfe_state.volout += cfe_state.actual_et_from_rain_m_per_timestep

        # ________________________________________________
        # Calculate the soil moisture deficit
        cfe_state.soil_reservoir_storage_deficit_m = (
            cfe_state.soil_params["smcmax"]
            * cfe_state.soil_params[
                "D"
            ]  # Ryoko modified here from "cfe_state.soil_params['smcmax'] * cfe_state.soil_params['D']"
            - cfe_state.soil_reservoir["storage_m"]
        )  # Schaake partitioning function 3  (Ogden's document)

        # ________________________________________________
        # Calculates infiltration excess overland flow
        self.Schaake_partitioning_scheme(cfe_state)
        # self.Xinanjiang_partitioning_scheme(cfe_state)
        cfe_state.vol_sch_runoff += cfe_state.surface_runoff_depth_m
        cfe_state.vol_sch_runoff_IOF += cfe_state.surface_runoff_depth_m

        # ________________________________________________
        # Calculates saturation excess overland flow
        # If the infiltration is more than the soil moisture deficit, additional SOF occurs
        if cfe_state.soil_reservoir_storage_deficit_m < cfe_state.infiltration_depth_m:
            self.diff_inf = (
                cfe_state.infiltration_depth_m
                - cfe_state.soil_reservoir_storage_deficit_m
            )
            cfe_state.vol_sch_runoff += self.diff_inf
            cfe_state.surface_runoff_depth_m += self.diff_inf
            cfe_state.vol_sch_runoff_SOF += self.diff_inf
            cfe_state.infiltration_depth_m = cfe_state.soil_reservoir_storage_deficit_m

        # Final infiltration value
        cfe_state.vol_sch_infilt += cfe_state.infiltration_depth_m

        # ________________________________________________
        # Solve ODE for the soil reservoir
        # Calculates primary_flux, secondary_flux, AET, and infiltration simultaneously
        cfe_state.actual_et_from_soil_m_per_timestep = 0
        self.soil_reservoir_flux_calc(cfe_state, cfe_state.soil_reservoir)

        cfe_state.flux_perc_m = cfe_state.primary_flux_m
        cfe_state.flux_lat_m = cfe_state.secondary_flux_m
        cfe_state.vol_et_from_soil += cfe_state.actual_et_from_soil_m_per_timestep
        cfe_state.vol_et_to_atm += cfe_state.actual_et_from_soil_m_per_timestep
        cfe_state.volout += cfe_state.actual_et_from_soil_m_per_timestep

        # ________________________________________________
        # Calculates groundwater storage deficit
        cfe_state.gw_reservoir_storage_deficit_m = (
            cfe_state.gw_reservoir["storage_max_m"]
            - cfe_state.gw_reservoir["storage_m"]
        )

        # ________________________________________________
        # Calculates saturation excess overland flow
        # When the groundwater storage is full, the overflowing amount goes to direct runoff
        if cfe_state.flux_perc_m > cfe_state.gw_reservoir_storage_deficit_m:
            self.diff_perc = (
                cfe_state.flux_perc_m - cfe_state.gw_reservoir_storage_deficit_m
            )
            cfe_state.flux_perc_m = cfe_state.gw_reservoir_storage_deficit_m
            cfe_state.vol_sch_runoff += self.diff_perc
            cfe_state.surface_runoff_depth_m += self.diff_perc
            cfe_state.vol_sch_runoff_SOF += self.diff_perc
            cfe_state.vol_sch_infilt -= self.diff_perc

        # Finalize the percolation and lateral flow
        cfe_state.gw_reservoir["storage_m"] += cfe_state.flux_perc_m
        cfe_state.vol_to_gw += cfe_state.flux_perc_m
        cfe_state.vol_soil_to_gw += cfe_state.flux_perc_m
        cfe_state.vol_soil_to_lat_flow += (
            cfe_state.flux_lat_m
        )  # TODO add this to nash cascade as input
        cfe_state.volout += cfe_state.flux_lat_m

        # ________________________________________________
        # Solve groundwater reservoir

        self.groundwater_reservoir_flux_calc(cfe_state, cfe_state.gw_reservoir)
        cfe_state.flux_from_deep_gw_to_chan_m = np.clip(
            cfe_state.gw_primary_flux_m, 0, cfe_state.gw_reservoir["storage_m"]
        )
        cfe_state.gw_reservoir["storage_m"] -= cfe_state.flux_from_deep_gw_to_chan_m
        # if cfe_state.gw_reservoir["storage_m"] < 0:
        #     print(cfe_state.gw_reservoir["storage_m"])
        cfe_state.vol_from_gw += cfe_state.flux_from_deep_gw_to_chan_m
        cfe_state.volout += cfe_state.flux_from_deep_gw_to_chan_m

        # ________________________________________________
        # SUBROUTINE
        # giuh_runoff_m = f(Schaake_output, giuh_ordinates, runoff_queue_m_per_timestep)
        self.convolution_integral(cfe_state)
        cfe_state.vol_out_giuh += cfe_state.flux_giuh_runoff_m
        cfe_state.volout += cfe_state.flux_giuh_runoff_m

        # ________________________________________________
        # SUBROUTINE
        self.nash_cascade(cfe_state)
        cfe_state.vol_in_nash += cfe_state.flux_lat_m
        cfe_state.vol_out_nash += cfe_state.flux_nash_lateral_runoff_m
        # (?) missing to add vol_nash to vol_out

        # ________________________________________________
        cfe_state.flux_Qout_m = (
            cfe_state.flux_giuh_runoff_m
            + cfe_state.flux_nash_lateral_runoff_m
            + cfe_state.flux_from_deep_gw_to_chan_m
        )
        cfe_state.total_discharge = (
            cfe_state.flux_Qout_m
            * cfe_state.catchment_area_km2
            * 1000000.0
            / cfe_state.time_step_size
        )

        # ________________________________________________
        cfe_state.current_time_step += 1
        cfe_state.current_time += pd.Timedelta(value=cfe_state.time_step_size, unit="s")

        return

    # __________________________________________________________________________________________________________
    def nash_cascade(self, cfe_state):
        """
        Solve for the flow through the Nash cascade to delay the
        arrival of the lateral flow into the channel
        """
        Q = np.zeros(cfe_state.num_lateral_flow_nash_reservoirs)

        for i in range(cfe_state.num_lateral_flow_nash_reservoirs):
            Q[i] = cfe_state.K_nash * cfe_state.nash_storage[i]

            cfe_state.nash_storage[i] -= Q[i]

            if i == 0:
                cfe_state.nash_storage[i] += cfe_state.flux_lat_m

            else:
                cfe_state.nash_storage[i] += Q[i - 1]

        cfe_state.flux_nash_lateral_runoff_m = Q[
            cfe_state.num_lateral_flow_nash_reservoirs - 1
        ]
        # print(cfe_state.nash_storage)
        return

    # __________________________________________________________________________________________________________
    def convolution_integral(self, cfe_state):
        """
        This function solves the convolution integral involving N GIUH ordinates.

        Inputs:
            Schaake_output_runoff_m
            num_giuh_ordinates
            giuh_ordinates
        Outputs:
            runoff_queue_m_per_timestep
        """

        #        cfe_state.runoff_queue_m_per_timestep[-1] = 0

        for i in range(cfe_state.num_giuh_ordinates):
            cfe_state.runoff_queue_m_per_timestep[i] += (
                cfe_state.giuh_ordinates[i] * cfe_state.surface_runoff_depth_m
            )

        cfe_state.flux_giuh_runoff_m = cfe_state.runoff_queue_m_per_timestep[0]

        # __________________________________________________________________
        # shift all the entries in preperation for the next timestep

        for i in range(1, cfe_state.num_giuh_ordinates):
            cfe_state.runoff_queue_m_per_timestep[i - 1] = (
                cfe_state.runoff_queue_m_per_timestep[i]
            )

        cfe_state.runoff_queue_m_per_timestep[-1] = 0

        return

    # __________________________________________________________________________________________________________
    def et_from_rainfall(self, cfe_state):
        """
        iff it is raining, take PET from rainfall first.  Wet veg. is efficient evaporator.
        """

        if cfe_state.timestep_rainfall_input_m > 0.0:
            if cfe_state.verbose:
                print("Raining")
            if (
                cfe_state.timestep_rainfall_input_m
                > cfe_state.potential_et_m_per_timestep
            ):
                cfe_state.actual_et_from_rain_m_per_timestep = (
                    cfe_state.potential_et_m_per_timestep
                )
                cfe_state.timestep_rainfall_input_m -= (
                    cfe_state.actual_et_from_rain_m_per_timestep
                )

            else:
                # cfe_state.actual_et_from_rain_m_per_timestep = cfe_state.potential_et_m_per_timestep
                cfe_state.actual_et_from_rain_m_per_timestep = (
                    cfe_state.timestep_rainfall_input_m
                )
                cfe_state.timestep_rainfall_input_m = 0.0

        cfe_state.reduced_potential_et_m_per_timestep = (
            cfe_state.potential_et_m_per_timestep
            - cfe_state.actual_et_from_rain_m_per_timestep
        )

        # else:
        # print('not raining')
        # print(cfe_state.reduced_potential_et_m_per_timestep)

        return

    def soil_reservoir_flux_calc(self, cfe_state, reservoir):
        """
        This function solves the soil moisture mass balance.

        Inputs:
            reservoir
        Outputs:
            primary_flux_m (percolation)
            secondary_flux_m (lateral flow)
            actual_et_from_soil_m_per_timestep (et_from_soil)
        """

        # Initialization
        y0 = [reservoir["storage_m"]]
        if cfe_state.time_step_size == 86400:
            t = np.linspace(0, 1, 24 * 5)
        elif cfe_state.time_step_size == 3600:
            t = np.array(0, 1, 5)

        # Solve and ODE
        sol = odeint(
            conceptual_reservoir_flux_calc,
            y0,
            t,
            args=(
                reservoir["storage_threshold_primary_m"],
                reservoir["storage_max_m"],
                reservoir["coeff_primary"],
                reservoir["coeff_secondary"],
                cfe_state.reduced_potential_et_m_per_timestep,
                cfe_state.infiltration_depth_m,
                cfe_state.soil_params["wltsmc"]
                * cfe_state.soil_params["D"],  # wilting point in meter
            ),
            tfirst=True,
            Dfun=jac,
            rtol=None,
            atol=None,
        )

        # Finalize results
        ts_concat = t
        ys_concat = np.concatenate(sol, axis=0)

        # Calculate fluxes
        t_proportion = np.diff(ts_concat)
        ys_avg = np.convolve(ys_concat, np.ones(2), "valid") / 2

        lateral_flux = np.zeros(ys_avg.shape)
        perc_lat_switch = ys_avg - reservoir["storage_threshold_primary_m"] > 0
        lateral_flux = reservoir["coeff_secondary"] * np.clip(
            (ys_avg - reservoir["storage_threshold_primary_m"])
            / (reservoir["storage_max_m"] - reservoir["storage_threshold_primary_m"]),
            0,
            1,
        )
        lateral_flux_frac = lateral_flux * t_proportion

        perc_flux = np.zeros(ys_avg.shape)
        perc_flux = reservoir["coeff_primary"] * np.clip(
            (ys_avg - reservoir["storage_threshold_primary_m"])
            / (reservoir["storage_max_m"] - reservoir["storage_threshold_primary_m"]),
            0,
            1,
        )
        perc_flux_frac = perc_flux * t_proportion

        et_from_soil = np.zeros(ys_avg.shape)
        et_from_soil = cfe_state.reduced_potential_et_m_per_timestep * np.clip(
            (ys_avg - cfe_state.soil_params["wltsmc"] * cfe_state.soil_params["D"])
            / (
                reservoir["storage_threshold_primary_m"]
                - cfe_state.soil_params["wltsmc"] * cfe_state.soil_params["D"]
            ),
            0,
            1,
        )
        et_from_soil_frac = et_from_soil * t_proportion

        infilt_to_soil = np.repeat(cfe_state.infiltration_depth_m, ys_avg.shape)
        infilt_to_soil_frac = infilt_to_soil * t_proportion

        # Scale fluxes
        sum_outflux = lateral_flux_frac + perc_flux_frac + et_from_soil_frac
        if sum_outflux.any() == 0:
            flux_scale = 0
            if cfe_state.infiltration_depth_m > 0:
                # To account for mass balance error by ODE
                final_storage_m = y0[0] + cfe_state.infiltration_depth_m
            else:
                final_storage_m = y0[0]
        else:
            flux_scale = (
                (ys_concat[0] - ys_concat[-1]) + np.sum(infilt_to_soil_frac)
            ) / np.sum(sum_outflux)
            final_storage_m = ys_concat[-1]

        scaled_lateral_flux = lateral_flux_frac * flux_scale
        scaled_perc_flux = perc_flux_frac * flux_scale
        scaled_et_flux = et_from_soil_frac * flux_scale

        # Pass the results
        cfe_state.primary_flux_m = math.fsum(scaled_perc_flux)
        cfe_state.secondary_flux_m = math.fsum(scaled_lateral_flux)
        cfe_state.actual_et_from_soil_m_per_timestep = math.fsum(scaled_et_flux)
        reservoir["storage_m"] = final_storage_m

        """
        # Comment out because this section raises Runtime error, as dS_soil_reservoir is extremely small
        # dS based on Soil reservoir
        dS_soil_reservoir = ys_concat[-1]-ys_concat[0]
        # dS based on fluxes
        dS_fluxes = cfe_state.infiltration_depth_m - cfe_state.primary_flux_m - cfe_state.secondary_flux_m - cfe_state.actual_et_from_soil_m_per_timestep
        if ((dS_soil_reservoir - dS_fluxes) / dS_soil_reservoir) >= 0.01:
            warnings.warn(f'Mass balance error is more than 1%. \n dS({ys_concat[-1]-ys_concat[0]}) = I({cfe_state.infiltration_depth_m}) - Perc({cfe_state.primary_flux_m}) - Lat({cfe_state.secondary_flux_m}) - AET({cfe_state.actual_et_from_soil_m_per_timestep})')
        """

        # __________________________________________________________________________________________________________

    ########## SINGLE OUTLET EXPONENTIAL RESERVOIR ###############
    ##########                -or-                 ###############
    ##########    TWO OUTLET NONLINEAR RESERVOIR   ###############
    def groundwater_reservoir_flux_calc(self, cfe_state, reservoir):
        """
        This function calculates the flux from a linear, or nonlinear
        conceptual reservoir with one or two outlets, or from an
        exponential nonlinear conceptual reservoir with only one outlet.
        In the non-exponential instance, each outlet can have its own
        activation storage threshold.  Flow from the second outlet is
        turned off by setting the discharge coeff. to 0.0.
        """

        # exponential nonliner conceptual reservoir with only one outlet
        if reservoir["is_exponential"] == True:
            flux_exponential = (
                np.exp(  # Equation 12 (Ogden's document).
                    reservoir["exponent_primary"]
                    * reservoir["storage_m"]
                    / reservoir["storage_max_m"]
                )
                - 1
            )  # NWM do not subtracts 1 from the formulation
            cfe_state.gw_primary_flux_m = reservoir["coeff_primary"] * flux_exponential
            cfe_state.secondary_flux_m = 0.0
            return

            # else:
            #     # linear/nonlinear conceptual reservoir with one/two outlets

            #     cfe_state.primary_flux_m = 0.0

            #     storage_above_threshold_m = (
            #         reservoir["storage_m"] - reservoir["storage_threshold_primary_m"]
            #     )  # Equation 11 (Ogden's document).
            #     # print('storage above threshold: %s' % (storage_above_threshold_m))
            #     if storage_above_threshold_m > 0.0:
            #         storage_diff = (
            #             reservoir["storage_max_m"]
            #             - reservoir["storage_threshold_primary_m"]
            #         )  # Equation 11 (Ogden's document).
            #         storage_ratio = (
            #             storage_above_threshold_m / storage_diff
            #         )  # Equation 11 (Ogden's document).
            #         storage_power = np.power(storage_ratio, reservoir["exponent_primary"])

            #         cfe_state.primary_flux_m = reservoir["coeff_primary"] * storage_power

            #         if cfe_state.primary_flux_m > storage_above_threshold_m:
            #             cfe_state.primary_flux_m = storage_above_threshold_m

            #     cfe_state.secondary_flux_m = 0.0

            #     storage_above_threshold_m = (
            #         reservoir["storage_m"] - reservoir["storage_threshold_secondary_m"]
            #     )

            #     if storage_above_threshold_m > 0.0:
            #         storage_diff = (
            #             reservoir["storage_max_m"]
            #             - reservoir["storage_threshold_secondary_m"]
            #         )  # Equation 12 (Ogden's document).
            #         storage_ratio = storage_above_threshold_m / storage_diff
            #         storage_power = np.power(storage_ratio, reservoir["exponent_secondary"])

            #         cfe_state.secondary_flux_m = (
            #             reservoir["coeff_secondary"] * storage_power
            #         )
            #         if cfe_state.secondary_flux_m > (
            #             storage_above_threshold_m - cfe_state.primary_flux_m
            #         ):
            #             cfe_state.secondary_flux_m = (
            #                 storage_above_threshold_m - cfe_state.primary_flux_m
            #             )
            #             # print('all excess water went to primary flux')

            return

    def Xinanjiang_partitioning_scheme(self, cfe_state):
        """
        This module takes the water_input_depth_m and separates it into surface_runoff_depth_m
        and infiltration_depth_m by calculating the saturated area and runoff based on a scheme developed
        for the Xinanjiang model by Jaywardena and Zhou (2000). According to Knoben et al.
        (2019) "the model uses a variable contributing area to simulate runoff.  [It] uses
        a double parabolic curve to simulate tension water capacities within the catchment,
        instead of the original single parabolic curve" which is also used as the standard
        VIC fomulation.  This runoff scheme was selected for implementation into NWM v3.0.
        REFERENCES:
        1. Jaywardena, A.W. and M.C. Zhou, 2000. A modified spatial soil moisture storage
            capacity distribution curve for the Xinanjiang model. Journal of Hydrology 227: 93-113
        2. Knoben, W.J.M. et al., 2019. Supplement of Modular Assessment of Rainfall-Runoff Models
            Toolbox (MARRMoT) v1.2: an open-source, extendable framework providing implementations
            of 46 conceptual hydrologic models as continuous state-space formulations. Supplement of
            Geosci. Model Dev. 12: 2463-2480.
            https://github.com/wknoben/MARRMoT/blob/master/MARRMoT/Models/Model%20files/m_28_xinanjiang_12p_4s.m#L34
        -------------------------------------------------------------------------
        Written by RLM May 2021
        Adapted by JMFrame September 2021 for new version of CFE
        Further adapted by QiyueL August 2022 for python version of CFE
        ------------------------------------------------------------------------
        Inputs
        double  time_step_rainfall_input_m           amount of water input to soil surface this time step [m]
        double  field_capacity_m                     amount of water stored in soil reservoir when at field capacity [m]
        double  max_soil_moisture_storage_m          total storage of the soil moisture reservoir (porosity*soil thickness) [m]
        double  column_total_soil_water_m     current storage of the soil moisture reservoir [m]
        double  a_inflection_point_parameter  a parameter
        double  b_shape_parameter             b parameter
        double  x_shape_parameter             x parameter
            //
        Outputs
        double  surface_runoff_depth_m        amount of water partitioned to surface water this time step [m]
        double  infiltration_depth_m          amount of water partitioned as infiltration (soil water input) this time step [m]
        -------------------------------------------------------------------------
        """

        # partition the total soil water in the column between free water and tension water
        free_water_m = (
            cfe_state.soil_reservoir_storage_deficit_m
        )  # cfe_state.soil_reservoir['storage_m']- cfe_state.soil_reservoir['storage_threshold_primary_m'];

        if 0.0 < free_water_m:
            tension_water_m = cfe_state.soil_reservoir["storage_threshold_primary_m"]

        else:
            free_water_m = 0.0
            tension_water_m = cfe_state.soil_reservoir["storage_m"]

        # estimate the maximum free water and tension water available in the soil column
        max_free_water_m = (
            cfe_state.soil_reservoir["storage_max_m"]
            - cfe_state.soil_reservoir["storage_threshold_primary_m"]
        )
        max_tension_water_m = cfe_state.soil_reservoir["storage_threshold_primary_m"]

        # check that the free_water_m and tension_water_m do not exceed the maximum and if so, change to the max value
        if max_free_water_m < free_water_m:
            free_water_m = max_free_water_m

        if max_tension_water_m < tension_water_m:
            tension_water_m = max_tension_water_m

        """
            NOTE: the impervious surface runoff assumptions due to frozen soil used in NWM 3.0 have not been included.
            We are assuming an impervious area due to frozen soils equal to 0 (see eq. 309 from Knoben et al).

            The total (pervious) runoff is first estimated before partitioning into surface and subsurface components.
            See Knoben et al eq 310 for total runoff and eqs 313-315 for partitioning between surface and subsurface
            components.

            Calculate total estimated pervious runoff. 
            NOTE: If the impervious surface runoff due to frozen soils is added,
            the pervious_runoff_m equation will need to be adjusted by the fraction of pervious area.
        """

        a_Xinanjiang_inflection_point_parameter = -0.49
        b_Xinanjiang_shape_parameter = 10
        x_Xinanjiang_shape_parameter = 10

        # When soil is dry
        if (tension_water_m / max_tension_water_m) <= (
            0.5 - a_Xinanjiang_inflection_point_parameter
        ):
            pervious_runoff_m = cfe_state.timestep_rainfall_input_m * (
                np.power(
                    (0.5 - a_Xinanjiang_inflection_point_parameter),
                    (1.0 - b_Xinanjiang_shape_parameter),
                )
                * np.power(
                    (1.0 - (tension_water_m / max_tension_water_m)),
                    b_Xinanjiang_shape_parameter,
                )
            )

        # When soil is wet
        else:
            pervious_runoff_m = cfe_state.timestep_rainfall_input_m * (
                1.0
                - np.power(
                    (0.5 + a_Xinanjiang_inflection_point_parameter),
                    (1.0 - b_Xinanjiang_shape_parameter),
                )
                * np.power(
                    (1.0 - (tension_water_m / max_tension_water_m)),
                    (b_Xinanjiang_shape_parameter),
                )
            )

        # Separate the surface water from the pervious runoff
        ## NOTE: If impervious runoff is added to this subroutine, impervious runoff should be added to
        ## the surface_runoff_depth_m.

        cfe_state.surface_runoff_depth_m = pervious_runoff_m * (
            1.0
            - np.power(
                (1.0 - (free_water_m / max_free_water_m)), x_Xinanjiang_shape_parameter
            )
        )

        # The surface runoff depth is bounded by a minimum of 0 and a maximum of the water input depth.
        # Check that the estimated surface runoff is not less than 0.0 and if so, change the value to 0.0.
        if cfe_state.surface_runoff_depth_m < 0.0:
            cfe_state.surface_runoff_depth_m = 0.0

        # Check that the estimated surface runoff does not exceed the amount of water input to the soil surface.  If it does,
        # change the surface water runoff value to the water input depth.
        if cfe_state.surface_runoff_depth_m > cfe_state.timestep_rainfall_input_m:
            cfe_state.surface_runoff_depth_m = cfe_state.timestep_rainfall_input_m

        # Separate the infiltration from the total water input depth to the soil surface.
        cfe_state.infiltration_depth_m = (
            cfe_state.timestep_rainfall_input_m - cfe_state.surface_runoff_depth_m
        )

        if cfe_state.timestep_rainfall_input_m > 0:
            if cfe_state.verbose:
                print(
                    f"Soil storage is {(tension_water_m / max_tension_water_m)*100:.2f} % saturated"
                )
                print(
                    f"Runoff is {cfe_state.surface_runoff_depth_m/cfe_state.timestep_rainfall_input_m*100:.2f} % of P (={cfe_state.timestep_rainfall_input_m*1000:.2f}[mm/hr])"
                )

        return

    # __________________________________________________________________________________________________________
    #  SCHAAKE RUNOFF PARTITIONING SCHEME
    def Schaake_partitioning_scheme(self, cfe_state):
        """
        This subtroutine takes water_input_depth_m and partitions it into surface_runoff_depth_m and
        infiltration_depth_m using the scheme from Schaake et al. 1996.
        !--------------------------------------------------------------------------------
        modified by FLO April 2020 to eliminate reference to ice processes,
        and to de-obfuscate and use descriptive and dimensionally consistent variable names.

        inputs:
          timestep_d
          Schaake_adjusted_magic_constant_by_soil_type = C*Ks(soiltype)/Ks_ref, where C=3, and Ks_ref=2.0E-06 m/s
          column_total_soil_moisture_deficit_m (soil_reservoir_storage_deficit_m)
          water_input_depth_m (timestep_rainfall_input_m) amount of water input to soil surface this time step [m]
        outputs:
          surface_runoff_depth_m      amount of water partitioned to surface water this time step [m]
          infiltration_depth_m
        """

        # If the rainfall is happening
        if cfe_state.timestep_rainfall_input_m > 0:
            # If soil water deficit is negative (soil storage is full)
            if 0 > cfe_state.soil_reservoir_storage_deficit_m:
                # All rainfall input goes to runoff, and there is no infiltration
                cfe_state.surface_runoff_depth_m = cfe_state.timestep_rainfall_input_m
                cfe_state.infiltration_depth_m = 0.0
                if cfe_state.verbose:
                    print("SM storage is full")
            # If there is soil water deficit, infiltration is calculated
            else:
                if cfe_state.verbose:
                    print("SM storage is not full")
                schaake_exp_term = np.exp(
                    -cfe_state.Schaake_adjusted_magic_constant_by_soil_type
                    * cfe_state.timestep_d
                )

                Schaake_parenthetical_term = 1.0 - schaake_exp_term

                Ic = (
                    cfe_state.soil_reservoir_storage_deficit_m
                    * Schaake_parenthetical_term
                )

                Px = cfe_state.timestep_rainfall_input_m

                cfe_state.infiltration_depth_m = Px * (Ic / (Px + Ic))

                if cfe_state.verbose:
                    print(f"Infiltration percentage: {Ic / (Px + Ic)*100:.3f} (%)")

                # If the rainfall input exceeds infiltration capacity, the remainings go to runoff
                if 0.0 < (
                    cfe_state.timestep_rainfall_input_m - cfe_state.infiltration_depth_m
                ):
                    cfe_state.surface_runoff_depth_m = (
                        cfe_state.timestep_rainfall_input_m
                        - cfe_state.infiltration_depth_m
                    )

                else:
                    # If the rainfall is less than infiltration capacity, there is no runoff, and infiltration equals to rainfall
                    cfe_state.surface_runoff_depth_m = 0.0
                    cfe_state.infiltration_depth_m = cfe_state.timestep_rainfall_input_m

            if cfe_state.verbose:
                print(
                    f"Runoff is {cfe_state.surface_runoff_depth_m/cfe_state.timestep_rainfall_input_m*100:.2f} % of P (={cfe_state.timestep_rainfall_input_m*1000:.2f}[mm/hr])"
                )

        # If there is no rainfall, both infiltration and runoff are zero
        else:
            cfe_state.surface_runoff_depth_m = 0.0
            cfe_state.infiltration_depth_m = 0.0

        return

    # __________________________________________________________________________________________________________
    def is_fabs_less_than_epsilon(self, a, epsilon):
        if np.abs(a) < epsilon:
            return True

        else:
            return False
