from eto import ETo
import pandas as pd
import os
import math
import numpy as np

#
# https://eto.readthedocs.io/en/latest/package_references.html
# Only few user exsits, possibly swtich to https://github.com/pyet-org/pyet (only daily) or https://github.com/woodcrafty/PyETo


class FAO_PET:
    def __init__(self, lat=np.nan, lon=np.nan, elevation=0, forcing=None) -> None:
        """
        :param cfg:
        """

        # Get forcing
        self.forcing = forcing
        self.input_forcing = self.prepare_input(forcing)

        # Get CAMELS basin attributes
        self.lon = lon
        self.lat = lat
        self.elevation = elevation

    def calc_PET(self):
        PET = ETo(
            self.input_forcing,
            freq="D",
            lon=self.lon,
            TZ_lon=self.lon,
            z_msl=self.elevation,
            lat=self.lat,
        ).eto_fao()
        # Don't need lat: (only needed if R_s or R_n are not in df)
        # Don't need z_msl ((nly needed if P is not in df)
        # lon: The longitude of the met station (dec deg) (only needed if calculating ETo hourly)
        # TZ_lon: The longitude of the center of the time zone (dec deg) (only needed if calculating ETo hourly).
        # freq (str) – The Pandas time frequency string of the input and output. The minimum frequency is hours (H) and the maximum is month (M).
        PET = PET.fillna(0)  # mm/d? TODO: checkunit

        return PET

    def prepare_input(self, df):
        # Convert the time column to DateTime type
        df["time"] = pd.to_datetime(df["date"])
        df["date"] = pd.to_datetime(df["date"])

        # Set the time column as the index
        df.set_index("time", inplace=True)

        # Calculate relative humidity for each row
        # df["RH_mean"] = df.apply(self.calculate_relative_humidity, axis=1)

        # Actual vapor pressure
        # df["e_a"] = (
        #     df.apply(self.calculate_actual_vapor_pressure_Pa, axis=1)
        #     / self.cfg.conversions.to_kilo
        # )

        # Mean find speed
        # df["U_z"] = (df["wind_u"] + df["wind_v"]) / 2

        # Unit conversion
        # df["R_s"] = (
        #     df["shortwave_radiation"] * 0.0036
        # )  # self.cfg.conversions.day_to_sec / self.cfg.conversions.to_mega
        # df["P"] = df["pressure"] / self.cfg.conversions.to_kilo
        # df["T_mean"] = df["temperature"]

        input_forcing = (
            df  # df[["date", "R_s", "P", "T_mean", "e_a", "RH_mean", "U_z"]]
        )

        return input_forcing

    # References:
    # https://cran.r-project.org/web/packages/humidity/vignettes/humidity-measures.html#eq:1
    # https://www.fao.org/3/x0490e/x0490e07.htm
    # https://github.com/NOAA-OWP/evapotranspiration/blob/1e971ffe784ade3c7ab322ccce6f49f48e942e3d/src/pet.c

    def calculate_relative_humidity(self, row):
        air_sat_vap_press_Pa = self.calc_air_saturation_vapor_pressure_Pa(row)
        actual_vapor_pressure_Pa = self.calculate_actual_vapor_pressure_Pa(row)
        relative_humidity = actual_vapor_pressure_Pa / air_sat_vap_press_Pa
        return relative_humidity

    def calc_air_saturation_vapor_pressure_Pa(self, row):
        air_temperature_C = row["temperature"]
        air_sat_vap_press_Pa = 611.0 * math.exp(
            17.27 * air_temperature_C / (237.3 + air_temperature_C)
        )
        return air_sat_vap_press_Pa

    def calculate_actual_vapor_pressure_Pa(self, row):
        # https://cran.r-project.org/web/packages/humidity/vignettes/humidity-measures.html#eq:1
        q = row["specific_humidity"]  # Specific humidity
        p = row["pressure"]  # Atmospheric pressure in pascals
        actual_vapor_pressure_Pa = q * p / (0.622 + 0.378 * q)
        return actual_vapor_pressure_Pa
