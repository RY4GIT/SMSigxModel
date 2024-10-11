import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import datetime
from scipy.optimize import curve_fit
import math


def piecewise_linear(x, P0, P1, P2, P3):
    y0 = np.where(x - P2 > 0, P0 + P1 * x, P0 + P1 * P2)
    return np.where(x - (P2 + P3) > 0, P0 + P1 * (P2 + P3), y0)


def datetime_to_timestamp(t):
    """Convert timestamp in seconds"""
    tstamp_ns = t - np.full((len(t),), np.datetime64("1970-01-01T00:00:00Z"))
    tstamp_d = tstamp_ns.astype("timedelta64[D]")
    return tstamp_d.astype("int")  # Timestamp array in dates


def movmean(sm, windowsize=30, plot=False):
    """Detrend the sm using moving average using 1-year window"""
    # TODO: Make window size flexibl

    # Get the moving average
    halfwindow = int(np.rint(windowsize / 2))
    sm_movmean = np.convolve(sm, np.ones(windowsize) / (windowsize), mode="same")

    # Repeat the first and last observaions to prevent data loss
    sm_movmean[0:halfwindow] = sm_movmean[halfwindow]
    sm_movmean[len(sm_movmean) - halfwindow : -1] = sm_movmean[
        len(sm_movmean) - halfwindow
    ]
    sm_movmean[-1] = sm_movmean[len(sm_movmean) - halfwindow]

    if plot:
        plt.plot(sm, label="original")
        plt.plot(sm_movmean, label="movmean")
        plt.show()
        plt.legend()

    return sm_movmean


class SMSig:
    def __init__(self, t, sm, plot_results, verbose=False):

        # Print/verbose settings
        self.plot_results = plot_results
        self.verbose = verbose

        self.t = t  # datetime
        self.sm = movmean(sm)  # soil moisture values
        self.tt = pd.Series(self.sm, index=self.t)  # timetable

        # Define a vairable
        self.t_date = datetime_to_timestamp(self.t)

        # self.timestep = hourly # TODO: make it flexible later

    def get_dates_to_crop_timeseries(self, transition_type, n_year):

        try:
            subset = self.seasonal_cycle[
                self.seasonal_cycle["transition"] == transition_type
            ].copy()
            start_date = subset["start_date"].values[n_year]
            end_date = subset["end_date"].values[n_year]
        except:
            start_date = None
            end_date = None
        return start_date, end_date

    def crop_timeseries(self, start_date, end_date):
        # Crop the season with 1 month buffer
        mask = (self.tt.index >= start_date) & (self.tt.index <= end_date)
        return self.tt.loc[mask]

    def fit_model_to_seasonaltransition(self, cropped_series, ax=None):
        if (
            cropped_series.empty
            or cropped_series.isnull().mean() > 0.3
            or len(cropped_series) < 90
        ):
            return None  # Indicates that the data is not suitable for analysis

        x = np.arange(1, len(cropped_series) + 1)
        y = cropped_series.to_numpy()
        bounds = self.get_bounds(y)
        p0 = self.get_initial_params(y)

        popt, _ = curve_fit(piecewise_linear, x, y, p0=p0, bounds=bounds)

        if self.verbose:
            print(
                self.current_season_type,
                f"Shift={popt[0]:.2f} slope={popt[1]:.05f} start timing: {popt[2]:02f} end timing: {popt[3]:03f}",
            )

        if self.plot_results:
            self.plot_transition(ax, x, y, popt)

        return popt

    def plot_transition(self, ax, x, y, popt):
        lcolor = "#67a9cf"  # blue

        # Using the piecewise linear function passed as an argument
        ax.plot(x, piecewise_linear(x, *popt), "--", color=lcolor)
        ax.plot(
            x, y, "o", color=lcolor, markersize=3, alpha=0.5
        )  # Plot the data points
        ax.legend()

    def get_bounds(self, y):
        # Define bounds based on season type
        if self.current_season_type == "dry2wet":
            # (y-shift, slope, 1-st inflection point, transition duration)
            lower_bound = (
                0,
                0,
                0,
                self.current_config["bounds"]["transition_duration"][0],
            )
            upper_bound = (
                1,
                1,
                len(y) / 2,
                self.current_config["bounds"]["transition_duration"][1],
            )
            return lower_bound, upper_bound
        elif self.current_season_type == "wet2dry":
            lower_bound = (
                -1,
                -1,
                0,
                self.current_config["bounds"]["transition_duration"][0],
            )
            upper_bound = (
                1,
                0,
                len(y) / 2,
                self.current_config["bounds"]["transition_duration"][1],
            )
            return lower_bound, upper_bound

    def get_initial_params(self, y):
        # Define initial parameters based on season type
        if self.current_season_type == "dry2wet":
            return [
                min(y),
                0.01,
                self.current_config["initial_params"]["P3"],
                self.current_config["initial_params"]["P4"],
            ]
        elif self.current_season_type == "wet2dry":
            return [
                max(y),
                -0.01,
                self.current_config["initial_params"]["P3"],
                self.current_config["initial_params"]["P4"],
            ]

    def calc_seasontrans(self, seasonal_cycle, parameter_config):

        # __________________________________
        # Initialization

        self.seasonal_cycle = (
            seasonal_cycle  # (year - by - (transition_type, start_date, end_date))
        )
        n_years = math.ceil(len(seasonal_cycle) / 2)
        self.config = parameter_config

        # initialization
        season_types = ["dry2wet", "wet2dry"]

        seasontrans_date = np.empty((n_years, 4))
        seasontrans_date[:] = np.nan

        if self.plot_results:
            fig, axs = plt.subplots(
                n_years,
                2,
                figsize=(8, 4 * (n_years)),
            )

        # __________________________________
        ## Main execusion
        # Loop for seasonal types
        for j, season_type in enumerate(season_types):
            self.current_season_type = season_type
            self.current_config = self.config["transitions"][season_type]
            # Loop for water years
            for i in range(n_years):
                if self.verbose:
                    print(
                        f"Processing {self.seasonal_cycle['start_date'][i].year}:{self.current_season_type}"
                    )

                # __________________________________
                # Crop the timeseries
                start_date, end_date = self.get_dates_to_crop_timeseries(
                    self.current_season_type, i
                )
                if start_date is None and end_date is None:
                    # If both are None, then continue with the next iteration or part of the code.
                    print(
                        f"Check seasonal_cycle.csv: No data for {self.current_season_type}, {i}-th year"
                    )
                    continue
                seasonsm_tt = self.crop_timeseries(start_date, end_date)

                # __________________________________
                # Actual signature calculation
                if self.plot_results:
                    ax = axs[i, j]
                else:
                    ax = None

                try:
                    popt = self.fit_model_to_seasonaltransition(seasonsm_tt, ax=ax)

                    # __________________________________
                    # Get signatures in Julian dates
                    trans_start_result = seasonsm_tt.index[0] + datetime.timedelta(
                        days=popt[2]
                    )
                    trans_end_result = seasonsm_tt.index[0] + datetime.timedelta(
                        days=popt[2] + popt[3]
                    )
                    seasontrans_date[i, 2 * j] = trans_start_result.to_julian_date()
                    seasontrans_date[i, 1 + 2 * j] = trans_end_result.to_julian_date()
                except:
                    None

        return seasontrans_date  # Output in julian date
