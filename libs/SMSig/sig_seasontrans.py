import time
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import optimize
import datetime
from scipy.optimize import curve_fit
import os


def piecewise_linear(x, P0, P1, P2, P3):
    y0 = np.where(x - P2 > 0, P0 + P1 * x, P0 + P1 * P2)
    return np.where(x - (P2 + P3) > 0, P0 + P1 * (P2 + P3), y0)


def sine_func(x, A, phi, k):
    w = 2 * np.pi / 365
    return A * np.sin(w * x - phi) + k


def datetime_to_timestamp(t):
    # Convert timestamp in seconds
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
        # Get the base start/end days
        if transition_type == "dry2wet":
            trans_start0 = self.t_valley[n_year]
            trans_end0 = trans_start0 + datetime.timedelta(days=300)
            buffer_start = 15
            buffer_end = 0
        elif transition_type == "wet2dry":
            trans_start0 = self.t_valley[n_year] + datetime.timedelta(days=365 / 2)
            trans_end0 = self.t_valley[n_year + 1] + datetime.timedelta(days=0)
            buffer_start = 60
            buffer_end = 15
        return trans_start0, trans_end0, buffer_start, buffer_end

    def crop_timeseries(self, start_date, end_date, buffer_start, buffer_end):
        # Crop the season with 1 month buffer
        mask = (self.tt.index >= start_date - datetime.timedelta(days=buffer_start)) & (
            self.tt.index <= end_date + datetime.timedelta(days=buffer_end)
        )
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
            return (min(y) - 1, 0, 0, 15), (1, 3.0, len(y), len(y) + 60)
        elif self.current_season_type == "wet2dry":
            return (-1, -0.1, 0, 30), (1, 0, len(y), len(y) + 60)

    def get_initial_params(self, y):
        # Define initial parameters based on season type
        if self.current_season_type == "dry2wet":
            return [min(y), 0.1, 50, len(y) - 30]
        elif self.current_season_type == "wet2dry":
            return [max(y), -0.0001, 50, len(y) - 30]

    # # Example method using the configuration
    # def get_dates_to_crop_timeseries(self, transition_type, n_year):

    #     trans_start0 = self.t_valley[n_year]
    #     if transition_type == "wet2dry":
    #         trans_start0 += datetime.timedelta(
    #             days=transition_config["trans_duration_days"]
    #         )
    #     trans_end0 = trans_start0 + datetime.timedelta(
    #         days=transition_config["trans_duration_days"]
    #     )
    #     buffer_start = transition_config["buffer_start_days"]
    #     buffer_end = transition_config["buffer_end_days"]

    #     return trans_start0, trans_end0, buffer_start, buffer_end

    # def get_bounds(self, y):
    #     bounds_config = self.config["transitions"][self.current_season_type]["bounds"]
    #     lower = [
    #         x if isinstance(x, int) else eval(x.replace("len(y)", str(len(y))))
    #         for x in bounds_config["lower"]
    #     ]
    #     upper = [
    #         x if isinstance(x, int) else eval(x.replace("len(y)", str(len(y))))
    #         for x in bounds_config["upper"]
    #     ]
    #     return (lower, upper)

    # def get_initial_params(self, y):
    #     params_config = self.config["transitions"][self.current_season_type]
    #     initial_params = params_config["initial_params"]
    #     if self.current_season_type == "dry2wet":
    #         initial_params.insert(0, min(y))
    #     elif self.current_season_type == "wet2dry":
    #         initial_params.insert(0, max(y))
    #     # Adjust the last two parameters based on y's length
    #     initial_params[2] += len(y)
    #     initial_params[3] += len(y)
    #     return initial_params

    def calc_seasontrans(self, t_valley, parameter_config):
        self.t_valley = t_valley
        self.config = parameter_config

        # initialization
        season_types = ["dry2wet", "wet2dry"]

        seasontrans_date = np.empty((len(self.t_valley), 4))
        seasontrans_date[:] = np.nan

        if self.plot_results:
            fig, axs = plt.subplots(
                len(self.t_valley) - 1, 2, figsize=(8, 4 * (len(self.t_valley) - 1))
            )

        ## Main execusion
        for j, season_type in enumerate(season_types):
            self.current_season_type = season_type
            # Loop for water years
            for i in range(len(self.t_valley) - 1):
                if self.verbose:
                    print(f"Processing {t_valley[i].year}:{self.current_season_type}")

                # __________________________________
                # Crop the timeseries
                _start_date, _end_date, buffer_start, buffer_end = (
                    self.get_dates_to_crop_timeseries(self.current_season_type, i)
                )

                seasonsm_tt = self.crop_timeseries(
                    _start_date, _end_date, buffer_start, buffer_end
                )

                # __________________________________
                # Actual signature calculation
                if self.plot_results:
                    ax = axs[i, j]
                else:
                    ax = None
                popt = self.fit_model_to_seasonaltransition(seasonsm_tt, ax=ax)

                # __________________________________
                # Get signatures in dates

                trans_start_result = seasonsm_tt.index[0] + datetime.timedelta(
                    days=popt[2]
                )
                trans_end_result = seasonsm_tt.index[0] + datetime.timedelta(
                    days=popt[2] + popt[3]
                )

                seasontrans_date[i, 2 * j] = trans_start_result.to_julian_date()
                seasontrans_date[i, 1 + 2 * j] = trans_end_result.to_julian_date()

        return seasontrans_date  # Output in julian date

    def calc_sinecurve(self):
        # https://scipy-lectures.org/intro/scipy/auto_examples/plot_curve_fit.html
        # Get the sine curve information from insitu data
        params, _ = optimize.curve_fit(
            sine_func, xdata=self.t_date, ydata=self.sm, p0=[1, 0, 0.5]
        )
        phi = params[1]

        # Get the transition valley
        # note that phi sign is opposite from MATLAB code
        sine_n = int(np.round(len(self.t_date) / 365))
        sine_start0 = np.floor(self.t_date[0] / 365 - phi / 2 / np.pi)
        sine_start_v = 365 / 2 / np.pi * (2 * sine_start0 * np.pi + np.pi / 2 + phi)
        valley = np.arange(
            start=sine_start_v, step=365, stop=sine_start_v + 365 * (sine_n + 1)
        )
        t_valley = (
            np.full((len(valley),), np.datetime64("1970-01-01T00:00:00Z"))
            + np.full((len(valley),), np.timedelta64(1, "D")) * valley
        )
        t_valley = pd.Series(t_valley)

        # Plot and confirm
        if self.plot_results:
            plt.figure(figsize=(6, 4))
            plt.scatter(self.t_date, self.sm, label="Data")
            plt.plot(
                self.t_date,
                sine_func(self.t_date, params[0], params[1], params[2]),
                label="Fitted function",
                color="k",
            )
            plt.scatter(
                valley, sine_func(valley, params[0], params[1], params[2]), color="k"
            )
            plt.legend(loc="best")
            plt.show()

        # return the dates
        return t_valley

    def detrend(self):
        # detrend the sm using moving average using 1-year window
        windowsize = 365  # 1wk just for now
        halfwindow = int(np.rint(windowsize / 2))
        y = np.convolve(self.sm, np.ones(windowsize) / (windowsize), mode="same")

        # repeat the first and last observations to prevent data loss
        y[0:halfwindow] = y[halfwindow]
        y[len(y) - halfwindow : -1] = y[len(y) - halfwindow]
        y[-1] = y[len(y) - halfwindow]

        # Subtract moving average from the original series to detrend the data
        # Additive decomposition
        ts_dtr0 = self.sm - y

        # Add 0.5 to avoid negative SM value (do not trust the absolute value of SM anymore
        ts_dtr = ts_dtr0 + 0.5

        self.sm = ts_dtr
