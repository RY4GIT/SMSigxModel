import time
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import optimize
import datetime

def sine_func(x, A, phi, k):
    w = 2*np.pi/365
    return A * np.sin(w * x - phi) + k

class SMSig():

    def __init__(self, ts_time, ts_value):

        if sum(ts_value)/len(ts_value) > 1.5:
            # The data is likely to be in percentage. Convert to VSMC
            ts_value = ts_value/100

        # Get daily average of timeseries of data
        dti = pd.to_datetime(ts_time)
        ts = pd.Series(ts_value, index=dti)
        ts_avg = ts.resample('D').mean()

        # read time series of data
        self.time = ts_avg.index.to_numpy()
        self.ts_value = ts_avg.to_numpy()

        # Timestamp in seconds
        ts_timestamp_ns = self.time - np.full((len(self.time),), np.datetime64('1970-01-01T00:00:00Z'))
        ts_timestamp_s = ts_timestamp_ns.astype('timedelta64[D]')
        self.ts_timestamp = ts_timestamp_s.astype('int')

        # print(time, self.ts_value)
        # self.timestep = hourly # TODO: make it flexible later

    def regular(self):
        # make sure the ts_value is in regular interval
        None
        # return detrended ts_value

    def detrend(self):
        # detrend the ts_value using moving average using 1-year window
        windowsize = 365 # 1wk just for now
        halfwindow = int(np.rint(windowsize/2))
        y = np.convolve(self.ts_value, np.ones(windowsize)/(windowsize), mode='same')

        # repeat the first and last observaions to prevent data loss
        y[0:halfwindow] = y[halfwindow]
        y[len(y)-halfwindow:-1] = y[len(y)-halfwindow]
        y[-1] = y[len(y)-halfwindow]

        # Subtract moving average from the original series to detrend the data
        # Additive decomposition
        ts_dtr0 = self.ts_value - y

        # Add 0.5 to avoid negative SM value (do not trust the absolute value of SM anymore
        ts_dtr = ts_dtr0 + 0.5

        self.ts_value = ts_dtr
        # plot
        # print(self.ts_value.shape, self.time.shape)
        # plt.plot(self.ts_value, label='original')
        # plt.plot(y, label='movmean')
        # plt.plot(self.ts_dtr, label='detrended')
        # plt.show()
        # plt.legend()

    def calc_sinecurve(self):
        # https://scipy-lectures.org/intro/scipy/auto_examples/plot_curve_fit.html
        # Get the sine curve information from insitu data
        pseudo_idx = np.arange(start=0, step=1, stop=len(self.ts_value))
        # pseudo_idx = list(range(0, len(self.ts_value), 1))
        params, params_covariance = optimize.curve_fit(sine_func, xdata=self.ts_timestamp, ydata=self.ts_value,
                                                       p0=[1, 0, 0.5])
        phi = params[1]

        # Get the transition valley
        # note that phi sign is opposite from MATLAB code
        sine_n = int(np.round(len(self.ts_timestamp)/365))
        sine_start0 = np.floor(self.ts_timestamp[0]/365 - phi/2/np.pi)
        sine_start_v = 365/2/np.pi * (2*sine_start0*np.pi + np.pi/2 + phi)
        valley = np.arange(start=sine_start_v, step=365, stop = sine_start_v+ 365* (sine_n+1))
        self.t_valley = np.full((len(valley),), np.datetime64('1970-01-01T00:00:00Z')) + np.full((len(valley),), np.timedelta64(1,'D')) * valley

        # Plot and confirm
        plt.figure(figsize=(6, 4))
        plt.scatter(self.ts_timestamp, self.ts_value, label='Data')
        plt.plot(self.ts_timestamp, sine_func(self.ts_timestamp, params[0], params[1], params[2]),
                 label='Fitted function', color='k')
        plt.scatter(valley, sine_func(valley, params[0], params[1], params[2]), color='k')
        plt.legend(loc='best')
        plt.show()

        # return the dates

    def calc_fcwp(self):
        # TODO: Calculate fc and wp to constrain the SM ts_value better
        None

    def calc_seasontrans(self):
        None
        # return signatures

    def compare_results(self, sim, obs):
        None
        # return errors

