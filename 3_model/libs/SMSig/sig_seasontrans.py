import time
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import optimize
import datetime

def sine_func(x, A, phi, k):
    w = 2*np.pi/(365*24)
    return A * np.sin(w * x - phi) + k

class SMSig():

    def __init__(self, time, timeseries):
        # read time series of data
        self.time = time
        # self.time_epoch = self.time.timestamp()
        self.timeseries = timeseries
        if sum(timeseries)/len(timeseries) > 1.5:
            # The data is likely to be in percentage. Convert to VSMC
            self.timeseries = timeseries/100
        print(time, self.timeseries)
        # self.timestep = hourly # TODO: make it flexible later

    def regular(self):
        # make sure the timeseries is in regular interval
        None
        # return detrended timeseries

    def detrend(self):
        # detrend the timeseries using moving average using 1-year window
        windowsize = 24 * 365 # 1wk just for now
        halfwindow = int(np.rint(windowsize/2))
        y = np.convolve(self.timeseries, np.ones(windowsize)/(windowsize), mode='same')

        # repeat the first and last observaions to prevent data loss
        y[0:halfwindow] = y[halfwindow]
        y[len(y)-halfwindow:-1] = y[len(y)-halfwindow]
        y[-1] = y[len(y)-halfwindow]

        # Subtract moving average from the original series to detrend the data
        # Additive decomposition
        ts_dtr0 = self.timeseries - y

        # Add 0.5 to avoid negative SM value (do not trust the absolute value of SM anymore
        ts_dtr = ts_dtr0 + 0.5

        self.timeseries = ts_dtr
        # plot
        # print(self.timeseries.shape, self.time.shape)
        # plt.plot(self.timeseries, label='original')
        # plt.plot(y, label='movmean')
        # plt.plot(self.ts_dtr, label='detrended')
        # plt.show()
        # plt.legend()

    def calc_sinecurve(self):
        # https://scipy-lectures.org/intro/scipy/auto_examples/plot_curve_fit.html
        # Get the sine curve information from insitu data
        pseudo_idx = np.arange(start=0, step=1, stop=len(self.timeseries))
        # pseudo_idx = list(range(0, len(self.timeseries), 1))
        params, params_covariance = optimize.curve_fit(sine_func, xdata=pseudo_idx, ydata=self.timeseries,
                                                       p0=[1, 0, 0.5])
        phi = params[1]
        # Get the transiiton valley
        # sine_n = int(np.floor(len(self.time)/(365*24)))
        # sine_start0 = np.floor(self.time[0]/365 + phi/2/np.pi)
        # sine_start0 = np.floor(self.time[0]/365 + phi/2/np.pi)

        """
        plt.figure(figsize=(6, 4))
        plt.scatter(pseudo_idx, self.timeseries, label='Data')
        plt.plot(pseudo_idx, sine_func(pseudo_idx, params[0], params[1], params[2]),
                 label='Fitted function', color='k')
        plt.legend(loc='best')
        plt.show()
        """
        # return the dates

    def calc_fcwp(self):
        # TODO: Calculate fc and wp to constrain the SM timeseries better
        None

    def calc_seasontrans(self):
        None
        # return signatures

    def compare_results(self, sim, obs):
        None
        # return errors

