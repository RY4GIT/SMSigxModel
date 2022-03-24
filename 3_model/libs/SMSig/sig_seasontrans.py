import time
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import optimize
import datetime
from scipy.optimize import curve_fit

def piecewise_linear(x, P0, P1, P2, P3):
    y0 = np.where(x-P2>0, P0+P1*x, P0+P1*P2)
    return np.where(x-(P2+P3)>0, P0+P1*(P2+P3), y0)

def sine_func(x, A, phi, k):
    w = 2*np.pi/365
    return A * np.sin(w * x - phi) + k

def datetime_to_timestamp(ts_datetime):
    ts_timestamp_ns = ts_datetime - np.full((len(ts_datetime),),
                                                      np.datetime64('1970-01-01T00:00:00Z'))
    ts_timestamp_d = ts_timestamp_ns.astype('timedelta64[D]')
    return ts_timestamp_d

    # Define a vairable
    self.ts_time_d = ts_timestamp_d.astype('int')  # Timestamp array in seconds


class MyTakeStep(object):
   def __init__(self, stepsize=0.5):
       self.stepsize = stepsize
   def __call__(self, x):
       s = self.stepsize
       x[:] += np.random.uniform(-s, s, x[:].shape)
       x[2] += np.random.uniform(-5.*s, 5.*s)
       x[3] += np.random.uniform(-5.*s, 5.*s)
       return x

class SMSig():

    def __init__(self, ts_time, ts_value, plot_results):

        self.plot_results = plot_results

        # If the data is likely to be in percentage. Convert to VSMC
        if sum(ts_value)/len(ts_value) > 1.5:
            ts_value = ts_value/100

        # Aggregate the timeseries of data into daily
        dti = pd.to_datetime(ts_time)
        tt = pd.Series(ts_value, index=dti)
        tt_avg = tt.resample('D').mean()

        # Define variables
        self.tt                 = tt_avg # Timetable
        self.ts_time_datetime   = tt_avg.index.to_numpy() # Timestamp array in datetime
        self.ts_value           = tt_avg.to_numpy()   # Soil moisture value array

        # Convert timestamp in seconds
        ts_timestamp_ns = self.ts_time_datetime - np.full((len(self.ts_time_datetime),), np.datetime64('1970-01-01T00:00:00Z'))
        ts_timestamp_d = ts_timestamp_ns.astype('timedelta64[D]')

        # Define a vairable
        self.ts_time_d = ts_timestamp_d.astype('int')  # Timestamp array in dates

        # self.timestep = hourly # TODO: make it flexible later
        # TODO: Moving average for 7days or 30days

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

    def movmean(self):
        # detrend the ts_value using moving average using 1-year window
        windowsize = 7  # 1wk just for now
        halfwindow = int(np.rint(windowsize / 2))
        ts_value_movmean = np.convolve(self.ts_value, np.ones(windowsize) / (windowsize), mode='same')

        # repeat the first and last observaions to prevent data loss
        ts_value_movmean[0:halfwindow] = ts_value_movmean[halfwindow]
        ts_value_movmean[len(ts_value_movmean) - halfwindow:-1] = ts_value_movmean[len(ts_value_movmean) - halfwindow]
        ts_value_movmean[-1] = ts_value_movmean[len(ts_value_movmean) - halfwindow]

        self.ts_value = ts_value_movmean
        self.tt = pd.Series(self.ts_value, index=self.ts_time_datetime)

        """
        # plot
        plt.plot(self.ts_value, label='original')
        plt.plot(ts_value_movmean, label='movmean')
        plt.show()
        plt.legend()
        """

    def calc_sinecurve(self):
        # https://scipy-lectures.org/intro/scipy/auto_examples/plot_curve_fit.html
        # Get the sine curve information from insitu data
        pseudo_idx = np.arange(start=0, step=1, stop=len(self.ts_value))
        # pseudo_idx = list(range(0, len(self.ts_value), 1))
        params, params_covariance = optimize.curve_fit(sine_func, xdata=self.ts_time_d, ydata=self.ts_value,
                                                       p0=[1, 0, 0.5])
        phi = params[1]

        # Get the transition valley
        # note that phi sign is opposite from MATLAB code
        sine_n = int(np.round(len(self.ts_time_d)/365))
        sine_start0 = np.floor(self.ts_time_d[0]/365 - phi/2/np.pi)
        sine_start_v = 365/2/np.pi * (2*sine_start0*np.pi + np.pi/2 + phi)
        valley = np.arange(start=sine_start_v, step=365, stop = sine_start_v+ 365* (sine_n+1))
        t_valley = np.full((len(valley),), np.datetime64('1970-01-01T00:00:00Z')) + np.full((len(valley),), np.timedelta64(1,'D')) * valley
        t_valley = pd.Series(t_valley)

        # Plot and confirm
        if self.plot_results:
            plt.figure(figsize=(6, 4))
            plt.scatter(self.ts_time_d, self.ts_value, label='Data')
            plt.plot(self.ts_time_d, sine_func(self.ts_time_d, params[0], params[1], params[2]),
                     label='Fitted function', color='k')
            plt.scatter(valley, sine_func(valley, params[0], params[1], params[2]), color='k')
            plt.legend(loc='best')
            plt.show()

        # return the dates
        return t_valley

    def calc_fcwp(self):
        # TODO: Calculate fc and wp to constrain the SM ts_value better
        None

    def calc_seasontrans(self, t_valley):
        self.t_valley = t_valley
        # initialization
        P0_d2w = [0, 0.0005, 60, 120]
        P0_w2d = [0.5, -0.0005, 120, 80]
        trans_type = ["dry2wet", "wet2dry"]

        seasontrans_date = np.empty((len(self.t_valley), 4))
        seasontrans_date[:] = np.nan

        ## Main execusion
        for trans in range(len(trans_type)):

            # Loop for water years
            for i in range(len(self.t_valley)-1):
                """
                Crop the timeseries
                """

                # Get the base start/end days
                if trans_type[trans] == "dry2wet":
                    trans_start0 = self.t_valley[i]
                    trans_end0 = trans_start0 + datetime.timedelta(days=365/2)
                elif trans_type[trans] == "wet2dry":
                    trans_start0 = self.t_valley[i] + datetime.timedelta(days=365/2)
                    trans_end0 = self.t_valley[i+1]
                # print(trans_start0, trans_end0)

                # Crop the season with 1 month buffer
                mask = (self.tt.index >= trans_start0 - datetime.timedelta(days=30)) & (self.tt.index <= trans_end0 + datetime.timedelta(days=30))
                seasonsm = self.tt.loc[mask]
                seasonsm_value = seasonsm.to_numpy()

                """
                # TODO: Add more cropping treatment ... maybe not needed for this data. seasontrans is quite stable. 
                # If data has too much NaN, or timeseires is empty, do nothing
                if np.count_nonzero(np.isnan(seasonsm_value))/len(seasonsm_value) > 0.3 \
                        or seasonsm_value.size == 0:
                    None
                else:
                    # Try finding the actual wettest & driest point and crop based on it
                    # Get the half length of the timeseries
                    nhalf = int(np.floor(len(seasonsm)/2))
                    if trans_type[trans] == "dry2wet":
                        # the driest point should be happening in the first half of the timeseries
                        t_start = seasonsm[1:nhalf].idxmin()
                        # As the dry period is short, do not use the wettest point to crop the timeseries
                        t_end = trans_end0
                        # Create the mask
                        mask = (self.tt.index >= t_start - datetime.timedelta(days=30)) & (
                                self.tt.index <= t_end + datetime.timedelta(days=30))
                    elif trans_type[trans] == "wet2dry":
                        # The wettest point should be happening in the first half of the timeseries
                        t_start = seasonsm[1:nhalf].idxmax()
                        # The driest point should be happening in the later half of the timeseries
                        t_end = seasonsm[nhalf:-1].idxmax()
                        mask = (self.tt.index >= t_start - datetime.timedelta(days=45)) & (
                                self.tt.index <= t_end + datetime.timedelta(days=15))
                    # re-crop the timeseries with buffer
                    seasonsm = self.tt.loc[mask]
                    seasonsm_value = seasonsm.to_numpy()
                """

                """
                Actual signature calculation
                """

                # If the data has too much NaN, skip the analysis
                if np.count_nonzero(np.isnan(seasonsm_value))/len(seasonsm_value) > 0.3 \
                        or seasonsm_value.size == 0\
                        or len(seasonsm_value) < 90:
                    # print('data is not good')
                    seasontrans_date[i,:] = np.nan
                else:
                # IF the data looks good, exesute the analysis
                    # print('data is good')

                    x = np.arange(start=1, step=1, stop=len(seasonsm_value)+1)
                    y = seasonsm_value
                    if trans_type[trans] == "dry2wet":
                        P0 = P0_d2w
                        Plb = [-5,  0,   0,   1]
                        Pub = [1.5, 0.1, 150, 200]
                    elif trans_type[trans] == "wet2dry":
                        P0 = P0_w2d
                        Plb =  [-1,  -0.1, 0,      1]
                        Pub =  [2.0, 0,    150,    200]

                    popt, pcov = curve_fit(piecewise_linear, x, y, p0=P0)
                    Pfit = popt
                    print(Pfit)
                    # print(res.fun)

                    """
                    # Show response surfaces
                    # sample input range uniformly at 0.1 increments
                    p1axis = np.arange(Plb[1], Pub[1], 0.001)
                    p2axis = np.arange(Plb[3], Pub[3], 1)
                    # create a mesh from the axis
                    p1, p2 = np.meshgrid(p1axis, p2axis)
                    results = np.full(p1.shape, 0)
                    # compute targets
                    for i in range(p1.shape[0]):
                        for j in range(p1.shape[1]):
                            P_01 = [P0[0], p1[i][j],P0[2],  p2[i][j], P0[4], P0[5]]
                            results[i][j] = piecewise_linear_residuals(P_01, x, y)
                    # create a surface plot with the jet color scheme
                    figure = plt.figure()
                    axis = figure.gca(projection='3d')
                    axis.plot_surface(p1, p2, results, cmap='jet')
                    # show the plot
                    plt.show()
                    """

                    # If the wp and fc coincides, or transition is shorter than 7 days, reject it (optimization is likely to have failed)
                    # TODO: modify
                    if Pfit[3] < 7:# abs(Pfit[5]-Pfit[4])<1.0e-03 or Pfit[3] < 7:
                        # if abs(Pfit[5]-Pfit[4])<1.0e-03:
                            # None # print('FC and WP coincides')
                        # elif Pfit[3] < 7:
                            # None # print('Duration too short')
                        Pfit[:] = np.nan
                    else:
                        # Get signatures
                        trans_start_result = seasonsm.axes[0][0] + datetime.timedelta(days=Pfit[2])
                        trans_end_result = seasonsm.axes[0][0] + datetime.timedelta(days=Pfit[2]+Pfit[3])
                        # print(trans_start_result, trans_end_result)

                        # Save in the array
                        seasontrans_date[i,2*trans] = trans_start_result.to_julian_date()
                        seasontrans_date[i,1+2*trans] = trans_end_result.to_julian_date()

                        if self.plot_results:
                            plt.figure(figsize=(6, 4))
                            plt.plot(x, piecewise_linear(x, Pfit[0], Pfit[1], Pfit[2], Pfit[3]))
                            plt.plot(x,y)
                            plt.show()

        # print(seasontrans_date)

        # return signatures
        return seasontrans_date

    def compare_results(self, sim, obs):
        None
        # return errors

