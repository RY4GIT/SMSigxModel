import time
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import optimize
import datetime
from scipy.optimize import Bounds, basinhopping, curve_fit


def piecewise_linear(x, P0, P1, P2, P3):
    y0 = np.where(x-P2>0, P0+P1*x, P0+P1*P2)
    return np.where(x-(P2+P3)>0, P0+P1*(P2+P3), y0)


# goal

# condition
"""
P_d2w = [0, 0.0003, 80, 100]
P0_d2w = [0, 0.0005, 60, 120]
Plb = [-5, 0, 0, 1]
Pub = [1.5, 0.1, 150, 200]
"""

P_w2d = [0.1, -0.0003, 60, 100]
P0_w2d = [0.5, -0.0005, 120, 80]
Plb = [-5, -0.1, 0, 1]
Pub = (2, 0, 150, 200)
Plb = (-5, -0.1, 0, 1)

x_obs = np.arange(1,300,1)
y_obs = piecewise_linear(x_obs, P_w2d[0], P_w2d[1], P_w2d[2], P_w2d[3]) + np.random.normal(size=299) * 0.0002

bounds = Bounds(lb=Plb, ub=Pub)

popt, pcov = curve_fit(piecewise_linear, x_obs, y_obs, p0=P0_w2d, bounds=(Plb, Pub))

print("%.2f, %.5f, %.2f, %.2f" % (popt[0], popt[1], popt[2], popt[3]))

plt.plot(x_obs,y_obs, label='original')
plt.plot(x_obs, piecewise_linear(x_obs, popt[0], popt[1], popt[2], popt[3]), label='calculated')
# plt.plot(x_obs, y_obs2, label='goal')
plt.legend()