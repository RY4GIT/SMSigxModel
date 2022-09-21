# https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.solve_ivp.html#scipy.integrate.solve_ivp

from scipy.integrate import solve_ivp
import numpy as np
import matplotlib.pyplot as plt
import timeit
import nbkode

start = timeit.timeit()
storage_max_m = 0.8
storage_threshold_primary_m = 0.5
wltsmc = 0.3
y = 0.8
coeff_primary = 0.5
coeff_secondary = 0.4
PET = 0.5
infilt = 0.

# conceptual_reservoir_flux_calc(1,   (0.54463605+ 0.5)/2, storage_threshold_primary_m, storage_max_m, coeff_primary, coeff_secondary, infilt)

# ODE for Zone 1
def conceptual_reservoir_flux_calc(t, S, dS, p):
    """
    storage_threshold_primary_m = p[0]
    storage_max_m = p[1]
    coeff_primary = p[2]
    coeff_secondary = p[3]
    PET = p[4]
    infilt = p[5]
    wltsmc = p[6]
    """
    S_ = nb.carray(S, (1,))
    p_ = nb.carray(p, (1,))

    storage_above_threshold_m = S_ - p_[0]
    storage_diff = p_[1] - p_[0]
    storage_ratio = np.minimum(storage_above_threshold_m / storage_diff, 1)

    perc_lat_switch = np.multiply(S_ - p_[0] > 0, 1)
    ET_switch = np.multiply(S_ - p_[6] > 0, 1)
    # print(perc_lat_switch, ET_switch)

    storage_above_threshold_m_paw = S_ -  p_[6]
    storage_diff_paw = p_[0] -  p_[6]
    storage_ratio_paw = np.minimum(storage_above_threshold_m_paw/storage_diff_paw, 1) # Equation 11 (Ogden's document).

    dS = p_[5] -1 * perc_lat_switch * (p_[2] + p_[3]) * storage_ratio - ET_switch * p_[4] * storage_ratio_paw
    return dS

def jac(t, S, storage_threshold_primary_m, storage_max_m, coeff_primary, coeff_secondary, PET, infilt, wltsmc):
    storage_diff = storage_max_m - storage_threshold_primary_m

    perc_lat_switch = np.multiply(S - storage_threshold_primary_m > 0, 1)
    ET_switch = np.multiply((S - wltsmc > 0) and (S - storage_threshold_primary_m < 0), 1)

    storage_diff_paw = storage_threshold_primary_m - wltsmc

    dfdS = -1 * perc_lat_switch * (coeff_primary + coeff_secondary) * 1/storage_diff - ET_switch * PET * 1/storage_diff_paw
    return [dfdS]


# Initialization
ts = []
ys = []
sol_case = []
t0 = 0
y0 = [y]
t = np.linspace(0, 1, 11)
p = np.array([storage_threshold_primary_m, storage_max_m, coeff_primary, coeff_secondary, PET, infilt, wltsmc])
solver = nbkode.Heun3(conceptual_reservoir_flux_calc, t0, y0, params=p)
t, y = solver.run(t)
