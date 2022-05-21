# https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.solve_ivp.html#scipy.integrate.solve_ivp

from scipy.integrate import solve_ivp, odeint
import numpy as np
import matplotlib.pyplot as plt
import time


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
def conceptual_reservoir_flux_calc(t, S, storage_threshold_primary_m, storage_max_m, coeff_primary, coeff_secondary, PET, infilt, wltsmc):
    storage_above_threshold_m = S - storage_threshold_primary_m
    storage_diff = storage_max_m - storage_threshold_primary_m
    storage_ratio = np.minimum(storage_above_threshold_m / storage_diff, 1)

    perc_lat_switch = np.multiply(S - storage_threshold_primary_m > 0, 1)
    ET_switch = np.multiply(S - wltsmc > 0, 1)
    # print(perc_lat_switch, ET_switch)

    storage_above_threshold_m_paw = S - wltsmc
    storage_diff_paw = storage_threshold_primary_m - wltsmc
    storage_ratio_paw = np.minimum(storage_above_threshold_m_paw/storage_diff_paw, 1) # Equation 11 (Ogden's document).

    dS = infilt -1 * perc_lat_switch * (coeff_primary + coeff_secondary) * storage_ratio - ET_switch * PET * storage_ratio_paw
    return dS

def conceptual_reservoir_flux_calc_odeint(S, t, storage_threshold_primary_m, storage_max_m, coeff_primary, coeff_secondary, PET, infilt, wltsmc):
    storage_above_threshold_m = S - storage_threshold_primary_m
    storage_diff = storage_max_m - storage_threshold_primary_m
    storage_ratio = np.minimum(storage_above_threshold_m / storage_diff, 1)

    perc_lat_switch = np.multiply(S - storage_threshold_primary_m > 0, 1)
    ET_switch = np.multiply(S - wltsmc > 0, 1)
    # print(perc_lat_switch, ET_switch)

    storage_above_threshold_m_paw = S - wltsmc
    storage_diff_paw = storage_threshold_primary_m - wltsmc
    storage_ratio_paw = np.minimum(storage_above_threshold_m_paw/storage_diff_paw, 1) # Equation 11 (Ogden's document).

    dS = infilt -1 * perc_lat_switch * (coeff_primary + coeff_secondary) * storage_ratio - ET_switch * PET * storage_ratio_paw
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

start =time.perf_counter()
sol = solve_ivp(conceptual_reservoir_flux_calc,
                t_span=[t0, 1],
                y0=y0,
                args=(storage_threshold_primary_m, storage_max_m, coeff_primary, coeff_secondary, PET, infilt, wltsmc),
                dense_output=True,
                method='Radau')
end =time.perf_counter()


start2 = time.perf_counter()
t = np.linspace(0, 1, 6)
sol = odeint(
    conceptual_reservoir_flux_calc,
             y0,
             t,
             args=(storage_threshold_primary_m, storage_max_m, coeff_primary, coeff_secondary, PET, infilt, wltsmc),
             tfirst=True
             )
end2 = time.perf_counter()

start3 =time.perf_counter()
sol = solve_ivp(conceptual_reservoir_flux_calc,
                t_span=[t0, 1],
                y0=y0,
                args=(storage_threshold_primary_m, storage_max_m, coeff_primary, coeff_secondary, PET, infilt, wltsmc),
                dense_output=True,
                method='Radau',
                jac=jac)
end3 =time.perf_counter()


start4 = time.perf_counter()
t = np.linspace(0, 1, 6)
sol = odeint(
    conceptual_reservoir_flux_calc,
             y0,
             t,
             args=(storage_threshold_primary_m, storage_max_m, coeff_primary, coeff_secondary, PET, infilt, wltsmc),
             tfirst=True,
            Dfun=jac
             )
end4 = time.perf_counter()

print(f"solve_ivp took {end-start} sec")
print(f"solve_ivp with Jac took {end3-start3} sec")
print(f"odeint took {end2-start2}sec")
print(f"odeint with Jac took {end4-start4}sec")
