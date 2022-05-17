# https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.solve_ivp.html#scipy.integrate.solve_ivp

from scipy.integrate import solve_ivp
import numpy as np
import matplotlib.pyplot as plt

storage_threshold_primary_m = 0.5
storage_max_m = 0.5
wltsmc = 0.45
y0 = [0.6]
coeff_primary = 0.5
coeff_secondary = 0.4
PET = 0.5
infilt = 0.3

def conceptual_reservoir_flux_calc2(t, reservoir, storage_threshold_primary_m, wltsmc, infilt):
    storage_above_threshold_m_paw = reservoir - wltsmc
    storage_diff_paw = storage_threshold_primary_m - wltsmc
    storage_ratio_paw = storage_above_threshold_m_paw / storage_diff_paw  # Equation 11 (Ogden's document).
    dS = infilt - PET * storage_ratio_paw
    return dS

def event_wltsmc(t,y,storage_threshold_primary_m, wltsmc, infilt):
    return y[0] - wltsmc

sol = solve_ivp(conceptual_reservoir_flux_calc2,
                t_span=[0, 1], y0=y0,
                args=(storage_threshold_primary_m, wltsmc, infilt),
                events = event_wltsmc,
                dense_output=True,
                method = 'Radau')

t_proportion = np.diff(sol.t)

et_from_soil = PET * (sol.y[0][:-1] - wltsmc)/(storage_threshold_primary_m -wltsmc) * t_proportion
et_from_soil_frac = et_from_soil * t_proportion

infilt_to_soil = np.repeat(infilt,et_from_soil_frac.shape)
infilt_to_soil_frac = infilt_to_soil* t_proportion

plt.plot(sol.t[:-1], et_from_soil_frac,label =  'ET from soil', color =  'mediumseagreen')
plt.plot(sol.t[:-1], np.diff(-sol.y[0],axis=0),label =  'dStorage',  ls= '--',color='orangered')
plt.plot(sol.t[:-1], -infilt_to_soil_frac,label =  'infilt')

plt.xlabel('timestep')
plt.ylabel('Storage')
plt.title('Soil coceptual reservoir')
plt.legend()
plt.show()

flux_scale = (-np.diff(sol.y[0],axis=0)+infilt_to_soil_frac)/et_from_soil_frac
scaled_et_flux = et_from_soil_frac* flux_scale # this is more than evaporative

plt.plot(sol.t[:-1], scaled_et_flux,label =  'ET from soil', color =  'mediumseagreen')
plt.plot(sol.t[:-1], np.diff(-sol.y[0],axis=0),label =  'dStorage',  ls= '--',color='orangered')
plt.plot(sol.t[:-1], infilt_to_soil_frac,label =  'infilt')
# plt.plot(ts_concat,ys_concat,label =  'Storage', color='slateblue')
plt.xlabel('timestep')
plt.ylabel('Storage')
plt.title('Soil coceptual reservoir')
plt.legend()
plt.show()


