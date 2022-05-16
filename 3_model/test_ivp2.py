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
infilt = 0.5

def conceptual_reservoir_flux_calc2(t, reservoir, storage_threshold_primary_m, wltsmc):
    storage_above_threshold_m_paw = reservoir - wltsmc
    storage_diff_paw = storage_threshold_primary_m - wltsmc
    storage_ratio_paw = storage_above_threshold_m_paw / storage_diff_paw  # Equation 11 (Ogden's document).
    dS = - PET * storage_ratio_paw
    return dS

def event_wltsmc(t,y,storage_threshold_primary_m, wltsmc):
    return y[0] - wltsmc

sol = solve_ivp(conceptual_reservoir_flux_calc2,
                t_span=[0, 1], y0=y0,
                args=(storage_threshold_primary_m, wltsmc),
                events = event_wltsmc,
                dense_output=True,
                method = 'Radau')

t_proportion = np.diff(sol.t)

et_from_soil = PET * (sol.y[0] - wltsmc)/(storage_threshold_primary_m -wltsmc) * t_proportion
et_from_soil_frac = et_from_soil[:-1] * t_proportion

plt.plot(sol.t[:-1], et_from_soil_frac,label =  'ET from soil', color =  'mediumseagreen')
plt.plot(sol.t[:-1], np.diff(-sol.y[0],axis=0),label =  'dStorage',  ls= '--',color='orangered')
# plt.plot(ts_concat,ys_concat,label =  'Storage', color='slateblue')
plt.xlabel('timestep')
plt.ylabel('Storage')
plt.title('Soil coceptual reservoir')
plt.legend()
plt.show()

flux_scale = np.diff(-sol.y[0],axis=0)/et_from_soil_frac
scaled_et_flux = et_from_soil_frac* flux_scale

plt.plot(ts_concat[:-1], scaled_et_flux,label =  'ET from soil', color =  'mediumseagreen')
plt.plot(ts_concat[:-1], np.diff(-sol.y[0],axis=0),label =  'dStorage',  ls= '--',color='orangered')
# plt.plot(ts_concat,ys_concat,label =  'Storage', color='slateblue')
plt.xlabel('timestep')
plt.ylabel('Storage')
plt.title('Soil coceptual reservoir')
plt.legend()
plt.show()

"""
sol = solve_ivp(conceptual_reservoir_flux_calc,
                t_span = [0, 1],
                y0 = y0,
                args = (storage_threshold_primary_m,
                        storage_max_m,
                        coeff_primary,
                        coeff_secondary),
                dense_output=True)


sol = solve_ivp(conceptual_reservoir_flux_calc,
                t_span = [0, 1], y0 = y0,
                args = (storage_threshold_primary_m,0.8,coeff_primary,coeff_secondary))

t = sol.t
y = sol.y[0]
plt.plot(t, y)
plt.xlabel('timestep')
plt.ylabel('Storage')
plt.title('Soil coceptual reservoir')
plt.show()

# add Evapo in the dS, if not a same form as Perc and Lat, use breakpoint to calculate three fluxes

dS = np.diff(y)
p_flux = -1 * dS * coeff_primary/(coeff_primary+coeff_secondary)
s_flux = -1 * dS * coeff_secondary/(coeff_primary+coeff_secondary)
ds_cum = y[-1] - y[0]
p_flux = -1 * ds_cum * coeff_primary/(coeff_primary+coeff_secondary)
s_flux = -1 * ds_cum * coeff_secondary/(coeff_primary+coeff_secondary)

# plt.plot(t[1:],p_flux)
# plt.plot(t[1:],s_flux)



storage_threshold_primary_m = 0.5
y0 = [0.4]

def conceptual_reservoir_flux_calc(t, reservoir, storage_threshold_primary_m, storage_max_m, coeff_primary, coeff_secondary):
    storage_above_threshold_m = reservoir - storage_threshold_primary_m
    storage_diff = storage_max_m - storage_threshold_primary_m
    storage_ratio = storage_above_threshold_m / storage_diff  # Equation 11 (Ogden's document).
    dS = (coeff_primary + coeff_secondary) * storage_ratio
    return reservoir - dS

def event(t, reservoir, storage_threshold_primary_m=storage_threshold_primary_m):
    return reservoir[0] - storage_threshold_primary_m

storage_event = lambda t,y:event(t,y)
storage_event.terminal = True
storage_event.direction = -1
sol = solve_ivp(lambda t, y:conceptual_reservoir_flux_calc(t, y, storage_threshold_primary_m=storage_threshold_primary_m,
                                                           storage_max_m=0.8,
                                                           coeff_primary=1,
                                                           coeff_secondary=1),
                t_span = [0, 1], y0=y0,
                dense_output=False, method='Radau', events=storage_event)

# https://stackoverflow.com/questions/50701425/how-to-use-if-statement-in-a-differential-equation-scipy
# https://stackoverflow.com/questions/55614547/how-to-pass-parameters-to-event-functions-inscipy-integrate-solve-ivp
plt.plot(sol.t, sol.y[0])
plt.xlabel('t')
plt.legend('reservoir', shadow=True)
plt.title('Soil coceptual reservoir')
plt.show()



"""




