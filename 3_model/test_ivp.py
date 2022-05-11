# https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.solve_ivp.html#scipy.integrate.solve_ivp

from scipy.integrate import solve_ivp
import numpy as np
import matplotlib.pyplot as plt

"""
def lotkavolterra(t, z, a, b, c, d):
    x, y = z
    return [a*x - b*x*y, -c*y + d*x*y]

sol = solve_ivp(lotkavolterra, [0, 15], [10, 5], args=(1.5, 1, 3, 1),
                dense_output=True)
t = np.linspace(0, 15, 300)
z = sol.sol(t)
plt.plot(t, z.T)
plt.xlabel('t')
plt.legend(['x', 'y'], shadow=True)
plt.title('Lotka-Volterra System')
plt.show()
"""

storage_threshold_primary_m = 0.7
storage_max_m = 0.8
y0 = [0.6]
coeff_primary = 0.5
coeff_secondary = 0.5

def conceptual_reservoir_flux_calc(t, reservoir, storage_threshold_primary_m, storage_max_m, coeff_primary, coeff_secondary):
    storage_above_threshold_m = reservoir - storage_threshold_primary_m
    storage_diff = storage_max_m - storage_threshold_primary_m
    storage_ratio = storage_above_threshold_m / storage_diff  # Equation 11 (Ogden's document).
    dS = -1* (coeff_primary + coeff_secondary) * storage_ratio
    return dS

sol = solve_ivp(lambda t, y:conceptual_reservoir_flux_calc(t, y, storage_threshold_primary_m, storage_max_m, coeff_primary, coeff_secondary),
                t_span = [0, 1], y0 = y0,
                dense_output=False)

sol = solve_ivp(conceptual_reservoir_flux_calc,
                t_span = [0, 1], y0 = y0,
                args = (storage_threshold_primary_m,0.8,coeff_primary,coeff_secondary))
t = sol.t
y = sol.y[0]
dS = np.diff(y)
p_flux = -1 * dS * coeff_primary/(coeff_primary+coeff_secondary)
s_flux = -1 * dS * coeff_secondary/(coeff_primary+coeff_secondary)
ds_cum = y[-1] - y[0]
p_flux = -1 * ds_cum * coeff_primary/(coeff_primary+coeff_secondary)
s_flux = -1 * ds_cum * coeff_secondary/(coeff_primary+coeff_secondary)
plt.plot(t, y)
# plt.plot(t[1:],p_flux)
# plt.plot(t[1:],s_flux)
plt.xlabel('timestep')
plt.ylabel('Storage')
plt.title('Soil coceptual reservoir')
plt.show()

"""
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




