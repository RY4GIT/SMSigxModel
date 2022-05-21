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

"""
sol = solve_ivp(conceptual_reservoir_flux_calc,
                t_span=[t0, 1],
                y0=y0,
                args=(storage_threshold_primary_m, storage_max_m, coeff_primary, coeff_secondary, PET, infilt, wltsmc),
                dense_output=True,
                method='LSODA',
                max_step=0.5,
                jac=jac,
                lband=0,
                uband=0)

# finalize results
ts_concat = sol.t
ys_concat = sol.y[0]

# Calculate fluxes
t_proportion = np.diff(ts_concat)
ys_avg =  np.convolve(ys_concat,np.ones(2),'valid')/2
ts_avg = np.convolve(ts_concat,np.ones(2),'valid')/2

lateral_flux = np.zeros(ys_avg.shape)
perc_lat_switch = ys_avg - storage_threshold_primary_m > 0
lateral_flux[perc_lat_switch] = coeff_secondary * np.minimum((ys_avg[perc_lat_switch] - storage_threshold_primary_m)/(storage_max_m-storage_threshold_primary_m),1)
lateral_flux_frac = lateral_flux* t_proportion

perc_flux = np.zeros(ys_avg.shape)
perc_flux[perc_lat_switch] = coeff_primary *  np.minimum((ys_avg[perc_lat_switch] - storage_threshold_primary_m)/(storage_max_m-storage_threshold_primary_m),1)
perc_flux_frac = perc_flux* t_proportion

et_from_soil = np.zeros(ys_avg.shape)
ET_switch = ys_avg - wltsmc > 0
# et_from_soil[sol_case_concat] = 0
et_from_soil[ET_switch] = PET * np.minimum((ys_avg[ET_switch] - wltsmc)/(storage_threshold_primary_m -wltsmc), 1)
# et_from_soil_frac = et_from_soil[:-1]  * t_proportion
et_from_soil_frac = et_from_soil * t_proportion

infilt_to_soil = np.repeat(infilt, ys_avg.shape)
# infilt_to_soil_frac = infilt_to_soil[:-1] * t_proportion
infilt_to_soil_frac = infilt_to_soil* t_proportion

# Plot fluxes before scaling
plt.subplot(2,1,1)
plt.plot(ts_avg, -lateral_flux_frac/t_proportion, label = 'Lateral flux', color = 'darkturquoise')
plt.plot(ts_avg, -perc_flux_frac/t_proportion, label = 'Percolation flux', color = 'royalblue')
plt.plot(ts_avg, -et_from_soil_frac/t_proportion,label =  'ET from soil', color =  'mediumseagreen')
plt.plot(ts_avg, +infilt_to_soil_frac/t_proportion, label =  'Infilt to soil', color = 'blue')
plt.plot(ts_avg, np.diff(ys_concat)/t_proportion,label =  'dStorage', color='orange', ls= '-.')
plt.plot(ts_avg, (-lateral_flux_frac-perc_flux_frac-et_from_soil_frac+infilt_to_soil_frac)/t_proportion, label = 'Flux sum', ls= '--', color =  'slateblue')
plt.xlabel('timestep [h]')
plt.ylabel('Flux [m/h]')
plt.title('Soil coceptual reservoir')
plt.legend()
plt.subplot(2,1,2)
plt.plot(ts_concat, ys_concat,label =  'Storage', color='orange', ls= '-')
plt.xlabel('timestep')
plt.ylabel('Storage')
plt.legend()
plt.show()

# Scale fluxes
sum_outflux = lateral_flux_frac+perc_flux_frac+et_from_soil_frac
if sum_outflux.any() == 0:
    flux_scale = 0
else:
    flux_scale = np.zeros(infilt_to_soil_frac.shape)
    flux_scale[sum_outflux!=0] = (np.diff(-ys_concat,axis=0)[sum_outflux!=0]+infilt_to_soil_frac[sum_outflux!=0])/sum_outflux[sum_outflux!=0]
    flux_scale[sum_outflux == 0] = 0
scaled_lateral_flux = lateral_flux_frac * flux_scale
scaled_perc_flux = perc_flux_frac* flux_scale
scaled_et_flux = et_from_soil_frac* flux_scale

# Plot scaled fluxes
plt.subplot(2,1,1)
plt.plot(ts_avg, -scaled_lateral_flux/t_proportion, label = 'Lateral flux', color = 'darkturquoise')
plt.plot(ts_avg, -scaled_perc_flux/t_proportion, label = 'Percolation flux', color = 'royalblue')
plt.plot(ts_avg, -scaled_et_flux/t_proportion,label =  'ET from soil', color =  'mediumseagreen')
plt.plot(ts_avg, +infilt_to_soil_frac/t_proportion, label =  'Infilt to soil', color = 'blue')
plt.plot(ts_avg, np.diff(ys_concat,axis=0)/t_proportion,label =  'dStorage', color='orange', ls= '-.')
plt.plot(ts_avg, (-scaled_lateral_flux-scaled_perc_flux-scaled_et_flux+infilt_to_soil_frac)/t_proportion, label = 'Flux sum', ls= '--', color =  'slateblue')
# plt.plot(ts_concat,ys_concat,label =  'Storage', color='slateblue')
plt.xlabel('timestep [h]')
plt.ylabel('Flux [m/h]')
plt.title('Soil coceptual reservoir')
plt.legend()
plt.subplot(2,1,2)
plt.plot(ts_concat, ys_concat,label =  'Storage', color='orange', ls= '-')
plt.xlabel('timestep')
plt.ylabel('Storage')
plt.legend()
plt.show()

print(f"Lateral flux: {sum(scaled_lateral_flux)}")
print(f"Percolation flux: {sum(scaled_perc_flux)}")
print(f"ET from soil: {sum(scaled_et_flux)}")
print(f"Infiltration to soil: {sum(infilt_to_soil_frac)}")

end = timeit.timeit()
print(f"took {end-start}sec")

"""