# https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.solve_ivp.html#scipy.integrate.solve_ivp

from scipy.integrate import solve_ivp
import numpy as np
import matplotlib.pyplot as plt

storage_max_m = 0.8
storage_threshold_primary_m = 0.5
wltsmc = 0.3
y = 0.5
coeff_primary = 0.5
coeff_secondary = 0.4
PET = 0.5
infilt = 0.05

# conceptual_reservoir_flux_calc(1,   (0.54463605+ 0.5)/2, storage_threshold_primary_m, storage_max_m, coeff_primary, coeff_secondary, infilt)

# ODE for Zone 1
def conceptual_reservoir_flux_calc(t, S, storage_threshold_primary_m, storage_max_m, coeff_primary, coeff_secondary, infilt):
    storage_above_threshold_m = S - storage_threshold_primary_m
    storage_diff = storage_max_m - storage_threshold_primary_m
    storage_ratio = np.minimum(storage_above_threshold_m / storage_diff, 1)
    dS = infilt -1* (coeff_primary + coeff_secondary) * storage_ratio - PET
    return dS

def event_thresh(t,y,storage_threshold_primary_m, storage_max_m, coeff_primary, coeff_secondary, infilt):
    return y[0] - storage_threshold_primary_m
event_thresh.terminal = True
event_thresh.direction = -1

# ODE for Zone 2
def conceptual_reservoir_flux_calc2(t, S, storage_threshold_primary_m, wltsmc, infilt):
    storage_above_threshold_m_paw = S - wltsmc
    storage_diff_paw = storage_threshold_primary_m - wltsmc
    storage_ratio_paw = np.minimum(storage_above_threshold_m_paw/storage_diff_paw, 1) # Equation 11 (Ogden's document).
    dS = infilt - PET * storage_ratio_paw
    return dS

def event_wltsmc(t,y,storage_threshold_primary_m, wltsmc, infilt):
    return y[0] - wltsmc

def event_thresh2(t,y,storage_threshold_primary_m, wltsmc, infilt):
    return y[0] - storage_threshold_primary_m

event_wltsmc.terminal = True
event_wltsmc.direction = -1
event_thresh2.terminal = True
event_thresh2.direction = +1

# ODE for Zone 3
def conceptual_reservoir_flux_calc3(t, reservoir, wltsmc, infilt):
    dS = infilt + reservoir*0
    return dS

def event_wltsmc2(t, y, wltsmc, infilt):
    return y[0] - wltsmc

event_wltsmc2.terminal = True
event_wltsmc2.direction = +1

# Initialization
ts = []
ys = []
sol_case = []
t0 = 0
y0 = [y]

# Solve ODE
while t0 < 1:
    if y0[0] > storage_threshold_primary_m:
        case = 1
        sol = solve_ivp(conceptual_reservoir_flux_calc,
                        t_span = [t0, 1],
                        y0 = y0,
                        args = (storage_threshold_primary_m,0.8,coeff_primary,coeff_secondary, infilt),
                        events = event_thresh,
                        dense_output=True,
                        method = 'Radau',
                        max_step = 0.5)

    elif (y0[0] <= storage_threshold_primary_m) and (y0[0]  >= wltsmc):
        case = 2
        sol = solve_ivp(conceptual_reservoir_flux_calc2,
                        t_span=[t0, 1],
                        y0=y0,
                        args=(storage_threshold_primary_m, wltsmc, infilt),
                        events = (event_wltsmc, event_thresh2),
                        dense_output=True,
                        method = 'Radau',
                        max_step = 0.5)

    elif y0[0]  < wltsmc:
        case = 3
        sol = solve_ivp(conceptual_reservoir_flux_calc3,
                        t_span=[t0, 1],
                        y0=y0,
                        args=(wltsmc, infilt),
                        events = event_wltsmc2,
                        dense_output=True,
                        method = 'Radau',
                        max_step = 0.5)

    # record results
    ts.append(sol.t[:-1])  # so that it do not overlap with the 2nd solution
    ys.append(sol.y[0][:-1])
    sol_case.append(np.repeat(case, sol.y.size-1))

    # init next loop
    y0 = sol.y[:,-1].copy()
    t0 = sol.t[-1]

# finalize results
ts.append(np.array([1]))
ys.append(np.array([sol.y[0][-1]]))
ts_concat = np.concatenate(ts, axis=0)
ys_concat = np.concatenate(ys, axis=0)
sol_case_concat = np.concatenate(sol_case, axis=0)

# Calculate fluxes
t_proportion = np.diff(ts_concat)
ys_avg =  np.convolve(ys_concat,np.ones(2),'valid')/2
ts_avg = np.convolve(ts_concat,np.ones(2),'valid')/2

lateral_flux = np.zeros(ys_avg.shape)
lateral_flux[sol_case_concat==1] = coeff_secondary * np.minimum((ys_avg[sol_case_concat==1] - storage_threshold_primary_m)/(storage_max_m-storage_threshold_primary_m),1)
# lateral_flux_frac = lateral_flux[:-1] * t_proportion
lateral_flux_frac = lateral_flux* t_proportion

perc_flux = np.zeros(ys_avg.shape)
perc_flux[sol_case_concat==1] = coeff_primary *  np.minimum((ys_avg[sol_case_concat==1] - storage_threshold_primary_m)/(storage_max_m-storage_threshold_primary_m),1)
# perc_flux_frac = perc_flux[:-1]  * t_proportion
perc_flux_frac = perc_flux* t_proportion

et_from_soil = np.repeat(PET, ys_avg.shape)
et_from_soil[sol_case_concat==3] = 0
et_from_soil[sol_case_concat==2] = PET * np.minimum((ys_avg[sol_case_concat==2] - wltsmc)/(storage_threshold_primary_m -wltsmc), 1)
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