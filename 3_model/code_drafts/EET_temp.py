import numpy as np
from numpy import matlib as mb
import matplotlib.pyplot as plt
import os

"""
                             mu   mu_star     sigma  mu_star_conf
bb                    -0.053354  0.063245  0.118680      0.010499
satdk                  0.702110  0.702110  1.033641      0.089968
satpsi                -0.059611  0.063160  0.121317      0.011089
slop                   0.000000  0.000000  0.000000      0.000000
smcmax                 0.122701  0.126175  0.229837      0.020349
wltsmc                -0.102437  0.106028  0.206397      0.017851
D                      0.736485  0.736744  0.879889      0.080392
coeff_secondary        0.000000  0.000000  0.000000      0.000000
exponent_secondary     0.000000  0.000000  0.000000      0.000000
max_gw_storage         0.000000  0.000000  0.000000      0.000000
Cgw                    0.000000  0.000000  0.000000      0.000000
expon                  0.000000  0.000000  0.000000      0.000000
K_nash                 0.000000  0.000000  0.000000      0.000000
refkdt                 0.000000  0.000000  0.000000      0.000000
trigger_z_m           -0.172239  0.173228  0.555678      0.048789
fc_atm_press_fraction  0.571610  0.572939  0.895315      0.080978
"""
Si = {
    'names': ['bb',
              'satdk',
              'satpsi',
              'slop',
              'smcmax',
              'wltsmc',
              'D',
              'coeff_secondary',
              'exponent_secondary',
              'max_gw_storage',
              'Cgw',
              'expon',
              'K_nash',
              'refkdt',
              'trigger_z_m',
              'fc_atm_press_fraction'
              ],
    'mu_star': [0.063245,
                0.702110,
                0.063160,
                0.000000,
                0.126175,
                0.106028,
                0.736744,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.000000,
                0.173228,
                0.572939],
    'sigma': [0.118680,
              1.033641,
              0.121317,
              0.000000,
              0.229837,
              0.206397,
              0.879889,
              0.000000,
              0.000000,
              0.000000,
              0.000000,
              0.000000,
              0.000000,
              0.000000,
              0.555678,
              0.895315]
}

problem = {
    'num_vars': 16,
    'names': ['bb',
              'satdk',
              'satpsi',
              'slop',
              'smcmax',
              'wltsmc',
              'D',
              'coeff_secondary',
              'exponent_secondary',
              'max_gw_storage',
              'Cgw',
              'expon',
              'K_nash',
              'refkdt',
              'trigger_z_m',
              'fc_atm_press_fraction'
              ],
    'bounds': [[2, 15],  # bb
               [0, 1],  # satdk
               [0.02, 0.78],  # satpsi
               [0, 1],  # slop
               [0.33, 0.7],  # smcmax
               [0, 0.5],  # wltsmc
               [0.01, 2],  # D
               [0.01, 3],  # coeff_secondary
               [1, 8],  # exponent_secondary
               [10, 250],  # max_gw_storage
               [0.01, 3],  # Cgw
               [1, 8],  # expon
               [0, 1],  # K_nash
               [0.1, 4],  # refkdt
               [0.01, 0.87],  # trigger_z_m
               [0.01, 0.33]  # fc_atm_press_fraction
               ]
}



# Options for the graphic
pltfont = {'fontname': 'DejaVu Sans', 'fontsize': 15}  # font for axes
pltfont_leg = {'family': 'DejaVu Sans', 'size': 15}  # font for legend

ms = 10  # Marker size
col = np.array([[228, 26, 28], [55, 126, 184], [77, 175, 74],
                [152, 78, 163], [255, 127, 0]]) / 256
A = len(col)
L = int(np.ceil(problem['num_vars'] / A))
clrs = mb.repmat(col, L, 1)

# Plot elementary effects and std's for Morris analysis
fig = plt.figure()

# First plot EEs mean & std as circles:
for i in range(len(Si['mu_star'])):
    plt.plot(
        Si['mu_star'][i], Si['sigma'][i], 'ok', markerfacecolor=clrs[i],
        markersize=ms, markeredgecolor='k'
    )
    plt.text(Si['mu_star'][i] + 0.02, Si['sigma'][i] + 0.02, Si['names'][i], fontsize=9)

# Create legend:
plt.legend(Si['names'], bbox_to_anchor=(1.01, 1), loc='upper left', prop=pltfont_leg)

plt.xlabel('Mean of EEs', **pltfont)
plt.ylabel('Standard deviation of EEs', **pltfont)
plt.grid(linestyle='--')
plt.xticks(**pltfont)
plt.yticks(**pltfont)
fig.set_size_inches(7, 7)

out_fn = 'test_EET.png'
out_path = '..\\4_out\\Mahurangi'
out_path_plot = os.path.join(out_path, 'plot_SALib')
if not os.path.exists(out_path_plot):
    # Create a new directory because it does not exist
    os.makedirs(out_path_plot)
plt.savefig(os.path.join(out_path_plot, out_fn), dpi=600, format='png')
plt.show()
