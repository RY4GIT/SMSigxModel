import os
import spotpy
import pandas as pd
import numpy as np
from numpy import matlib as mb
import scipy
import matplotlib.pyplot as plt
import seaborn as sns
import json
import warnings

from SALib.sample import saltelli
from SALib.analyze import sobol
from SALib.sample import morris as morris_s
from SALib.analyze import morris as morris_a
from SALib.sample import fast_sampler
from SALib.analyze import fast

# %matplotlib inline
import itertools
from math import pi
from matplotlib.legend_handler import HandlerPatch

class SALib_CFE():

    def __init__(self, cfe_instance=None, config_path='', method_SALib=None, like_measure='NashSutcliffe', var_measure='Flow', out_path=None):
        self.cfe_instance = cfe_instance
        self.method_SALib = method_SALib
        self.out_path = out_path
        self.like_measure = like_measure
        self.var_measure = var_measure

        with open(config_path, 'r') as outfile:
            problem = json.load(outfile)
        self.problem = problem

    def run(self):
        if self.method_SALib['method'] == "stability_test":
            # sample

            # Array with n's to use
            nsamples = np.arange(10, 300, 10)

            # run a model
            S1_estimates, ST_estimates, nsample_record = run_cfes_for_stability_test(
                problem = self.problem,
                cfe_instance = self.cfe_instance,
                nrun = nsamples,
                like_measure=self.like_measure,
                var_measure=self.var_measure
            )

            out_path = 'G:/Shared drives/Ryoko and Hilary/SMSigxModel/analysis/4_out/sensitivity_analysis/Mahurangi/sensitivity_stability'
            S1_estimates.to_csv(os.path.join(out_path, 'S1.csv'))
            ST_estimates.to_csv(os.path.join(out_path, 'ST.csv'))
            nsample_record.to_csv(os.path.join(out_path, 'nsample.csv'))

        if self.method_SALib['method'] == "Sobol":
            # sample
            n = self.method_SALib['n']
            self.param_values = saltelli.sample(self.problem, n, calc_second_order=True)

            # run a model
            self.Y = run_cfes(
                problem = self.problem,
                cfe_instance = self.cfe_instance,
                param_values= self.param_values,
                nrun = n*(2*self.problem['num_vars']+2),
                like_measure=self.like_measure,
                var_measure=self.var_measure
            )

            # evaluation
            self.Si = sobol.analyze(self.problem, self.Y, calc_second_order=True, print_to_console=False)
            print(self.Si)

        if self.method_SALib['method'] == "FAST":
            # sample
            n = 8
            m = 1
            self.param_values = fast_sampler.sample(self.problem, N=n, M=m)

            # run a model
            self.Y = run_cfes(
                problem=self.problem,
                cfe_instance=self.cfe_instance,
                param_values=self.param_values,
                nrun=n * (2 * self.problem['num_vars'] + 2),
                like_measure=self.like_measure
            )

            # evaluation
            self.Si = fast.analyze(self.problem, self.Y, M=4, num_resamples=10, conf_level=0.95, print_to_console=False, seed=None)
            print(self.Si)

        if self.method_SALib['method'] == "Morris":
            N = self.method_SALib['N']
            n_levels = self.method_SALib['n_levels']
            self.param_values = morris_s.sample(self.problem, N=N, num_levels=n_levels)

            # run a model
            self.Y = run_cfes(
                problem = self.problem,
                cfe_instance = self.cfe_instance,
                param_values=self.param_values,
                nrun = N * (self.problem['num_vars']+1),
                like_measure = self.like_measure,
                var_measure = self.var_measure
            )

            # evaluation
            self.Si = morris_a.analyze(self.problem, self.param_values, self.Y, print_to_console=True)
            print(self.Si['mu'])

        file = open(os.path.join(self.out_path, "param_bounds.txt"), "w")
        file.write("%s" % str(self.problem))
        file.close

    def plot(self, plot_type=None):
        # Add dotty plot module. Either here or in the above method
        # https://pynetlogo.readthedocs.io/en/latest/_docs/SALib_ipyparallel.html

        if self.method_SALib['plot'] == "EET":
            # Options for the graphic
            pltfont = {'fontname': 'DejaVu Sans', 'fontsize': 15}  # font for axes
            pltfont_leg = {'family': 'DejaVu Sans', 'size': 15}  # font for legend
            ms = 10  # Marker size
            col = np.array([[228, 26, 28], [55, 126, 184], [77, 175, 74],
                            [152, 78, 163], [255, 127, 0]]) / 256
            A = len(col)
            L = int(np.ceil(self.problem['num_vars'] / A))
            clrs = mb.repmat(col, L, 1)

            # Plot elementary effects and std's for Morris analysis
            fig = plt.figure()

            # First plot EEs mean & std as circles:
            for i in range(len(self.Si['mu_star'])):
                plt.plot(
                    self.Si['mu_star'][i], self.Si['sigma'][i], 'ok', markerfacecolor=clrs[i],
                    markersize=ms, markeredgecolor='k'
                )

            # Create legend:
            plt.legend(self.Si['names'], loc='best', prop=pltfont_leg)

            plt.xlabel('Mean of EEs', **pltfont)
            plt.ylabel('Standard deviation of EEs', **pltfont)
            plt.grid(linestyle='--')
            plt.xticks(**pltfont)
            plt.yticks(**pltfont)
            fig.set_size_inches(7, 7)

            out_fn = 'test_EET.png'

        if self.method_SALib['plot'] == "dotty":
            # Dotty plots for any types of sampled parameters
            nrow = 1
            ncol = 3
            fig, ax = plt.subplots(nrow, ncol, sharey=True)
            y = self.Y
            for i, a in enumerate(ax.flatten()):
                x = self.param_values[:, i]
                sns.regplot(
                    x, y, ax=a, ci=None, color='k', scatter_kws={'alpha': 0.2, 's': 4, 'color': 'gray'}
                )
                pearson = scipy.stats.pearsonr(x, y)
                a.annotate("r: {:6.3f}".format(pearson[0]), xy=(0.15, 0.85), xycoords='axes fraction', fontsize=12)
                a.set_ylabel('NSE value')
                if divmod(i, ncol)[1] > 0:
                    a.get_yaxis().set_visible(False)
                a.set_xlabel(self.problem['names'][i])
                a.set_ylim([0, 1.1 * np.max(y)])

            fig.set_size_inches(3*ncol, 4*nrow, forward=True)
            fig.subplots_adjust(wspace=0.2, hspace=0.3)

            out_fn = 'test_dotty.png'

        if self.method_SALib['plot'] == "STS1":
            # Bar plots for Sobol analysis (total order & 1st order indices)
            if self.method_SALib['method'] == 'Sobol':
                Si_filter = {k: self.Si[k] for k in ['ST', 'ST_conf', 'S1', 'S1_conf']}
            elif self.method_SALib['method'] == 'FAST':
                Si_filter = {k: self.Si[k] for k in ['ST', 'S1']}
            Si_df = pd.DataFrame(Si_filter, index=self.problem['names'])

            fig, ax = plt.subplots(1)
            indices = Si_df[['S1', 'ST']]
            if self.method_SALib['method'] == 'Sobol':
                err = Si_df[['S1_conf', 'ST_conf']]
                indices.plot.bar(yerr=err.values.T, ax=ax)
            elif self.method_SALib['method'] == 'FAST':
                indices.plot.bar(ax=ax)
            fig.set_size_inches(8, 4)

            out_fn = 'test_STS1.png'

        if self.method_SALib['plot'] == "radial":
            # TODO: need a debug
            # Radial plot for Sobol analysis
            # IndexError: index 5 is out of bounds for axis 0 with size 5
            # https://pynetlogo.readthedocs.io/en/latest/_docs/SALib_ipyparallel.html

            criterion = 'ST'
            threshold = 0.01
            max_linewidth_s2 = 15  # 25*1.8
            max_s_radius = 0.3

            # prepare data
            # use the absolute values of all the indices
            # sobol_indices = {key:np.abs(stats) for key, stats in sobol_indices.items()}

            # dataframe with ST and S1
            sobol_stats = {key: self.Si[key] for key in ['ST', 'S1']}
            sobol_stats = pd.DataFrame(sobol_stats, index=self.problem['names'])

            smax = sobol_stats.max().max()
            smin = sobol_stats.min().min()

            # dataframe with s2
            s2 = pd.DataFrame(self.Si['S2'], index=self.problem['names'],
                              columns=self.problem['names'])
            s2[s2 < 0.0] = 0.  # Set negative values to 0 (artifact from small sample sizes)
            s2max = s2.max().max()
            s2min = s2.min().min()

            names = self.problem['names']
            n = len(names)
            ticklocs = np.linspace(0, 2 * pi, n)
            locs = ticklocs[0:-1]

            filtered_names, filtered_locs = filter(self.Si, names, locs,
                                                   criterion, threshold)

            # setup figure
            fig = plt.figure()
            ax = fig.add_subplot(111, polar=True)
            ax.grid(False)
            ax.spines['polar'].set_visible(False)
            ax.set_xticks(ticklocs)

            ax.set_xticklabels(names)
            ax.set_yticklabels([])
            ax.set_ylim(top=1.4)
            legend(ax)

            # plot ST
            plot_circles(ax, filtered_locs, filtered_names, max_s_radius,
                         sobol_stats['ST'], smax, smin, 'w', 'k', 1, 9)

            # plot S1
            plot_circles(ax, filtered_locs, filtered_names, max_s_radius,
                         sobol_stats['S1'], smax, smin, 'k', 'k', 1, 10)

            # plot S2
            for name1, name2 in itertools.combinations(zip(filtered_names, filtered_locs), 2):
                name1, loc1 = name1
                name2, loc2 = name2

                weight = s2.loc[name1, name2]
                lw = 0.5 + max_linewidth_s2 * normalize(weight, s2min, s2max)
                ax.plot([loc1, loc2], [1, 1], c='darkgray', lw=lw, zorder=1)

            out_fn = 'test_radial.png'
            sns.set_style('whitegrid')

        out_path_plot = os.path.join(self.out_path, 'plot_SALib')
        if not os.path.exists(out_path_plot):
            # Create a new directory because it does not exist
            os.makedirs(out_path_plot)
        plt.savefig(os.path.join(out_path_plot, out_fn), dpi=600, format='png')
        # plt.show()


def run_cfes(problem, cfe_instance, param_values, nrun, like_measure, var_measure):
    Y = np.zeros([param_values.shape[0]])
    for i, X in enumerate(param_values):
        print('{} of {}'.format(i+1, nrun))
        Y[i] = salib_cfe_interface(X=X, param_names=problem['names'], myCFE=cfe_instance, like_measure=like_measure, var_measure=var_measure)
    return Y

def run_cfes_for_stability_test(problem, cfe_instance, nrun, like_measure, var_measure):
    S1_estimates = np.zeros([problem['num_vars'], len(nrun)])
    ST_estimates = np.zeros([problem['num_vars'], len(nrun)])
    nsample_record = np.zeros(len(nrun))
    for i in range (len(nrun)):
        print('Major run: n={} of n_Total={}'.format(nrun[i], nrun[-1]))
        sampleset = saltelli.sample(problem, nrun[i], calc_second_order=False)
        Y = np.zeros([sampleset.shape[0]])
        nsample_record[i] = len(sampleset)
        for j, X in enumerate(sampleset):
            print('Subrun: {} of {}'.format(j+1, len(sampleset)))
            Y[j] = salib_cfe_interface(X=X, param_names=problem['names'], myCFE=cfe_instance, like_measure=like_measure, var_measure=var_measure)
        results = sobol.analyze(problem, Y, calc_second_order=False, print_to_console=False)
        ST_estimates[:, i] = results['ST']
        S1_estimates[:, i] = results['S1']
    return S1_estimates, ST_estimates, nsample_record


def salib_cfe_interface(X, param_names, myCFE, like_measure, var_measure):

    # write the randomly-generated parameters to the config json file
    with open(myCFE.cfg_file) as data_file:
        cfe_cfg = json.load(data_file)

    for i, param_name in enumerate(param_names):
        if param_name in ['depth', 'bb', 'mult', 'satdk', 'satpsi', 'slop', 'smcmax', 'wltsmc', 'D']:
            cfe_cfg["soil_params"][param_name] = X[i]
        else:
            cfe_cfg[param_name] = X[i]

    with open(myCFE.cfg_file, 'w') as out_file:
        json.dump(cfe_cfg, out_file)

    # Here the model is actualy started with a unique parameter combination that it gets from spotpy for each time the model is called
    myCFE.initialize()
    myCFE.run_unit_test(plot=False, print_fluxes=False,warm_up=True)

    sim = myCFE.cfe_output_data[["Time", var_measure]]
    sim.loc[:, "Time"] = pd.to_datetime(sim["Time"], format="%Y-%m-%d %H:%M:%S")
    sim = sim.set_index("Time")

    # Get the comparison data
    data = myCFE.unit_test_data
    obs = data[["Time", var_measure]]
    obs.loc[:, "Time"] = pd.to_datetime(obs["Time"], format="%m/%d/%Y %H:%M") # "%d-%b-%Y %H:%M:%S"
    obs = obs.set_index("Time")

    if obs.index[0] != sim.index[0]:
        diff_time = obs.index[0] - sim.index[0]
        warnings.warn("The start of observation and simulation time is different by %s" % diff_time)

    if obs.index[-1] != sim.index[-1]:
        diff_time = obs.index[-1] - sim.index[-1]
        warnings.warn("The end of observation and simulation time is different by %s" % diff_time)

    df = pd.merge_asof(sim, obs, on = "Time")

    sim_synced = df[var_measure+"_x"]
    obs_synced = df[var_measure+"_y"]

    # Calculate objective metrics
    if like_measure == "NashSutcliffe":
        like = spotpy.objectivefunctions.nashsutcliffe(evaluation=obs_synced, simulation=sim_synced)

    return like

def normalize(x, xmin, xmax):
    return (x-xmin)/(xmax-xmin)


def plot_circles(ax, locs, names, max_s, stats, smax, smin, fc, ec, lw,
                 zorder):
    s = np.asarray([stats[name] for name in names])
    s = 0.01 + max_s * np.sqrt(normalize(s, smin, smax))

    fill = True
    for loc, name, si in zip(locs, names, s):
        if fc=='w':
            fill=False
        else:
            ec='none'

        x = np.cos(loc)
        y = np.sin(loc)

        circle = plt.Circle((x,y), radius=si, ec=ec, fc=fc, transform=ax.transData._b,
                            zorder=zorder, lw=lw, fill=True)
        ax.add_artist(circle)


def filter(sobol_indices, names, locs, criterion, threshold):
    if criterion in ['ST', 'S1', 'S2']:
        data = sobol_indices[criterion]
        data = np.abs(data)
        data = data.flatten() # flatten in case of S2
        # TODO:: remove nans

        filtered = ([(name, locs[i]) for i, name in enumerate(names) if
                     data[i]>threshold])
        filtered_names, filtered_locs = zip(*filtered)
    elif criterion in ['ST_conf', 'S1_conf', 'S2_conf']:
        raise NotImplementedError
    else:
        raise ValueError('unknown value for criterion')

    return filtered_names, filtered_locs

from matplotlib.legend_handler import HandlerPatch
class HandlerCircle(HandlerPatch):
    def create_artists(self, legend, orig_handle,
                       xdescent, ydescent, width, height, fontsize, trans):
        center = 0.5 * width - 0.5 * xdescent, 0.5 * height - 0.5 * ydescent
        p = plt.Circle(xy=center, radius=orig_handle.radius)
        self.update_prop(p, orig_handle, legend)
        p.set_transform(trans)
        return [p]

def legend(ax):
    some_identifiers = [plt.Circle((0,0), radius=5, color='k', fill=False, lw=1),
                        plt.Circle((0,0), radius=5, color='k', fill=True),
                        plt.Line2D([0,0.5], [0,0.5], lw=8, color='darkgray')]
    ax.legend(some_identifiers, ['ST', 'S1', 'S2'],
              loc=(1,0.75), borderaxespad=0.1, mode='expand',
              handler_map={plt.Circle: HandlerCircle()})