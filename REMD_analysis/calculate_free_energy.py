import os
import sys
import time
import glob
import math
import pickle
import pymbar
import natsort
import warnings
import argparse 
import alchemlyb
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib import rc
from matplotlib import cm
from rich.progress import track
from pymbar.timeseries import statistical_inefficiency
from alchemlyb.parsing.gmx import extract_u_nk, extract_dHdl
from alchemlyb.preprocessing.subsampling import equilibrium_detection
from alchemlyb.estimators import TI, BAR, MBAR

def initialize(args):
    parser = argparse.ArgumentParser(
        description='This code performs free energy calculations given dhdl file(s) from GROMACS.')
    parser.add_argument('-d',
                        '--dir',
                        type=str,
                        default='.',
                        help='The directory where the dhdl files are.')
    parser.add_argument('-T',
                        '--temp',
                        type=float,
                        default=298.15,
                        help='The simulation temperature in Kelvin.')
    parser.add_argument('-b',
                        '--bounds',
                        nargs='+',
                        default=[None, None],
                        help='The lower and upper bounds (in ps) for truncating the data.')
    parser.add_argument('-n',
                        '--n_replicas',
                        type=int,
                        required=True,
                        help='The number of replicas. A number of 1 assumes expanded ensemble.')
    parser.add_argument('-s',
                        '--spacing',
                        type=int,
                        default=1,
                        help='The number of steps between samples to used in the timeseries.')

    args_parse = parser.parse_args(args)

    return args_parse 

def preprocess_data(folders, temp, spacing, bounds=[None, None]):
    """
    Preprocesses data with the following steps:
    - For each state, extract the dHdl and u_nk data from the *dhdl.xvg files and concatenate them
    - Truncate the equilibrium region and subsample the dataset with the statistical inefficiency
    - Concatenate dHdl/u_nk data of all replicas

    Parameters
    ----------
    folders : list
        A list of folders where the dhdl files reside.
    temp : float
        The simulation temperature in Kelvin.
    bounds : list
        The lower and uppor bounds for truncating the data

    Returns
    -------
    dHdl_data : pd.Dataframe
        The preprocessed dHdl data that can serve as the input to free energy estimators.
    u_nk_data : pd.Dataframe
        The preprocessed u_nk data that can serve as the input to free energy estimators.
    """
    dHdl_data, u_nk_data = [], []
    for i in range(len(folders)):
        if len(folders) > 1:
            logger(f'\nData preprocessing for state {i}')
        files = glob.glob(os.path.join(folders[i], '*dhdl*xvg*'))
        files = natsort.natsorted(files, reverse=False)

        print(f'  Collecting data for the following files: {", ".join(files)}')
        print('  Subsampling the dHdl and u_nk data ...')
        dHdl = alchemlyb.concat([extract_dHdl(xvg, T=temp) for xvg in files])
        dHdl_series = subsampling.dhdl2series(dHdl)
        dHdl, dHdl_series = subsampling._prepare_input(dHdl, dHdl_series, drop_duplicates=True, sort=True)
        dHdl = subsampling.slicing(dHdl, lower=bounds[0], upper=bounds[1], step=spacing)
        dHdl_series = subsampling.slicing(dHdl_series, lower=bounds[0], upper=bounds[1], step=spacing)
        t, statinef, Neff_max = detect_equilibration(dHdl_series.values)
        print(f'    Adopted spacing: {spacing: .0f}')
        print(f'   {t / len(dHdl_series) * 100: .1f}% of the dHdl data was in the equilibrium region and therfore discarded.')
        print(f'    Statistical inefficiency of dHdl: {statinef: .1f}')
        print(f'    Number of effective samples: {Neff_max}')
        dHdl_series_equil, dHdl_equil = dHdl_series[t:], dHdl[t:]
        indices = subsample_correlated_data(dHdl_series_equil, g=statinef)
        preprocessed_dHdl = dHdl_equil.iloc[indices]

        u_nk = alchemlyb.concat([extract_u_nk(xvg, T=temp) for xvg in files])
        u_nk_series = subsampling.u_nk2series(u_nk)
        u_nk, u_nk_series = subsampling._prepare_input(u_nk, u_nk_series, drop_duplicates=True, sort=True)
        u_nk = subsampling.slicing(u_nk, lower=bounds[0], upper=bounds[1], step=spacing)
        u_nk_series = subsampling.slicing(u_nk_series, lower=bounds[0], upper=bounds[1], step=spacing)
        t, statinef, Neff_max = detect_equilibration(u_nk_series.values)
        print(f'    Adopted spacing: {spacing: .0f}')
        print(f'   {t / len(u_nk_series) * 100: .1f}% of the dHdl data was in the equilibrium region and therfore discarded.')
        print(f'    Statistical inefficiency of dHdl: {statinef: .1f}')
        print(f'    Number of effective samples: {Neff_max}')
        u_nk_series_equil, u_nk_equil = u_nk_series[t:], u_nk[t:]
        indices = subsample_correlated_data(u_nk_series_equil, g=statinef)
        preprocessed_u_nk = u_nk_equil.iloc[indices]

        dHdl_data.append(preprocessed_dHdl)
        u_nk_data.append(preprocessed_u_nk)
    
    # Finally, concatenate dHdl/u_nk data from different replicas as needed
    dHdl_data = alchemlyb.concat(dHdl_data)
    u_nk_data = alchemlyb.concat(u_nk_data)

    logger('Pickling the preprocessed dHdl/u_nk data ...')
    with open('dHdl_u_nk_data.pickle', 'wb') as handle:
        pickle.dump([dHdl, u_nk], handle, protocol=pickle.HIGHEST_PROTOCOL)

    return dHdl_data, u_nk_data

def free_energy_calculation(dHdl, u_nk):
    logger('Fitting TI on dHdl ...')
    ti = TI().fit(dHdl)

    logger('Fitting BAR on u_nk ...')
    bar = BAR().fit(u_nk)

    logger('Fitting MBAR on u_nk ...\n')
    try:
        mbar_stop = False
        mbar = MBAR().fit(u_nk)
    except pymbar.utils.ParameterError():
        mbar_stop = True
        logger("\sum_n W_nk is not equal to 1, probably due to insufficient overlap between states.")
        logger("Stop using MBAR ...")

    logger("------ Results based on the whole dataset ------")
    logger(f"TI: {ti.delta_f_.iloc[0, -1]:.4f} +/- {ti.d_delta_f_.iloc[0, -1]:.4f} kT")
    logger(f"BAR: {bar.delta_f_.iloc[0, -1]:.4f} +/- unknown kT")

    if mbar_stop is False:
        logger(f"MBAR: {mbar.delta_f_.iloc[0, -1]:.4f} +/- {mbar.d_delta_f_.iloc[0, -1]:.4f} kT")
        logger("------------------------------------------------")
        return ti, bar, mbar
    else:
        logger("------------------------------------------------")
        return ti, bar

def get_overlap_matrix(u_nk):
    # sort by state so that rows from same state are in contiguous blocks
    u_nk = u_nk.sort_index(level=u_nk.index.names[1:])

    groups = u_nk.groupby(level=u_nk.index.names[1:])
    N_k = [(len(groups.get_group(i)) if i in groups.groups else 0) for i in u_nk.columns]        

    MBAR = pymbar.mbar.MBAR(u_nk.T, N_k)
    overlap_matrix = np.array(MBAR.computeOverlap()['matrix'])

    return overlap_matrix

def plot_matrix(matrix, png_name, start_idx=0):
    sns.set_context(rc={
    'family': 'sans-serif',
    'sans-serif': ['DejaVu Sans'],
    'size': 5
    })

    K = len(matrix)
    plt.figure(figsize=(K / 3, K / 3))
    annot_matrix = np.zeros([K, K])   # matrix for annotating values

    mask = []
    for i in range(K):
        mask.append([])
        for j in range(len(matrix[0])):
            if matrix[i][j] < 0.005:            
                mask[-1].append(True)
            else:
                mask[-1].append(False)

    for i in range(K):
        for j in range(K):
            annot_matrix[i, j] = round(matrix[i, j], 2)

    x_tick_labels = y_tick_labels = np.arange(start_idx, start_idx + K)
    ax = sns.heatmap(matrix, cmap="YlGnBu", linecolor='silver', linewidth=0.25,
                    annot=annot_matrix, square=True, mask=matrix<0.005, fmt='.2f', cbar=False, xticklabels=x_tick_labels, yticklabels=y_tick_labels)
    ax.xaxis.tick_top()
    ax.tick_params(length=0)
    cmap = cm.get_cmap('YlGnBu')   # to get the facecolor
    ax.set_facecolor(cmap(0))      # use the brightest color (value = 0)
    for _, spine in ax.spines.items():
        spine.set_visible(True)    # add frames to the heat map
    plt.annotate('$\lambda$', xy=(0, 0), xytext=(-0.45, -0.20))
    plt.title('Overlap matrix', fontsize=10, weight='bold')
    plt.tight_layout(pad=1.0)

    plt.savefig(png_name, dpi=600)
    #plt.show()
    plt.close()

def logger(*args, **kwargs):
    print(*args, **kwargs)
    with open("result.txt", "a") as f:
        print(file=f, *args, **kwargs)

def main():
    # Suppress pandas FutureWarning
    warnings.simplefilter(action='ignore', category=FutureWarning)

    rc('font', **{
        'family': 'sans-serif',
        'sans-serif': ['DejaVu Sans'],
        'size': 6
    })
    # Set the font used for MathJax - more on this later
    rc('mathtext', **{'default': 'regular'})
    plt.rc('font', family='serif')

    t1 = time.time()
    args = initialize(sys.argv[1:])

    logger(f'Commandline: {" ".join(sys.argv)}')

    if os.path.isfile('dHdl_u_nk_data.pickle') is True:
        logger('Loading the preprocessed data dHdl and u_nk ...')
        with open('dHdl_u_nk_data.pickle', 'rb') as handle:
            data = pickle.load(handle)
        dHdl, u_nk = data[0], data[1]
    else:
        if args.n_replicas == 1:
            folders = [args.dir]
        else:
            folders = [f'{args.dir}/state_{i}' for i in range(args.n_replicas)]
        dHdl, u_nk = preprocess_data(folders, args.temp, args.spacing, args.bounds)

    logger('\nPerforming free energy calculations on the whole dataset ...')
    output = free_energy_calculation(dHdl, u_nk)

    logger("\nComputing and visualizing the overlap matrix ...")
    matrix = get_overlap_matrix(u_nk)
    plot_matrix(matrix, 'overlap_matrix.png')
    t2 = time.time()
    logger(f"Time elapsed: {t2 - t1:.0f} seconds.")
