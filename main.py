#!/usr/bin/python3
import os
import sys
import argparse

from itertools import combinations, product
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt

from datahandler import reader, seqDataHandler


def arg_parse(args):
    parser = argparse.ArgumentParser(description='Pass parameters for comparing 2 or more'
                                                 ' sequencing data sets with each other.')
    parser.add_argument(
        '--input_data',
        '-i',
        action='append',
        required=True,
        help='List with input sequencing signals as a '
             'relative path to your current directory.',
        type=str)
    parser.add_argument('--name', '-n', action='append', type=str, help='Names of the data sets')
    parser.add_argument('--norm', type=str, help='The applied normalisation method. If parameter is not'
                                                 'passed, no normalisation is applied. Otherwise choose'
                                                 'between remap and center. Remap scales data between 0 and 1,'
                                                 'with 1 being the maximum and 0 the minimum. Center shifts'
                                                 'the mean to 0 and sets std to 1.')
    parser.add_argument('--smoothing', action='append', type=int, help='Smoothing values for ')
    parser.add_argument('--thresh', type=float, help='Ratio of data that must not exceed a certain distance distance to'
                                                     'the reference signal that is at least required to treat them'
                                                     'as being essentially the same.')
    parser.add_argument('--save_plot', dest='save_plot', action='store_true', help='Should the plots be saved?')
    parser.add_argument('--save_prefix', type=str, help='Prefix that is added to every saved plot to give them '
                                                        'certain identifiers.')
    parser.add_argument('--num_lags', type=int, help='Maximal number of values that the signal is shifted for the MSE')
    parser.add_argument('--num_bins', type=int, help='Number of bins used for the KS-test. If no parameter is passed,'
                                                     'the data is binned according to the doane method (see '
                                                     'numpy documentation).')

    parsed_args = parser.parse_args(args)
    return parsed_args


def cross_mse(x, y, num_lags):
    mse = np.zeros(2 * num_lags + 1)
    error = (x - y)
    mse[num_lags] = error.dot(error) / float(error.size)
    for lag in range(1, num_lags+1):
        error_p = (x[lag:] - y[:-lag])
        error_m = (x[:-lag] - y[lag:])
        mse[num_lags + lag] = error_p.dot(error_p) / float(error_p.size)
        mse[num_lags - lag] = error_m.dot(error_m) / float(error_m.size)

    return mse


def cross_corr(x, y, num_lags):
    corr = np.zeros(2 * num_lags + 1)
    norm_quot = np.square(np.square(x)).sum() * np.square(np.square(y).sum())
    corr[num_lags] = x.dot(y) / norm_quot
    for lag in range(1, num_lags+1):
        corr[num_lags + lag] = x[lag:].dot(y[:-lag]) / norm_quot
        corr[num_lags - lag] = x[:-lag].dot(y[lag:]) / norm_quot

    return corr


def kolmogorow_smirnow(x, y, binning=10):
    if binning is None:
        _, binning = np.histogram(x, bins='doane')
    if type(binning) == int:
        _, binning = np.histogram(x, bins=binning)

    b_x = np.digitize(x, bins=binning)
    b_y = np.digitize(y, bins=binning)
    cdf_x = stats.norm.cdf(b_x)
    cdf_y = stats.norm.cdf(b_y)

    return stats.ks_2samp(cdf_x, cdf_y)


def main():
    arguments = arg_parse(sys.argv[1:])

    step = 0.01

    num_lags = 100 if arguments.num_lags is None else arguments.num_lags
    thresh = 0.9 if arguments.thresh is None else arguments.thresh
    smoothing = [None if x == 'None' else int(x) for x in arguments.smoothing] \
        if arguments.smoothing is not None else None
    save_plots = False if arguments.save_plot is None else arguments.save_plot
    save_prefix = '' if arguments.save_prefix is None else arguments.save_prefix
    names = list(range(len(arguments.input_data))) if arguments.name is None else arguments.name
    normalisation = arguments.norm
    num_bins = arguments.num_bins

    bw_files = []
    for file_path in arguments.input_data:
        bw_files.append(
            reader.load_big_file(
                name=file_path,
                rel_path=''
            )
        )

    print('########### Applied normalisation method: %s' % normalisation)
    all_values, chrom_start = seqDataHandler.get_values(bw_list=bw_files)
    if normalisation is not None:
        if normalisation.lower() == 'remap':
            print('########### Apply remap normalisation')
            all_values = seqDataHandler.remap_norm_all(all_values)
        elif normalisation.lower() == 'center':
            print('########### Apply center normalisation')
            all_values = seqDataHandler.center_norm_all(all_values)
        else:
            raise ValueError('Invalid normalisation method selected. Use remap or center or do not '
                             'pass anything')
    if smoothing is not None:
        print('########### Apply smoothing')
        all_values = seqDataHandler.smooth_all(all_values, smooth_list=smoothing)

    thresh_list = []
    mse_list = []
    mse_centre_list = []
    ks_list = []

    print('########### Number of lags used: %s' % num_lags)
    print('########### Used threshold: %s' % thresh)
    print('########### Number of bins or binning method: %s' % (num_bins if num_bins is not None else'doane'))
    prod = list(combinations(range(len(all_values)), 2))
    for org, refer in prod:
        diff = np.abs(all_values[org] - all_values[refer])
        thresh_list.append([(np.asarray(diff) < theta).sum() / float(len(diff)) for theta in np.arange(0, 1, step)])
        mse = cross_mse(all_values[org], all_values[refer], num_lags=num_lags)
        mse_list.append(mse)
        d, p = kolmogorow_smirnow(all_values[org], all_values[refer], binning=num_bins)
        ks_list.append(p)
        mse_centre_list.append(mse[num_lags])

    fig_diff, ax_diff = plt.subplots(figsize=(10, 5))
    fig_mse, ax_mse = plt.subplots(figsize=(10, 5))
    for diff_thresh, mse, (org, ref) in zip(thresh_list, mse_list, prod):
        ax_diff.plot(np.arange(0, 1, step), diff_thresh, label='Comparison\n%s, %s' % (names[org], names[ref]))
        ax_mse.plot(np.arange(-len(mse)/2., len(mse)/2.), mse, label='Comparison\n%s, %s' % (names[org], names[ref]))

    ax_diff.plot(
        np.arange(0, 1, step),
        np.repeat(thresh, len(thresh_list[0])),
        label='Minimum similarity'
    )
    ax_diff.set_title('Ratio of the data that differs maximal by $\\theta$')
    ax_diff.set_xlabel('$\\theta$')
    ax_diff.set_ylabel('Ratio')
    box = ax_diff.get_position()
    ax_diff.set_position([box.x0, box.y0, box.width * 0.8, box.height])
    ax_diff.legend(title='Legend', bbox_to_anchor=(1.05, 1))
    fig_diff.tight_layout(rect=[0, 0, 0.85, 1])

    ax_mse.set_title('MSE as a function of shift $\\tau$')
    ax_mse.set_xlabel('$\\tau$')
    ax_mse.set_ylabel('MSE')
    box = ax_mse.get_position()
    ax_mse.set_position([box.x0, box.y0, box.width * 0.8, box.height])
    ax_mse.legend(title='Legend', bbox_to_anchor=(1.05, 1))
    fig_mse.tight_layout(rect=[0, 0, 0.85, 1])

    fig_heat, ax_heat = plt.subplots()
    heatmap = np.zeros((len(all_values), len(all_values)))
    mask_tri = np.triu(np.ones(heatmap.shape, dtype=bool))
    heatmap[np.logical_and(mask_tri, ~np.eye(heatmap.shape[0], dtype=bool))] = np.asarray(ks_list)
    heatmap.T[np.logical_and(mask_tri, ~np.eye(heatmap.shape[0], dtype=bool))] = np.asarray(ks_list)
    heat_plt = ax_heat.imshow(heatmap, cmap='seismic')
    ax_heat.set_xticks(np.arange(len(names)))
    ax_heat.set_yticks(np.arange(len(names)))
    ax_heat.set_xticklabels(names)
    ax_heat.set_yticklabels(names)
    fig_heat.colorbar(heat_plt)
    plt.setp(ax_heat.get_xticklabels(), rotation=45, ha='right', rotation_mode='anchor')

    for r, c in product(range(heatmap.shape[0]), range(heatmap.shape[1])):
        text = ax_heat.text(r, c, '%.2f' % heatmap[r, c], ha='center', va='center', color='w')

    ax_heat.set_title('KS p-value Heatmap')
    fig_heat.tight_layout()

    if not save_plots:
        plt.show()
    else:
        curr_dir = os.getcwd()
        fig_diff.savefig('%s/%s_diff.png' % (curr_dir, save_prefix))
        fig_mse.savefig('%s/%s_mse.png' % (curr_dir, save_prefix))
        fig_heat.savefig('%s/%s_heat.png' % (curr_dir, save_prefix))
        plt.clf()


if __name__ == '__main__':
    main()
