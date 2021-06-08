#!/usr/bin/python3
"""Chromosomal Data Comparison

This script provides a command line interface to compare sequencing data. It combines a variety of different tools,
inter alia mean-squared with a lag function, cross correlation, and Kolmogorov-Smirnow.
Run
```
python3 main.py --help
```
to print all the command line options.
"""
import os
import sys
import argparse
import multiprocessing

from itertools import combinations, product
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt

from datahandler import reader, seqDataHandler


def print_status(per):
    """
    Print process bar
    :param per: Percent
    :type per: float
    :return: None
    """
    sys.stdout.write('\r')
    sys.stdout.write("[%-100s] %d%%" % ('=' * int(per), per))
    sys.stdout.flush()


def arg_parse(args):
    """
    Parse command line arguments
    :param args: Argument list
    :type args: list
    :return: Parsed argiments
    """
    parser = argparse.ArgumentParser(description='Pass parameters for comparing 2 or more'
                                                 ' sequencing data sets with each other.')
    parser.add_argument('--input_data', '-i', action='append', required=True, type=str,
                        help='List with input sequencing signals as a relative path to your current directory.')
    parser.add_argument('--name', '-n', action='append', type=str,
                        help='Names of the data sets')
    parser.add_argument('--norm', type=str,
                        help='The applied normalisation method. If parameter is not passed, no normalisation is'
                             'applied. Otherwise choose between remap and center. Remap scales data between 0 and 1,'
                             'with 1 being the maximum and 0 the minimum. Center shifts the mean to 0 and sets'
                             'std to 1.')
    parser.add_argument('--bed', type=str,
                        help='If passed, data is partitioned according to bed file definitions.')
    parser.add_argument('--smoothing', action='append', type=int,
                        help='Smoothing values')
    parser.add_argument('--thresh', type=float,
                        help='Ratio of data that must not exceed a certain distance distance to the reference signal '
                             'that is at least required to treat them as being essentially the same.')
    parser.add_argument('--save_plot', dest='save_plot', action='store_true',
                        help='Should the plots be saved?')
    parser.add_argument('--save_prefix', type=str,
                        help='Prefix that is added to every saved plot to give them certain identifiers.')
    parser.add_argument('--num_lags', type=int,
                        help='Maximal number of values that the signal is shifted for the MSE')
    parser.add_argument('--num_bins', type=int,
                        help='Number of bins used for the KS-test. If no parameter is passed, the data is binned'
                             'according to the doane method (see numpy documentation).')

    parsed_args = parser.parse_args(args)
    return parsed_args


def cross_mse(x, y, num_lags):
    """
    Mean-squared error with a lag function
    :param x: Data array 1
    :type x: numpy.array
    :param y: Data array 2
    :type y: numpy.array
    :param num_lags: Number of lags to the left and right that are to be plotted
    :type num_lags: int
    :return: Array with the mean-squared error sorted from the smallest lag (i.e. -num_lags) to the larges lag
    (i.e. num_lags)
    """
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
    """
    Calculate the linear cross correlation
    :param x: Data array 1
    :type x: numpy.array
    :param y: Data array 2
    :type y: numpy.array
    :param num_lags: Number of lags to the left and right that are to be plotted
    :type num_lags: int
    :return: Array with the linear cross correlation sorted from the smallest lag (i.e. -num_lags) to the larges lag
    (i.e. num_lags)
    """
    corr = np.zeros(2 * num_lags + 1)
    norm_quot = np.square(np.square(x)).sum() * np.square(np.square(y).sum())
    corr[num_lags] = x.dot(y) / norm_quot
    for lag in range(1, num_lags+1):
        corr[num_lags + lag] = x[lag:].dot(y[:-lag]) / norm_quot
        corr[num_lags - lag] = x[:-lag].dot(y[lag:]) / norm_quot

    return corr


def kolmogorow_smirnow(x, y, binning=10):
    """
    Calculate Kolmogorow-Smirnow test
    :param x: Data array 1
    :type x: numpy.array
    :param y: Data array 2
    :type y: numpy.array
    :param binning: Number of bins
    :type binning: int or str
    :return: Largest distance between the cumulative distribution functions, p-value
    """
    if binning is None:
        _, binning = np.histogram(x, bins='doane')
    if isinstance(binning, int) or isinstance(binning, np.integer):
        _, binning = np.histogram(x, bins=binning)

    b_x = np.digitize(x, bins=binning)
    b_y = np.digitize(y, bins=binning)
    cdf_x = stats.norm.cdf(b_x)
    cdf_y = stats.norm.cdf(b_y)

    return stats.ks_2samp(cdf_x, cdf_y)


def parallel_diff(x, y):
    """
    Wrapper function for running absolute distance calculation in parallel
    :param x: Data array 1
    :type x: numpy.array
    :param y: Data array 2
    :type y: numpy.array
    :return: Array with absolute distance
    """
    return np.abs(x - y)


def parallel_thresh(x, theta):
    """
    Wrapper function for running thresholding in parallel. Determine the ratio of values that is below the threshold.
    :param x: Data array
    :type x: numpy.array
    :param theta: threshold value
    :type theta: float
    :return: Ratio of values lower than the ratio
    """
    return (np.asarray(x) < theta).sum() / float(len(x))


def parallel_mse(x, y, num_lags):
    """
    Wrapper function for running mean-squared error calculations in parallel
    :param x: Data array 1
    :type x: numpy.array
    :param y: Data array 2
    :type y: numpy.array
    :param num_lags: Number of shifts to the left and the right
    :type num_lags: int
    :return: Mean-squared error with lag function (see docstring  of cross_mse)
    """
    return cross_mse(x, y, num_lags=num_lags)


def parallel_ks(x, y, num_bins):
    """
    Wrapper function for running KS-Test in parallel
    :param x: Data array 1
    :type x: numpy.array
    :param y: Data array 2
    :type y: numpy.array
    :param num_bins: Number of bins
    :type num_bins: int or str
    :return: p-value
    """
    return kolmogorow_smirnow(np.asarray(x), np.asarray(y), binning=num_bins)[1]


def main():
    """
    Main function which is executed through the command line interface
    :return: None
    """
    arguments = arg_parse(sys.argv[1:])

    step = 0.01

    num_lags = 100 if arguments.num_lags is None else arguments.num_lags
    thresh = 0.9 if arguments.thresh is None else arguments.thresh
    smoothing = [None if x == 'None' else int(x) for x in arguments.smoothing] \
        if arguments.smoothing is not None else None
    save_plots = False if arguments.save_plot is None else arguments.save_plot
    save_prefix = '' if arguments.save_prefix is None else arguments.save_prefix
    names = list(range(len(arguments.input_data))) if arguments.name is None else arguments.name
    bed = arguments.bed
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

    if bed:
        bed = reader.load_bam_bed_file(bed, rel_path='', is_abs_path=True)

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

    if bed:
        print('########### Segment data according to bed file')
        all_values, _ = seqDataHandler.annotate_all(all_values, bed_ref=bed, chrom_start=chrom_start)
        max_len = len(sorted(all_values[0], reverse=True, key=len)[0])
        for i in range(len(all_values[0])):
            len_diff = max_len - len(all_values[0][i])
            rp = len_diff // 2
            lp = len_diff - rp
            for sig in range(len(all_values)):
                all_values[sig][i] = np.pad(all_values[sig][i], (lp, rp), 'constant', constant_values=0)
    else:
        all_values = [all_values[i][np.newaxis, ...] for i in range(len(all_values))]

    thresh_list = []
    mse_list = []
    mse_centre_list = []
    ks_list = []

    print('########### Number of lags used: %s' % num_lags)
    print('########### Used threshold: %s' % thresh)
    print('########### Number of bins or binning method: %s' % (num_bins if num_bins is not None else'doane'))
    print('########### Calculate similarity')
    prod = list(combinations(range(len(all_values)), 2))
    thresh_steps = np.arange(0, 1, step)
    with multiprocessing.Pool(np.maximum(multiprocessing.cpu_count() - 1, 0)) as parallel:
        for num, (org, refer) in enumerate(prod):
            print_status(100 * (num / float(len(prod))))
            diff = parallel.starmap(parallel_diff, zip(all_values[org], all_values[refer]))
            thresh_pair = []
            for d in diff:
                thresh_pair.append(parallel.starmap(
                    parallel_thresh,
                    zip(np.repeat(d, thresh_steps.size).reshape((d.size, thresh_steps.size)), thresh_steps)
                ))
            thresh_list.append(thresh_pair)
            mse_pair = parallel.starmap(
                parallel_mse,
                zip(all_values[org], all_values[refer], np.repeat(num_lags, len(all_values[org])))
            )
            mse_list.append(mse_pair)
            ks_pair = parallel.starmap(
                parallel_ks,
                zip(all_values[org], all_values[refer], np.repeat(num_bins, len(all_values[org])))
            )
            ks_list.append(ks_pair)
            mse_centre_list.append([mse_pair[i][num_lags] for i in range(len(all_values[org]))])

    print_status(100)
    fig_diff, ax_diff = plt.subplots(figsize=(10, 5))
    fig_mse, ax_mse = plt.subplots(figsize=(10, 5))

    for diff_thresh, mse, (org, ref) in zip(thresh_list, mse_list, prod):
        line = ax_diff.plot(
            thresh_steps,
            np.mean(diff_thresh, axis=0),
            label='Comparison\n%s, %s' % (names[org], names[ref])
        )
        ax_diff.fill_between(
            thresh_steps,
            np.maximum(np.mean(diff_thresh, axis=0) - np.std(diff_thresh, axis=0), 0),
            np.mean(diff_thresh, axis=0) + np.std(diff_thresh, axis=0),
            color=line[0].get_color(),
            alpha=.2
        )

        ax_mse.plot(
            np.arange(-len(mse[0])/2., len(mse[0])/2.),
            np.mean(mse, axis=0),
            label='Comparison\n%s, %s' % (names[org], names[ref]),
            color=line[0].get_color()
        )
        ax_mse.fill_between(
            np.arange(-len(mse[0]) / 2., len(mse[0]) / 2.),
            np.maximum(np.mean(mse, axis=0) - np.std(mse, axis=0), 0),
            np.mean(mse, axis=0) + np.std(mse, axis=0),
            color=line[0].get_color(),
            alpha=.2
        )

    ax_diff.plot(
        np.arange(0, 1, step),
        np.repeat(thresh, np.arange(0, 1, step).size),
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
    heatmap_idx = np.zeros((len(all_values), len(all_values)))
    mask_tri = np.triu(np.ones(heatmap.shape, dtype=bool))
    heatmap[np.logical_and(mask_tri, ~np.eye(heatmap.shape[0], dtype=bool))] = np.mean(ks_list, axis=1)
    heatmap_idx[np.logical_and(mask_tri, ~np.eye(heatmap.shape[0], dtype=bool))] = np.arange(len(ks_list))
    heatmap.T[np.logical_and(mask_tri, ~np.eye(heatmap.shape[0], dtype=bool))] = np.mean(ks_list, axis=1)
    heatmap_idx.T[np.logical_and(mask_tri, ~np.eye(heatmap.shape[0], dtype=bool))] = np.arange(len(ks_list))
    heat_plt = ax_heat.imshow(heatmap, cmap='seismic')
    ax_heat.set_xticks(np.arange(len(names)))
    ax_heat.set_yticks(np.arange(len(names)))
    ax_heat.set_xticklabels(names)
    ax_heat.set_yticklabels(names)
    fig_heat.colorbar(heat_plt)
    plt.setp(ax_heat.get_xticklabels(), rotation=45, ha='right', rotation_mode='anchor')

    for r, c in product(range(heatmap.shape[0]), range(heatmap.shape[1])):
        ks_values = ks_list[int(heatmap_idx[r, c])]
        text = '%.2f\n(+-%.2f)' % (heatmap[r, c], 0.0 if r == c else np.std(ks_values))
        ax_heat.text(r, c, text, ha='center', va='center', color='w')

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
