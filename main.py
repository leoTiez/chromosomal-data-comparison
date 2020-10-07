#!/usr/bin/python3
import numpy as np
import matplotlib.pyplot as plt

import dataHandler as dh


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


def main():
    calc_per_transcript = False
    plot_transcription_diffs = False
    use_pos_test = True
    step_size = 0.01
    explication_tresh = 0.9

    own_bed = dh.load_bam_bed_file(
        name='TSS_TES_steinmetz_jacquier.mRNA.bed',
        rel_path='data/own/reference'
    )
    own_t0_1 = dh.load_big_file(
        name='L3_41_UV3b_CPD_T0.BOWTIE.SacCer3.pe.bin1.RPM.rmdup.bamCoverage.bw',
        rel_path='data/own/GenomeCoverage'
    )
    own_t0_2 = dh.load_big_file(
        name='L3_41_UV3b_CPD_T0.BOWTIE.SacCer3.pe.bin1.RPM.rmdup.bamCoverage.bw',
        rel_path='data/own/GenomeCoverage'
    )
    ref_t0_min = dh.load_big_file(
        name='GSM2109560_UV_0hr_A2_dipy_bkgd_minus.bw',
        rel_path='data/reference_chromosomal_landscape'
    )
    ref_t0_plus = dh.load_big_file(
        name='GSM2109560_UV_0hr_A2_dipy_bkgd_plus.bw',
        rel_path='data/reference_chromosomal_landscape'
    )

    if not use_pos_test:
        gen_mapping, means, stds, all_values = dh.normalise_over_annotation(
            [own_t0_1, ref_t0_min, ref_t0_plus],
            own_bed,
            smoothing=[None, 200, 200]
        )
    else:
        gen_mapping, means, stds, all_values = dh.normalise_over_annotation(
            [own_t0_1, own_t0_2],
            own_bed,
            smoothing=[None, None]
        )

    own_data = all_values[0]
    reference_data = all_values[1] + all_values[2] if not use_pos_test else all_values[1]
    all_diffs_trans = []
    if calc_per_transcript:
        zipped = zip(gen_mapping[0], gen_mapping[1], gen_mapping[2]) \
            if not use_pos_test else zip(gen_mapping[0], gen_mapping[1])
        for own_gen, ref_gen in zipped:
            ref_total = (ref_gen[0] + ref_gen[1]) if not use_pos_test else ref_gen
            diff = own_gen - ref_total
            all_diffs_trans.extend(np.abs(diff).tolist())
            corr = cross_corr(own_gen, ref_total, num_lags=200)

            if plot_transcription_diffs:
                fig, axs = plt.subplots(2, 1)
                fig.tight_layout(h_pad=2.0)
                axs[0].plot(range(own_gen.size), own_gen, label='Own gene')
                # axs[0].plot(range(ref_gen_m.size), ref_gen_m, label='Reference gene -')
                # axs[0].plot(range(ref_gen_p.size), ref_gen_p, label='Reference gene +')
                axs[0].plot(range(ref_total.size), ref_total, label='Reference gene total')
                axs[0].set_xlabel('Position')
                axs[0].set_ylabel('Normalised data value')
                axs[0].legend()

                axs[1].plot(np.arange(-corr.size / 2., corr.size / 2.), corr)
                axs[1].set_title('Cross-correlation as a function of shift $\\tau$')
                axs[1].set_xlabel('$\\tau$')
                axs[1].set_ylabel('Cross-correlation')

                plt.show()

    total_diffs = np.abs(own_data - reference_data)

    ratios_trans = []
    ratios_total = []
    for thresh in np.arange(0, 1., step_size):
        ratios_total.append((np.asarray(total_diffs) < thresh).sum() / float(len(total_diffs)))
        if calc_per_transcript:
            ratios_trans.append((np.asarray(all_diffs_trans) < thresh).sum() / float(len(all_diffs_trans)))

    min_thresh_total = np.arange(0, 1., step_size)[np.asarray(ratios_total) >= explication_tresh][0]
    if calc_per_transcript:
        min_thresh_trans = np.arange(0, 1., step_size)[np.asarray(ratios_trans) >= explication_tresh][0]

    num_subplots = 2 if calc_per_transcript else 1
    fig, axs = plt.subplots(1, num_subplots, figsize=(10, 20))
    ax = axs[0] if calc_per_transcript else axs

    ax.plot(np.arange(0, 1., step_size), ratios_total, label='Differences below the threshold')
    ax.set_xticks([0, min_thresh_total, 1])
    ax.set_yticks([0, explication_tresh, 1])
    ax.plot(
        np.arange(0, 1., step_size),
        np.repeat(explication_tresh, len(ratios_total)),
        label='Minimum for accountability'
    )
    ax.set_xlabel('Threshold')
    ax.set_ylabel('Ratio within the similarity')
    ax.set_title('Ratio of all differences that are closer\nthan the defined threshold')
    ax.get_xticklabels()[1].set_color("red")
    ax.get_yticklabels()[1].set_color("red")

    if calc_per_transcript:
        axs[1].plot(np.arange(0, 1., step_size), ratios_trans, label='Differences below the threshold')
        axs[1].plot(
            np.arange(0, 1., step_size),
            np.repeat(explication_tresh, len(ratios_trans)),
            label='Minimum for accountability'
        )
        axs[1].set_xticks([0, min_thresh_trans, 1])
        axs[1].set_yticks([0, explication_tresh, 1])
        axs[1].set_xlabel('Threshold')
        axs[1].set_ylabel('Ratio within the similarity')
        axs[1].set_title('Ratio of the differences for transcribed genes\nthat are closer than the defined threshold')
        axs[1].get_xticklabels()[1].set_color("red")
        axs[1].get_yticklabels()[1].set_color("red")

    handles, labels = ax.get_legend_handles_labels()
    fig.legend(handles, labels, loc='lower center')
    plt.show()

    mse = cross_mse(own_data, reference_data, num_lags=100)
    plt.plot(np.arange(-mse.size / 2., mse.size / 2.), mse)
    plt.title('MSE as a function of shift $\\tau$')
    plt.xlabel('$\\tau$')
    plt.ylabel('MSE')
    plt.show()

    corr = cross_corr(own_data, reference_data, num_lags=1000)
    plt.plot(np.arange(-corr.size/2., corr.size/2.), corr)
    plt.title('Cross-correlation as a function of shift $\\tau$')
    plt.xlabel('$\\tau$')
    plt.ylabel('Cross-correlation')
    plt.show()


if __name__ == '__main__':
    main()


