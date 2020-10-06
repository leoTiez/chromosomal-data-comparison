#!/usr/bin/python3
import numpy as np
import matplotlib.pyplot as plt

import dataHandler as dh


def main():
    own_bed = dh.load_bam_bed_file(
        name='TSS_TES_steinmetz_jacquier.mRNA.bed',
        rel_path='data/own/reference'
    )
    ref_bed = dh.load_bam_bed_file(
        name='GSM2109559_NoUV_A1_damage.bed',
        rel_path='data/reference_chromosomal_landscape'
    )
    own_t0 = dh.load_big_file(
        name='L3_32_UV4_CPD_T0.BOWTIE.SacCer3.pe.bin1.RPM.bamCoverage.bw',
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

    gen_mapping, means, stds, all_values = dh.normalise_over_annotation(
        [own_t0, ref_t0_min, ref_t0_plus],
        own_bed,
        smoothing=[None, 200, 200]
    )

    # TODO this does not work. Does the NoUV bed file not represent their annotations?
    # TODO the wig files have been manually converted (with an official saccer genome.size file) --> Issue?
    # ref_gen_mapping, ref_means, ref_stds, _ = dl.normalise_over_annotation([ref_t0_min, ref_t0_plus], ref_bed)

    all_diffs_trans = []
    plot_transcription_diffs = False
    for own_gen, ref_gen_m, ref_gen_p in zip(gen_mapping[0], gen_mapping[1], gen_mapping[2]):
        if plot_transcription_diffs:
            plt.plot(range(own_gen.size), own_gen, label='Own gene')
            plt.plot(range(ref_gen_m.size), ref_gen_m, label='Reference gene -')
            plt.plot(range(ref_gen_p.size), ref_gen_p, label='Reference gene +')
            plt.plot(range(ref_gen_p.size), ref_gen_m + ref_gen_p, label='Reference gene total')
            plt.legend()
            plt.show()
        all_diffs_trans.extend(np.abs(own_gen - (ref_gen_m + ref_gen_p)).tolist())

    ratios_trans = []
    ratios_total = []
    total_diffs = np.abs(all_values[0] - (all_values[1] + all_values[2]))
    step_size = 0.01
    for thresh in np.arange(0, 1., step_size):
        ratios_total.append((np.asarray(total_diffs) < thresh).sum() / float(len(total_diffs)))
        ratios_trans.append((np.asarray(all_diffs_trans) < thresh).sum() / float(len(all_diffs_trans)))

    explication_tresh = 0.85
    min_thresh_total = np.arange(0, 1., step_size)[np.asarray(ratios_total) >= explication_tresh][0]
    min_thresh_trans = np.arange(0, 1., step_size)[np.asarray(ratios_trans) >= explication_tresh][0]

    fig, axs = plt.subplots(1, 2, figsize=(10, 20))
    axs[0].plot(np.arange(0, 1., step_size), ratios_trans, label='Differences below the threshold')
    axs[0].plot(
        np.arange(0, 1., step_size),
        np.repeat(explication_tresh, len(ratios_trans)),
        label='Minimum for accountability'
    )
    axs[0].set_xticks([0, min_thresh_trans, 1])
    axs[0].set_xlabel('Threshold')
    axs[0].set_ylabel('Ratio within the similarity')
    axs[0].set_title('Ratio of the differences for transcribed genes\nthat are closer than the defined threshold')
    axs[0].get_xticklabels()[1].set_color("red")

    axs[1].plot(np.arange(0, 1., step_size), ratios_total, label='Differences below the threshold')
    axs[1].set_xticks([0, min_thresh_total, 1])
    axs[1].plot(
        np.arange(0, 1., step_size),
        np.repeat(explication_tresh, len(ratios_total)),
        label='Minimum for accountability'
    )
    axs[1].set_xlabel('Threshold')
    axs[1].set_ylabel('Ratio within the similarity')
    axs[1].set_title('Ratio of all differences that are closer\nthan the defined threshold')
    axs[1].get_xticklabels()[1].set_color("red")

    handles, labels = axs[-1].get_legend_handles_labels()
    fig.legend(handles, labels, loc='lower center')
    plt.show()


if __name__ == '__main__':
    main()


