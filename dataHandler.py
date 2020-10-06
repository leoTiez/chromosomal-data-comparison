#!/usr/bin/python3
import os

import roman
import numpy as np
import pyBigWig
from pybedtools import BedTool


def load_bam_bed_file(name, rel_path='data'):
    cur_dir = os.getcwd()
    abs_path = cur_dir + '/' + rel_path + '/' + name
    file = BedTool(abs_path)
    return file


def load_big_file(name, rel_path='data'):
    cur_dir = os.getcwd()
    abs_path = cur_dir + '/' + rel_path + '/' + name

    file = pyBigWig.open(abs_path)
    return file


def normalise_over_annotation(bw_list, bed_ref, smoothing = [None, None, None]):
    if len(bw_list) == 0:
        raise ValueError('List with bigwig objects must not be empty.')

    all_values = []
    bw_gen_mapping = []
    for _ in bw_list:
        all_values.append([])
        bw_gen_mapping.append([])

    counter = 0
    chrom_start = {}
    for chrom_num, _ in enumerate(bw_list[0].chroms().values()):
        roman_num = roman.toRoman(chrom_num+1)
        if chrom_num == 16:
            roman_num = 'M'
        length = int(bw_list[0].chroms('chr%s' % roman_num))
        chrom_start['chr%s' % roman_num] = counter
        counter += length
        for num, bw in enumerate(bw_list):
            values = np.asarray(bw.values('chr%s' % roman_num, 0, length))
            all_values[num].extend(np.nan_to_num(values, nan=0.0).tolist())

    means = []
    stds = []
    for num, smooth in enumerate(smoothing):
        all_values[num] = np.asarray(all_values[num])
        if smooth is not None:
            all_values[num] = np.convolve(all_values[num], np.ones(smooth), mode='same')
        means.append(all_values[num].mean())
        stds.append(all_values[num].std())
        all_values[num] = (all_values[num] - means[-1]) / stds[-1]

    for int_num, interval in enumerate(bed_ref):
        for num, (values, gen_mapping) in enumerate(zip(all_values, bw_gen_mapping)):
            start = chrom_start[interval[0]]
            values = values[start + int(interval[1]): start + int(interval[2])]
            values = np.nan_to_num(values, copy=False, nan=0.)
            gen_mapping.append(values)

    return bw_gen_mapping, means, stds, all_values




