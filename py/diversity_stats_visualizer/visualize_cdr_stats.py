import os
import sys
import operator

import matplotlib as mplt
mplt.use('Agg')

import utils

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

from Bio.Seq import Seq

def visualize_region_lengths(labeling_df, region, region_name, output_fname, log):
    region_seq = list(labeling_df[region])
    region_len = [len(s) for s in region_seq if len(s) > 1]
    plt.figure()
    plt.hist(region_len)
    plt.xlabel(region_name + ' length (nt)', fontsize = 14)
    plt.ylabel('# ' + region_name + 's', fontsize = 14)
    utils.output_figure(output_fname, region_name + " length distribution", log)

############################################################################
def visualize_length_abundance_dist(labeling_df, region, region_name, output_fname, log):
    region_seq = list(labeling_df[region])
    region_dict = dict()
    for seq in region_seq:
        if seq not in region_dict:
            region_dict[seq] = 0
        region_dict[seq] += 1
    abun = [] #np.array()
    lens = [] #np.array()
    for seq in region_dict:
        if region_dict[seq]  == 1:
            continue
        abun.append(region_dict[seq])
        lens.append(len(seq))
    abun = np.asarray(abun)
    lens = np.asarray(lens)
    f, ax = plt.subplots()
    sns.jointplot(abun, lens, size = 6)
    #plt.xlabel(region_name + ' abundance', fontsize = 14)
    #ax.xaxis.set_label_position('top')
    #plt.ylabel(region_name + ' length', fontsize = 14)
    #plt.xticks(fontsize = 14)
    #plt.yticks(fontsize = 14)
    #plt.xlim(-1, abun.max() + 1)
    #plt.ylim(-1, lens.max() + 1)
    utils.output_figure(output_fname, region_name + " joint distribution of abundances & lengths", log)

############################################################################
def get_region_largest_group(region_seq):
    len_dict = dict()
    for s in region_seq:
        if len(s) not in len_dict:
            len_dict[len(s)] = list()
        len_dict[len(s)].append(s)
    max_len = 0
    max_group = 0
    for l in len_dict:
        if len(len_dict[l]) > max_group:
            max_group = len(len_dict[l])
            max_len = l
    return len_dict[max_len]

def get_nucls_lists(max_group):
    region_len = len(max_group[0])
    nucl_dict = dict()
    nucl_dict['A'] = [0] * region_len
    nucl_dict['C'] = [0] * region_len
    nucl_dict['G'] = [0] * region_len
    nucl_dict['T'] = [0] * region_len
    for s in max_group:
        for i in range(0, len(s)):
            if s[i] in nucl_dict:
                nucl_dict[s[i]][i] += 1
    for i in range(0, region_len):
        sum = nucl_dict['A'][i] + nucl_dict['C'][i] + nucl_dict['G'][i] + nucl_dict['T'][i]
        nucl_dict['A'][i] = float(nucl_dict['A'][i]) / float(sum) * 100
        nucl_dict['C'][i] = float(nucl_dict['C'][i]) / float(sum) * 100
        nucl_dict['G'][i] = float(nucl_dict['G'][i]) / float(sum) * 100
        nucl_dict['T'][i] = float(nucl_dict['T'][i]) / float(sum) * 100
    nucl_dict['A'] = np.array(nucl_dict['A'])
    nucl_dict['C'] = np.array(nucl_dict['C'])
    nucl_dict['G'] = np.array(nucl_dict['G'])
    nucl_dict['T'] = np.array(nucl_dict['T'])
    return nucl_dict

def visualize_largest_region_nucls(labeling_df, region, region_name, output_fname, log):
    region_seq = list(labeling_df[region])
    max_group = get_region_largest_group(region_seq)
    if len(max_group) == 0:
        return
    nucl_dict = get_nucls_lists(max_group)
    x = range(0, len(max_group[0]))
    x_l = [str(i) for i in range(1, len(max_group[0]) + 1)]
    acgt = nucl_dict['A'] + nucl_dict['C'] + nucl_dict['G'] + nucl_dict['T']
    cgt = nucl_dict['C'] + nucl_dict['G'] + nucl_dict['T']
    gt = nucl_dict['G'] + nucl_dict['T']
    #sns.set_color_codes("pastel")
    sns.set_color_codes("muted")
    f, ax = plt.subplots(figsize=(15, 6))
    sns.barplot(x=x, y=acgt, label="A", color = 'b')
    sns.barplot(x=x, y=cgt, label="C", color = 'g')
    sns.barplot(x=x, y=gt, label="G", color = 'r')
    sns.barplot(x=x, y=nucl_dict['T'], label="T", color = 'orange')
    plt.ylim(0, 115)
    ax.legend(ncol = 4, loc="upper center", frameon=True, fontsize = 16)
    plt.xlabel(region_name + ' position (nt)', fontsize = 16)
    plt.ylabel('Nucleotide %', fontsize = 16)
    plt.xticks(x, x_l, fontsize = 14)
    plt.yticks(fontsize = 14)
    utils.output_figure(output_fname, region_name + " nucleotide distribution", log)

############################################################################
amino_acids = ['A', 'G', 'L', 'R', 'W', 'N', 'V', 'I', 'P', 'F', 'Y', 'C', 'T', 'S', 'M', 'Q', 'K', 'H', 'D', 'E', '*']

def get_aa_colors():
    aa_colors = []
    for aa in amino_acids:
        hydrophoby = utils.hydrophoby_dict[aa]
        rel_value = float(hydrophoby - utils.hydro_min) / float(utils.hydro_max - utils.hydro_min)
        aa_colors.append(plt.get_cmap('bwr')(rel_value))
    return aa_colors

def visualize_largest_group_aa_variability(labeling_df, region, region_name, output_fname, log):
    region_seq = list(labeling_df[region])
    max_group = get_region_largest_group(region_seq)
    if len(max_group) == 0:
        return
    group_len = len(max_group[0])
    if group_len % 3 != 0:
        print "Largest " + region_name + " is not out-of-frame"
        return
    aa_len = group_len / 3 
    aa_seqs = [Seq(cdr).translate(to_stop=True) for cdr in max_group]
    aa_dict = {'Position' : [], 'Hidrophobicity' : []}
    for aa_seq in aa_seqs:
        aa_row = [utils.hydrophoby_dict[aa] for aa in aa_seq]
        for i in range(len(aa_row)):
            aa_dict['Position'].append(i + 1)
            aa_dict['Hidrophobicity'].append(aa_row[i])
    plt.figure()
    sns.barplot(x = 'Position', y = 'Hidrophobicity', data = aa_dict, order = range(1, aa_len + 1), color = 'blue')
    plt.xlabel('Position (aa)', fontsize = 14)
    plt.ylabel('Hidrophobicity', fontsize = 14)
    plt.ylim(min(utils.hydrophoby_dict.values()) - 10, max(utils.hydrophoby_dict.values()) + 10)
    utils.output_figure(output_fname, region_name + " aa variability", log)

def output_cdr_stats_for_locus(locus_df, locus_name, column_name, region_name, output_dir, log):
    length_fname = os.path.join(output_dir, locus_name + "_" + region_name + "_length")
    visualize_region_lengths(locus_df, column_name, locus_name + " " + region_name, length_fname, log)
    nucl_fname = os.path.join(output_dir, locus_name + "_" + region_name + "_nucls")
    visualize_largest_region_nucls(locus_df, column_name, locus_name + " " + region_name, nucl_fname, log)
    aa_fname = os.path.join(output_dir, locus_name + "_" + region_name + "_aa")
    visualize_largest_group_aa_variability(locus_df, column_name, locus_name + " " + region_name, aa_fname, log)
    return length_fname, nucl_fname, aa_fname

def output_cdrs_stats_for_locus(vj_df, locus_name, output_config):
    locus_df = vj_df.loc[vj_df['Chain_type'] == locus_name]
    num_records = len(vj_df['Read_name'])
    num_locus_records = len(locus_df['Read_name'])
    if float(num_locus_records) / float(num_records) < .05 or num_locus_records < 10:
        output_config.Log().info("Output contains very low number (" + str(num_locus_records) + ") of " + locus_name + " records. Drawing plots was skipped")
        return
    output_config.Log().info("Visualization of CDR statistics for " + locus_name + " locus")
    region_column_names = ['CDR1_nucls', 'CDR2_nucls', 'CDR3_nucls']
    region_names = ['CDR1', 'CDR2', 'CDR3']
    for col_name, region_name in zip(region_column_names, region_names):
        f1, f2, f3 = output_cdr_stats_for_locus(locus_df, locus_name, col_name, region_name, output_config.cdr_plot_dir, output_config.Log())
        output_config.AddCDRNames(locus_name, region_name, f1, f2, f3)

def main(df_fname, output_config):
    if not os.path.exists(df_fname):
        output_config.Log().info("File containing CDR details " + df_fname + " was not found")
        sys.exit(1)
    vj_df = pd.read_table(df_fname, delim_whitespace = True)
    if len(vj_df['Read_name']) == 0:
        output_config.Log().info("CDR data-frame contains 0 records. CDR visualization will be skipped")
        return
    loci = ['IGH', 'IGK', 'IGL']
    for l in loci:
        output_cdrs_stats_for_locus(vj_df, l, output_config)

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print "Invalid input parameters"
        print "python visualize_cdr_stats.py cdr_details.txt output_config"
        sys.exit(1)
    main(sys.argv[1], sys.argv[2])
