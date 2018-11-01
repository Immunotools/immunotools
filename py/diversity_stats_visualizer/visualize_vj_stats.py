import os
import sys
import operator
import warnings

import utils

import matplotlib as mplt
mplt.use('Agg')

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def GetBaseName(gene_name):
    return gene_name.split('*')[0]

class VJMatrix:
    def __init__(self, vj_df):
        self.vj_df = vj_df
        self._CreateVJDicts()
        self._CreateVJMatrix()

    def _CreateVJDicts(self):
        self.vj_dict = dict()
        self.v_dict = dict()
        self.j_dict = dict()
        for i in range(len(self.vj_df)):
            base_v = GetBaseName(self.vj_df['V_hit'][i])
            base_j = GetBaseName(self.vj_df['J_hit'][i])
            if (base_v, base_j) not in self.vj_dict:
                self.vj_dict[(base_v, base_j)] = 0
            self.vj_dict[(base_v, base_j)] += 1
            if base_v not in self.v_dict:
                self.v_dict[base_v] = 0
            self.v_dict[base_v] += 1
            if base_j not in self.j_dict:
                self.j_dict[base_j] = 0
            self.j_dict[base_j] += 1
        self.sorted_vs = sorted(self.v_dict.keys())
        self.sorted_js = sorted(self.j_dict.keys())

    def _FilterLowAbundantVGenes(self):
        self.used_vs = []
        self.abundant_vj_matrix = []
        for i in range(len(self.vj_matrix)):
            sum_percentages = sum(self.vj_matrix[i])
            if sum_percentages < 0.1:
                continue
            self.abundant_vj_matrix.append(self.vj_matrix[i])
            self.used_vs.append(self.sorted_vs[i])

    def _CreateVJMatrix(self):
        self.vj_matrix = []
        for v in self.sorted_vs:
            self.vj_matrix.append([0] * len(self.sorted_js))
        for vj in self.vj_dict:
            v_index = self.sorted_vs.index(vj[0])
            j_index = self.sorted_js.index(vj[1])
            self.vj_matrix[v_index][j_index] = float(self.vj_dict[vj]) / len(self.vj_df) * 100
        self._FilterLowAbundantVGenes()

    def OutputHeatmap(self, output_fname, log):
        plt.figure(figsize = (10, 15))
        sns.heatmap(np.array(self.abundant_vj_matrix), xticklabels = self.sorted_js, yticklabels = self.used_vs, cmap = 'jet')
        plt.yticks(rotation = 0, fontsize = 10)
        plt.xticks(rotation = 90, fontsize = 10)
        utils.output_figure(output_fname, "VJ heatmap for the most abundant VJ combinations", log)

    def OutputVUsage(self, output_fname, log):
        plt.figure(figsize = (10, 8))
        perc_list = [float(self.v_dict[v]) / len(self.vj_df) * 100 for v in self.sorted_vs]
        plt.bar(range(len(self.sorted_vs)), perc_list)
        plt.xticks(range(len(self.sorted_vs)), self.sorted_vs, rotation = 90, fontsize = 10)
        plt.ylabel('% of sequences')
        utils.output_figure(output_fname, "Usage of V genes", log)

    def OutputJUsage(self, output_fname, log):
        plt.figure(figsize = (10, 8))
        perc_list = [float(self.j_dict[j]) / len(self.vj_df) * 100 for j in self.sorted_js]
        plt.bar(range(len(self.sorted_js)), perc_list)
        plt.xticks(range(len(self.sorted_js)), self.sorted_vs, rotation = 90, fontsize = 10)
        plt.ylabel('% of sequences')
        utils.output_figure(output_fname, "Usage of J genes", log)

def visualize_vj_stats(labeling_df, output_dir, log):
    vj_matrix = VJMatrix(labeling_df)
    vj_matrix.OutputHeatmap(os.path.join(output_dir, 'vj_heatmap'), log)
    vj_matrix.OutputVUsage(os.path.join(output_dir, 'v_usage'), log)
    vj_matrix.OutputJUsage(os.path.join(output_dir, 'j_usage'), log)

############################################################################
def checkout_output_dir_fatal(output_dir, log):
    if not os.path.exists(output_dir):
        log.info("ERROR: Directory " + output_dir + " was not found")
        sys.exit(1)

def checkout_output_dir(output_dir):
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)

def main(argv):
    warnings.filterwarnings('ignore')
    if len(argv) != 5:
        print "Invalid input parameters"
        print "python visualize_vj_stats.py cdr_details.txt output_dir plot_dir logger"
        return
    vj_df = pd.read_table(argv[1], delim_whitespace = True)
    output_dir = argv[2]
    plot_dir = argv[3]
    log = argv[4] #utils.get_logger_by_arg(plot_dir, "diversity_analyzer_vis")
    checkout_output_dir_fatal(output_dir, log)
    checkout_output_dir(plot_dir)
    visualize_vj_stats(vj_df, plot_dir, log)

if __name__ == "__main__":
    main(sys.argv)
