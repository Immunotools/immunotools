import os
import sys
import operator
import warnings

import utils
import plot_output_config

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
            base_v = utils.GetBaseName(self.vj_df['V_hit'][i])
            base_j = utils.GetBaseName(self.vj_df['J_hit'][i])
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
        plt.xticks(range(len(self.sorted_js)), self.sorted_js, rotation = 90, fontsize = 10)
        plt.ylabel('% of sequences')
        utils.output_figure(output_fname, "Usage of J genes", log)

def visualize_vj_stats(labeling_df, output_config):
    vj_matrix = VJMatrix(labeling_df)
    vj_matrix.OutputHeatmap(output_config.vj_usage, output_config.Log())
    vj_matrix.OutputVUsage(output_config.v_usage, output_config.Log())
    vj_matrix.OutputJUsage(output_config.j_usage, output_config.Log())

############################################################################
def main(cdr_details, output_config):
    warnings.filterwarnings('ignore')
    vj_df = pd.read_table(cdr_details, delim_whitespace = True)
    visualize_vj_stats(vj_df, output_config)

if __name__ == "__main__":
    main(sys.argv)
