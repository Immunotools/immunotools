import os
import sys
import operator

import igraph 
from igraph import *

import matplotlib as mplt
mplt.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.patches as patches
import seaborn as sns
sns.set()

from Bio import SeqIO
from Bio.Seq import Seq

import utils
import dataset

class AminoAcidDict:
    def __init__(self):
        self.aa_dict = dict() # aa seq -> list of seq IDs
        self.id_aa_dict = dict() # seq D -> aa seq
        self.aa_index_map = dict() # aa seq -> index in aa_list
        self.aa_list = [] # refactor: too many storages of AA seqs

    def AddLineage(self, full_length_lineage):
        self.full_length_lineage = full_length_lineage
        self._InitAADictFromLineage()

    def _InitAADictFromLineage(self):
        for seq in self.full_length_lineage.FullLengthSeqIdIter():
            aa_seq = str(Seq(seq.seq).translate())
            self._AddAminoAcidSeq(seq, aa_seq)

    def AddClonalTree(self, clonal_tree):
        self.clonal_tree = clonal_tree
        self.full_length_lineage = clonal_tree.FullLengthLineage()
        self._InitAADictFromTree()

    def _InitAADictFromTree(self):
        for v in self.clonal_tree.VertexIter():
            seq = self.clonal_tree.GetSequenceByVertex(v)
            aa_seq = str(Seq(seq.seq).translate())
            self._AddAminoAcidSeq(seq, aa_seq)

    def _AddAminoAcidSeq(self, nucl_seq, aa_seq):
        self.id_aa_dict[nucl_seq.id] = aa_seq
        if aa_seq not in self.aa_dict:
            self.aa_dict[aa_seq] = []
        self.aa_dict[aa_seq].append(nucl_seq.id)
        if aa_seq not in self.aa_index_map:
            self.aa_list.append(aa_seq)
            self.aa_index_map[aa_seq] = len(self.aa_list) - 1 

    def __iter__(self):
        for aa in self.aa_dict:
            yield aa

    def GetAAMultiplicity(self, aa):
        return len(self.aa_dict[aa])

    def GetIdsByAA(self, aa):
        return self.aa_dict[aa]

    def GetAAById(self, seq_id):
        return self.id_aa_dict[seq_id]

    def GetIndexByAA(self, aa):
        return self.aa_index_map[aa]

    def GetAAByIndex(self, aa_index):
        return self.aa_list[aa_index]

class AADiversityAnalyzer:
    def __init__(self, full_length_lineage):
        self.full_length_lineage = full_length_lineage
        self.aa_dict = AminoAcidDict(self.full_length_lineage)
        self._FindMostAbundantAALength()
    
    def _FindMostAbundantAALength(self):
        len_dict = dict()
        for aa in self.aa_dict:
            if len(aa) not in len_dict:
                len_dict[len(aa)] = 0
            len_dict[len(aa)] += 1
        self.most_freq_len = max(len_dict.iteritems(), key=operator.itemgetter(1))[0]

    def _AddCDRsOnAminoAcidPlot(self, ax, y_width):
        seq_id = -1
        tree = self.full_length_lineage.UndirectedClonalTree()
        for v in tree.VertexIter():
            v_aa_seq = self.aa_dict.GetAAById(self.full_length_lineage.GetFullLengthSequenceByIndex(v).id)
            if len(v_aa_seq) == self.most_freq_len:
                seq_id = self.full_length_lineage.GetFullLengthSequenceByIndex(v).id
                break
        cdr1_bounds = self.full_length_lineage.Dataset().GetCDR1BoundsBySeqName(seq_id)
        cdr1_bounds = (cdr1_bounds[0] / 3, cdr1_bounds[1] / 3)
        cdr2_bounds = self.full_length_lineage.Dataset().GetCDR2BoundsBySeqName(seq_id)
        cdr2_bounds = (cdr2_bounds[0] / 3, cdr2_bounds[1] / 3)
        cdr3_bounds = self.full_length_lineage.Dataset().GetCDR3BoundsBySeqName(seq_id)
        cdr3_bounds = (cdr3_bounds[0] / 3, cdr3_bounds[1] / 3)
        rect_1 = patches.Rectangle((cdr1_bounds[0], 0), cdr1_bounds[1] - cdr1_bounds[0], y_width, facecolor = 'red', alpha = 0.25)
        ax.add_patch(rect_1)
        rect_2 = patches.Rectangle((cdr2_bounds[0], 0), cdr2_bounds[1] - cdr2_bounds[0], y_width, facecolor = 'red', alpha = 0.25)
        ax.add_patch(rect_2)
        rect_3 = patches.Rectangle((cdr3_bounds[0], 0), cdr3_bounds[1] - cdr3_bounds[0], y_width, facecolor = 'red', alpha = 0.25)
        ax.add_patch(rect_3)

    def ComputeVariabilityPlot(self, output_fname):
        # for each aa position compute how many aa in it
        self.variability_sets = []
        for i in range(self.most_freq_len):
            self.variability_sets.append(set())
        num_used_seq = 0
        num_good_seq = 0
        for aa in self.aa_dict:
            if len(aa) != self.most_freq_len:
                continue
            num_good_seq += 1
            for i in range(self.most_freq_len):
                self.variability_sets[i].add(aa[i])
            num_used_seq += len(self.aa_dict[aa])
        self.variability = [len(self.variability_sets[i]) for i in range(self.most_freq_len)]
        fig, ax = plt.subplots(1)
        self._AddCDRsOnAminoAcidPlot(ax, max(self.variability) + 0.25)
        plt.bar(range(self.most_freq_len), self.variability)
        plt.ylim(1, max(self.variability) + 0.25)
        plt.xlabel('Amino-acid position')
        plt.ylabel('# amino-acids')
        plt.title(str(num_used_seq) + ' aa seq were used')
        utils.OutputPlotToPdf(output_fname)

    def OutputPublicAA(self, output_fname):
        tree = self.full_length_lineage.UndirectedClonalTree()
        frac_shared_aa = []
        for v in tree.VertexIter():
            if tree.VertexDegree(v) < 10:
                continue
            neighs = tree.GetVertexNeighs(v)
            v_aa_seq = self.aa_dict.GetAAById(self.full_length_lineage.GetFullLengthSequenceByIndex(v).id)
            num_shared_aa = 0
            for n in neighs:
                n_aa_seq = self.aa_dict.GetAAById(self.full_length_lineage.GetFullLengthSequenceByIndex(n).id)
                if v_aa_seq == n_aa_seq:
                    num_shared_aa += 1
            frac_shared_aa.append((float(num_shared_aa) / tree.VertexDegree(v), tree.VertexDegree(v)))
        return frac_shared_aa
        print("Fraction of shared aa seqs: " + str(frac_shared_aa))

    def OutputNumCodonsPerAAPosition(self, output_fname):
        tree = self.full_length_lineage.UndirectedClonalTree()
        most_frequent_aa = [''] * self.most_freq_len
        for i in range(self.most_freq_len):
            cur_aa_dict = dict()
            for v in tree.VertexIter():
                v_aa_seq = self.aa_dict.GetAAById(self.full_length_lineage.GetFullLengthSequenceByIndex(v).id)
                if len(v_aa_seq) != self.most_freq_len:
                    continue
                if v_aa_seq[i] not in cur_aa_dict:
                    cur_aa_dict[v_aa_seq[i]] = 0
                cur_aa_dict[v_aa_seq[i]] += 1
            most_frequent_aa[i] = max(cur_aa_dict.iteritems(), key=operator.itemgetter(1))[0]
        codons = [0] * self.most_freq_len
        for i in range(self.most_freq_len):
            cur_codons = set()
            for v in tree.VertexIter():
                v_nucl_seq = self.full_length_lineage.GetFullLengthSequenceByIndex(v).seq
                v_aa_seq = self.aa_dict.GetAAById(self.full_length_lineage.GetFullLengthSequenceByIndex(v).id)
                if len(v_aa_seq) != self.most_freq_len:
                    continue
                if v_aa_seq[i] != most_frequent_aa[i]:
                    continue
                cur_codons.add(v_nucl_seq[i * 3 : i * 3 + 3])
            codons[i] = len(cur_codons)
        #print codons
        fig, ax = plt.subplots(1)
        self._AddCDRsOnAminoAcidPlot(ax, max(codons) + 0.25)
        plt.bar(range(self.most_freq_len), codons)
        plt.ylim(1, max(codons) + 0.25)
        utils.OutputPlotToPdf(output_fname)

class AminoAcidTreeConstructor:
    def __init__(self, min_aa_mult = 20):
        self.trees = []
        self.weights = []
        self.lineage_ids = []
        self.min_aa_mult = min_aa_mult

    def AddLineage(self, full_length_lineage):
        aa_dict = AminoAcidDict(full_length_lineage)
        for aa in aa_dict:
            if aa_dict.GetAAMultiplicity(aa) < self.min_aa_mult:
                continue
            print("Processing aa with multiplicity " + str(aa_dict.GetAAMultiplicity(aa)) + '...') 
            seq_ids = aa_dict.GetIdsByAA(aa)
            print("V genes: " + str(set([full_length_lineage.Dataset().GetGeneHitBySeqName(seq_id, dataset.AnnotatedGene.V) for seq_id in seq_ids])))
            print("J genes: " + str(set([full_length_lineage.Dataset().GetGeneHitBySeqName(seq_id, dataset.AnnotatedGene.J) for seq_id in seq_ids])))
            seqs = [full_length_lineage.GetFullLengthSequenceByName(seq_id).seq for seq_id in seq_ids]
            hgraph, distance_dict = self._ConstructGraph(seqs)
            distances = []
            for i in range(len(seqs)):
                for j in range(i + 1, len(seqs)):
                    distances.append(distance_dict[(i, j)])
            spanning_tree = hgraph.spanning_tree(weights = distances, return_tree = True)
            self.trees.append(spanning_tree)
            tree_weights = [distance_dict[e] for e in spanning_tree.get_edgelist()]
            self.weights.append(tree_weights)
            self.lineage_ids.append(full_length_lineage.id())
            
    def _ConstructGraph(self, seqs):
        graph = Graph()
        graph.add_vertices(len(seqs))
        edges = []
        edge_weights = dict()
        for i in range(len(seqs)):
            for j in range(i + 1, len(seqs)):
                edges.append((i, j))
                edge_weights[(i, j)] = utils.HammingDistance(seqs[i], seqs[j])
        graph.add_edges(edges)
        return graph, edge_weights

    def OutputTrees(self, output_dir):
        for i in range(len(self.trees)):
            self._OutputTreeByIndex(i, os.path.join(output_dir, 'lineage_' + self.lineage_ids[i] + '_tree_' + str(i) + '_size_' + str(len(self.weights[i]) - 1)))

    def _OutputTreeByIndex(self, tree_ind, output_base):
        output_dot = output_base + '.dot'
        tree_edges = self.trees[tree_ind].get_edgelist()
        weights = self.weights[tree_ind]
        fh = open(output_dot, 'w')
        fh.write('digraph{\n')
        for i in range(len(tree_edges)):
            fh.write(str(tree_edges[i][0]) + ' -> ' + str(tree_edges[i][1]) + ' [dir = none, tooltip = ' + str(weights[i]) + ', penwidth = ' + str(weights[i]) + ']\n')
        fh.write('}\n')
        fh.close()
        os.system('fdp -Tsvg ' + output_dot + ' -o ' + output_base + '.svg')
