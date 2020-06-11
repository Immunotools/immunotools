import os
import sys
import operator
from datetime import datetime

from Bio import SeqIO

import igraph 
from igraph import *

import dataset
import amino_acid_utils
import utils
import clonal_tree_utils
import graph_utils

######################## SEQUENCE FILTER #############################
class AbundantVJFilter:
    def __init__(self, full_length_lineage):
        self.full_length_lineage = full_length_lineage
        self.most_freq_v = self._FindMostAbundantGene(dataset.AnnotatedGene.V)
        self.most_freq_j = self._FindMostAbundantGene(dataset.AnnotatedGene.J)

    def _FindMostAbundantGene(self, gene_type):
        gene_dict = dict() # gene name -> num sequences
        for seq in self.full_length_lineage.FullLengthSeqIdIter():
            gene_name = self.full_length_lineage.Dataset().GetGeneHitBySeqName(seq.id, gene_type)
            base_name = utils.GetBaseGeneName(gene_name)
            if base_name not in gene_dict:
                gene_dict[base_name] = 0
            gene_dict[base_name] += 1
        most_freq_gene = max(gene_dict.iteritems(), key=operator.itemgetter(1))[0]
#        print "Most frequent " + str(gene_type.name) + ' gene: ' + most_freq_gene + ' (' + str(gene_dict[most_freq_gene]) + ' sequences)'
        return max(gene_dict.iteritems(), key=operator.itemgetter(1))[0]

    def SequenceIsGood(self, seq):
        v_gene = self.full_length_lineage.Dataset().GetGeneHitBySeqName(seq.id, dataset.AnnotatedGene.V)
        j_gene = self.full_length_lineage.Dataset().GetGeneHitBySeqName(seq.id, dataset.AnnotatedGene.J)
        return utils.GetBaseGeneName(v_gene) == self.most_freq_v and utils.GetBaseGeneName(j_gene) == self.most_freq_j

class AbundantLengthFilter:
    def __init__(self, full_length_lineage):
        self.full_length_lineage = full_length_lineage
        self._InitLengthDict()

    def _InitLengthDict(self):
        self.len_dict = dict()
        for seq in self.full_length_lineage.FullLengthSeqIdIter():
            if len(seq) not in self.len_dict:
                self.len_dict[len(seq)] = 0
            self.len_dict[len(seq)] += 1
        self.most_freq_length = max(self.len_dict.iteritems(), key=operator.itemgetter(1))[0]
#        print 'LENGTH: ' + str(float(max(self.len_dict.values())) / sum(self.len_dict.values()) * 100)
#        print "Lengths of sequences:"
#        for it in sorted(self.len_dict.iteritems(), key=operator.itemgetter(1), reverse = True):
#            print str(it[0]) + ' nt: ' + str(it[1]) + ' sequences' 
#        print "Most frequent sequence length: " + str(self.most_freq_length) + ' (' + str(self.len_dict[self.most_freq_length]) + ' sequences)'

    def SequenceIsGood(self, seq):
        return len(seq) == self.most_freq_length

class TrivialFilter:
    def SequenceIsGood(self, seq):
        return True

class NaiveSequenceFilter:
    def __init__(self, full_length_lineage, min_num_shms = 3):
        self.full_length_lineage = full_length_lineage
        self.dataset = full_length_lineage.Dataset()
        self.min_num_shms = min_num_shms

    def SequenceIsGood(self, seq):
        v_shms = self.dataset.GetVSHMsOutsideCDR3(seq.id) #, dataset.AnnotatedGene.V)
        j_shms = self.dataset.GetJSHMsOutsideCDR3(seq.id) #, dataset.AnnotatedGene.J)
        return len(v_shms) + len(j_shms) > self.min_num_shms

class CustomFilter:
    def __init__(self, filters):
        self.filters = filters

    def SequenceIsGood(self, seq):
        for f in self.filters:
            if not f.SequenceIsGood(seq):
                return False
        return True    

#########################  SEQUENCE ITERATOR #########################
class AllSequenceIterator:
    def __init__(self, full_length_lineage):
        self.full_length_lineage = full_length_lineage
        self.seqs = []
        self.seqs.append([seq for seq in self.full_length_lineage.FullLengthSeqIdIter()])
        
    def __iter__(self):
        for seqs in self.seqs:
            yield seqs

class AbundantSequenceIterator:
    def __init__(self, full_length_lineage, min_abs_abundance):
        self.full_length_lineage = full_length_lineage
        self.dataset = self.full_length_lineage.Dataset()
        self.min_abs_abundance = min_abs_abundance

    def __iter__(self):
        abundant_seqs = []
        num_seqs = 0
        for seq in self.full_length_lineage.FullLengthSeqIdIter():
            num_seqs += 1
            if self.dataset.GetSeqMultiplicity(seq.id) < self.min_abs_abundance:
                continue
            abundant_seqs.append(seq)
        print "The number of sequences is reduced from " + str(num_seqs) + ' to ' + str(len(abundant_seqs))
        yield abundant_seqs
            
class AbundantAASequenceIterator:
    def __init__(self, full_length_lineage, min_rel_abundance, min_abs_abundance):
        self.full_length_lineage = full_length_lineage
        self.min_rel_abundance = min_rel_abundance # float number from 0 to 1
        self.min_abs_abundance = min_abs_abundance
        self.aa_dict = amino_acid_utils.AminoAcidDict()
        self.aa_dict.AddLineage(full_length_lineage)

    def __iter__(self):
        for aa in self.aa_dict:
            #if float(self.aa_dict.GetAAMultiplicity(aa)) / len(self.full_length_lineage) < self.min_rel_abundance or self.aa_dict.GetAAMultiplicity(aa) < self.min_abs_abundance:
            if self.aa_dict.GetAAMultiplicity(aa) < self.min_abs_abundance:
                continue
            seq_ids = self.aa_dict.GetIdsByAA(aa)
            seqs = [self.full_length_lineage.GetFullLengthSequenceByName(seq_id) for seq_id in seq_ids]
            yield seqs

class AllAbundantAAsIterator:
    def __init__(self, full_length_lineage, min_rel_abundance, min_abs_abundance):
        self.full_length_lineage = full_length_lineage
        self.min_rel_abundance = min_rel_abundance # float number from 0 to 1
        self.min_abs_abundance = min_abs_abundance
        self.aa_dict = amino_acid_utils.AminoAcidDict()
        self.aa_dict.AddLineage(full_length_lineage)
        abundant_aa_seqs = []
        for aa in self.aa_dict:
            if float(self.aa_dict.GetAAMultiplicity(aa)) / len(full_length_lineage) < self.min_rel_abundance or self.aa_dict.GetAAMultiplicity(aa) < min_abs_abundance:
                continue
#            print 'adding aa with multiplicity ' + str(self.aa_dict.GetAAMultiplicity(aa))
            seq_ids = self.aa_dict.GetIdsByAA(aa)
            seqs = [self.full_length_lineage.GetFullLengthSequenceByName(seq_id) for seq_id in seq_ids]
            abundant_aa_seqs.extend(seqs)
        self.seqs = []
        self.seqs.append(abundant_aa_seqs)

    def __iter__(self):
        for seqs in self.seqs:
            yield seqs

############################# EDGE COMPUTER #############################
class NaiveEdgeComputer:
    def ComputeEdges(self, seqs):
        edge_dict = dict()
        for i in range(len(seqs)):
            for j in range(i + 1, len(seqs)):
                edge_dict[(i, j)] = utils.HammingDistance(seqs[i].seq, seqs[j].seq)
        return edge_dict

class HGToolEdgeComputer:
    def __init__(self, output_dir, hg_running_line, hg_tau = 10):
        self.output_dir = output_dir
        self.hg_running_line = hg_running_line
        if not os.path.exists(self.output_dir):
            os.mkdir(self.output_dir)
        # algorithm params
        self.hg_tau = 10

    def ComputeEdges(self, seqs):
#        print "# sequences: " + str(len(seqs))
        current_time = str(datetime.now()).replace(' ', '_')
        fasta_fname = os.path.join(self.output_dir, 'lineage_' + current_time + '.fasta')
        self._OutputSeqsToFasta(seqs, fasta_fname)
        graph_fname = os.path.join(self.output_dir, 'lineage_' + current_time + '.graph')
        graph = self._ConstructHG(fasta_fname, graph_fname)
#        self.connected_components = graph.GetConnectedComponents()
#        component_sizes = [len(c) for c in self.connected_components]
#        print "Component sizes: " + str(sorted(component_sizes, reverse = True))
        largest_component = self._GetLargestConnectedComponent(graph)
#        print len(largest_component)
        return self._GetComponentEdges(graph, largest_component), float(len(largest_component)) / graph.NumVertices()

    def _OutputSeqsToFasta(self, seqs, fasta_fname):
        fasta_fh = open(fasta_fname, 'w')
        for s in seqs:
            fasta_fh.write('>' + s.id + '\n')
            fasta_fh.write(s.seq + '\n')
        fasta_fh.close()

    def _ConstructHG(self, fasta_fname, graph_fname):
        os.system(self.hg_running_line + ' -i ' + fasta_fname + ' -o ' + graph_fname + ' -k 20 --tau ' + str(self.hg_tau) + ' -T 0 -t 32 > /dev/null') # todo: refactor parameters
        return graph_utils.Graph(graph_fname)

    def _GetLargestConnectedComponent(self, graph):
        connected_components = graph.GetConnectedComponents()
        sorted_components = sorted(connected_components, key = lambda s : len(s), reverse = True)
        print "Connected component: " + str(float(len(sorted_components[0])) / graph.NumVertices())
        return sorted_components[0]

    def _GetComponentEdges(self, graph, component):
        component_edges = dict()
        component_set = set(component)
        for e in graph.EdgeIterator():
            if e[0] not in component_set:
                continue
            component_edges[e] = graph.GetEdgeWeight(e)
        return component_edges

######################## TREE CONSTRUCTOR #############################
class ClonalTreeConstructor:
    def __init__(self, full_length_lineage, seq_iterator, seq_filter, edge_computer, tree_computer, min_tree_size, min_component_frac):
        self.full_length_lineage = full_length_lineage
        self.seq_iterator = seq_iterator
        self.seq_filter = seq_filter
        self.edge_computer = edge_computer
        self.tree_computer = tree_computer
        self.min_tree_size = min_tree_size
        self.min_component_frac = min_component_frac
        self.clonal_trees = []
        self.comp_fractions = []
        self._ProcessLineage()

    def _ProcessLineage(self):
        for seqs in self.seq_iterator:
            filtered_seqs = [seq for seq in seqs if self.seq_filter.SequenceIsGood(seq)]
#            print '  ' + str(len(seqs)) + ' -> ' + str(len(filtered_seqs))
            if len(filtered_seqs) < self.min_tree_size:
                continue
            print "Computing edges of the Hamming graph..."
            edge_weights, component_fraction = self.edge_computer.ComputeEdges(filtered_seqs)
            if component_fraction < self.min_component_frac:
                print "WARN: lineage is likely combined by similar naive sequences and will be skipped"
                continue
            print "Computing MST of the Hamming graph..."
            tree_weights = self.tree_computer.ComputeSpanningTree(filtered_seqs, edge_weights)
            undirected_tree = clonal_tree_utils.UndirectedClonalTree(self.full_length_lineage, filtered_seqs)
            for e in tree_weights:
                undirected_tree.AddEdge(e, tree_weights[e])
            print str(undirected_tree.NumVertices()) + ' vertices in undirected clonal tree'
            if undirected_tree.NumVertices() == 0:
                print "WARN: tree is empty"
                continue
            root_finder = clonal_tree_utils.SimpleRootComputer(undirected_tree)
            directed_tree = clonal_tree_utils.DirectedClonalTree(undirected_tree, root_finder.GetRoot())
            self.clonal_trees.append(directed_tree)
            self.comp_fractions.append(component_fraction)

    def GetClonalTrees(self):
        return self.clonal_trees

    def GetComponentFractions(self):
        return self.comp_fractions
