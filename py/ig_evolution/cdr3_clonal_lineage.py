import os
import sys
import shutil
from enum import Enum
import pandas as pd
from Bio import SeqIO
import igraph 
from igraph import *
import dataset

############################## CLONAL LINEAGE CONSTRUCTOR ##############
class ClonalLineageConfig:
    hg_constructor = 'hg_utils/./construct_graph'
    hg_processor = 'python hg_utils/process_graph.py'
    cdr3_decomposer = 'python hg_utils/create_clonal_lineages.py'

class CDR3Lineage:
    def __init__(self, dataset, lineage_fasta_fname, lineage_graph_fname, lineage_index):
        self.dataset = dataset
        self._InitLineage(lineage_fasta_fname, lineage_graph_fname)
        self.graph = Graph()
        self.graph.add_vertices(len(self))
        self._InitEdges()
        self._InitFullLengthSeqs()
        self._InitId(lineage_index)
#        print "Sanity check: " + str(lineage_index)
#        for cdr3 in self.CDR3Iter():
#            print cdr3

    def _InitId(self, lineage_index):
        self.lineage_id = 'lineage' + str(lineage_index) #+ '_numCDR3s' + str(len(self)) + '_numAbs' + str(self.NumFullLengthSequences())

    def _LineIsVertex(self, l):
        return l.find('->') == -1

    def _AddVertex(self, l, vertices):
        splits = l.strip().split()
        vertices.append(int(splits[0]))

    def _AddEdge(self, l, edges):
        splits = l.strip().split()
        v1 = int(splits[0])
        v2 = int(splits[2])
        weight = int(splits[5][ : len(splits[5]) - len('];')])
        edges[(v1, v2)] = weight

    def _InitLineage(self, lineage_fasta_fname, lineage_graph_fname):
        self.dataset_indices = []
        for r in SeqIO.parse(lineage_fasta_fname, 'fasta'):
            cur_index = int(r.id.split('|')[0].split(':')[1])
            self.dataset_indices.append(cur_index)
        local_ind_vertices, local_ind_edges = self._ReadLineageGraph(lineage_graph_fname)
        if len(self.dataset_indices) != len(local_ind_vertices):
            print "ERROR: the number of vertices in " + lineage_graph_fname + ' does not match with the number of sequences in ' + lineage_fasta_fname 
        local_dataset_ind_dict = dict() # local index -> dataset index
        for i in range(len(self.dataset_indices)):
            local_dataset_ind_dict[local_ind_vertices[i]] = self.dataset_indices[i]
        self.edges = dict()
        for e in local_ind_edges:
            dataset_ind1 = local_dataset_ind_dict[e[0]]
            dataset_ind2 = local_dataset_ind_dict[e[1]]
            self.edges[(dataset_ind1, dataset_ind2)] = local_ind_edges[e]
   
    def _ReadLineageGraph(self, lineage_graph_fname):
        vertices = []
        edges = dict()
        lines = open(lineage_graph_fname, 'r').readlines()
        for i in range(1, len(lines) - 1):
            l = lines[i].strip()
            if self._LineIsVertex(l):
                self._AddVertex(l, vertices)
            else:
                self._AddEdge(l, edges)
        return vertices, edges

    def _InitEdges(self):
        index_dict = dict() # key: index from 0 to NUM_CDR3s_IN_LINEAGE to index in the cdr3 list
        index = 0
        for ind in self.dataset_indices:
            index_dict[ind] = index
            index += 1
        graph_edges = []
        distances = []
        for e in self.edges:
            v1 = index_dict[e[0]] #lineage_indices.index(e[0])
            v2 = index_dict[e[1]] #lineage_indices.index(e[1])
            distances.append(self.edges[e])
            graph_edges.append((v1, v2))
#        print "Adding edges..."
        self.graph.add_edges(graph_edges)
#        print "Finding spanning tree..."
        spanning_tree = self.graph.spanning_tree(weights = distances, return_tree = True)
        index_list = index_dict.keys() 
        self.tree_edges = []
        for e in spanning_tree.get_edgelist():
            self.tree_edges.append((self.dataset_indices[e[0]], self.dataset_indices[e[1]]))

    def _InitFullLengthSeqs(self):
        self.full_length_seqs = []
        for cdr3_ind in self.dataset_indices:
            self.full_length_seqs.extend(self.dataset.GetSeqIdsByCDR3Index(cdr3_ind))

    def __len__(self):
        return len(self.dataset_indices)

    def __hash__(self):
        return hash(self.lineage_id)

    def __eq__(self, other):
        return self.lineage_id == other.lineage_id

    def id(self):
        return self.lineage_id

    def __str__(self):
        return self.id()

    def __repr__(self):
        return self.id()

    def Trivial(self):
        return len(self) == 1

    def NumFullLengthSequences(self):
        return len(self.full_length_seqs)

    def FullLengthIter(self):
        for fl_id in self.full_length_seqs:
            yield self.dataset.GetSeqByName(fl_id)

    def GetFullLengthSequenceByIndex(self, seq_index):
        return self.dataset.GetSeqByName(self.full_length_seqs[seq_index])

    def CDR3Iter(self):
        for cdr3_ind in self.dataset_indices:
            yield self.dataset.GetCDR3ByIndex(cdr3_ind)

    def Dataset(self):
        return self.dataset

class CDR3LineageStorage:
    def __init__(self, dataset, graph_dir, lineage_dir):
        self.lineages = []
        self.dataset = dataset
        subgraph_dirs = os.listdir(graph_dir)
        lineage_index = 0
        for d in subgraph_dirs:
            subgraph_dir = os.path.join(graph_dir, d)
            graph_files = os.listdir(subgraph_dir)
            for f in graph_files:
                if f == 'graph_stats.txt':
                    continue
                lineage_graph = os.path.join(subgraph_dir, f)
                lineage_fasta = os.path.join(os.path.join(lineage_dir, d), f.split('.')[0] + '.fasta')
                self.lineages.append(CDR3Lineage(dataset, lineage_fasta, lineage_graph, lineage_index))
                lineage_index += 1
        print str(len(self.lineages)) + " lineages were extracted from " + graph_dir + ' & ' + lineage_dir

    def __iter__(self):
        for l in self.lineages:
            yield l

    def __len__(self):
        return len(self.lineages)

class CDR3LineageConstructor:
    def __init__(self, dataset, output_dirs):
        self.dataset = dataset
        self.output_dir = output_dirs['main_dir']
        # output files and dirs
        self.cdr3_fasta_dir = output_dirs['cdr3_fasta_dir'] #os.path.join(self.output_dir, 'cdr3_fasta')
        self.cdr3_graph_dir = output_dirs['cdr3_graph_dir'] #os.path.join(self.output_dir, 'cdr3_graphs')
        self.connected_comp_dir = output_dirs['cdr3_comp_dir'] #os.path.join(self.output_dir, "cdr3_connected_components")
        self.lineage_dir = output_dirs['cdr3_lin_dir'] #os.path.join(self.output_dir, 'cdr3_lineages')

    def Construct(self):
        self._CheckIntermediateFiles()
        self._OutputCDR3sByLengthLevels()
        self._ConstructCDR3GraphsByLengthLevel()
        self._ComposeCDR3GraphsIntoConnectedComponents()
        self._DecomposeCDR3s()
        self.cdr3_lineages = CDR3LineageStorage(self.dataset, self.connected_comp_dir, self.lineage_dir)
        return self.cdr3_lineages

    def _CreateNonExistingDir(self, dir_name):
        if not os.path.exists(dir_name):
            os.mkdir(dir_name)

    def _CheckIntermediateFiles(self):
        # TODO:
        self._CreateNonExistingDir(self.cdr3_fasta_dir)
        self._CreateNonExistingDir(self.cdr3_graph_dir)
        self._CreateNonExistingDir(self.connected_comp_dir)
        self._CreateNonExistingDir(self.lineage_dir)

    def _OutputCDR3sByLengthLevels(self):
        print "Writing CDR3s to FASTA files according to their lengths"
        length_level = 30
        length_tau = 3
        cdr_length_dict = dict()
        self.cdr3_index_map = dict() # index of CDR3 in FASTA -> index of CDR3 in dataset
        cdr3_index = 0
        for i in range(self.dataset.NumDistinctCDR3s()):
            if self.dataset.GetCDR3ByIndex(i) == '-':
                continue
            cdr3_length = len(self.dataset.GetCDR3ByIndex(i))            
            self.cdr3_index_map[cdr3_index] = i
            cdr3_index += 1
            if cdr3_length not in cdr_length_dict:
                cdr_length_dict[cdr3_length] = []
            cdr_length_dict[cdr3_length].append(i)
        start_level = min(cdr_length_dict.keys()) / length_level * length_level
        end_level = (max(cdr_length_dict.keys()) / length_level + 1) * length_level
        levels = range(start_level, end_level, length_level)
        fhandler_dict = dict()
        self.file_distance_dict = dict()
        for l in levels:
            fname = os.path.join(self.cdr3_fasta_dir, "cdr3s_" + str(l) + '_' + str(l + length_level) + '.fasta')
            self.file_distance_dict[fname] = (l / length_level + 1) * length_tau
            fhandler_dict[l] = open(fname, 'w')
        for cdr3_len in sorted(cdr_length_dict):
            cdr3_len_level = cdr3_len / length_level * length_level
            fh = fhandler_dict[cdr3_len_level]
            for cdr3_index in cdr_length_dict[cdr3_len]:
                fh.write('>INDEX:' + str(cdr3_index) + '|MULT:' + str(self.dataset.GetCDR3Multiplicity(cdr3_index)) + '\n')
                fh.write(self.dataset.GetCDR3ByIndex(cdr3_index) + '\n')
        for l in fhandler_dict:
            fhandler_dict[l].close()

    def _ConstructCDR3GraphsByLengthLevel(self):
        print "Constructing Hamming graph on CDR3s..."
        fasta_files = os.listdir(self.cdr3_fasta_dir)
        for f in fasta_files:
            full_path = os.path.join(self.cdr3_fasta_dir, f)
            print "  Constructing Hamming graph on " + full_path + ' with tau = ' + str(self.file_distance_dict[full_path])
            base_name = f.split('.')[0]
            os.system(ClonalLineageConfig.hg_constructor + ' ' + full_path + ' ' + os.path.join(self.cdr3_graph_dir, base_name + '.graph') + ' ' + str(self.file_distance_dict[full_path]) + ' > /dev/null')

    def _ComposeCDR3GraphsIntoConnectedComponents(self):
        print "Decomposing graph into connected components..."
        graph_files = os.listdir(self.cdr3_graph_dir)
        for f in graph_files:
            print "  Decomposing " + os.path.join(self.cdr3_graph_dir, f) + '...'
            basename = f.split('.')[0]
            component_subdir = os.path.join(self.connected_comp_dir, basename)
            os.mkdir(component_subdir)
            os.system(ClonalLineageConfig.hg_processor + ' ' + os.path.join(self.cdr3_graph_dir, f) + ' ' + component_subdir + ' > /dev/null')

    def _DecomposeCDR3s(self):
        print "Decomposing CDR3s into CDR3 clonal lineages..."
        self.lineage_dir = os.path.join(self.output_dir, 'cdr3_lineages')
        component_dirs = os.listdir(self.connected_comp_dir)
        for d in component_dirs:
            lineage_subdir = os.path.join(self.lineage_dir, d)
            cdr3_fasta = os.path.join(self.cdr3_fasta_dir, d + '.fasta')
            component_dir = os.path.join(self.connected_comp_dir, d)
            os.system(ClonalLineageConfig.cdr3_decomposer + ' ' + cdr3_fasta + ' ' + component_dir + ' ' + lineage_subdir + ' > /dev/null')

    def _RemoveExtraFiles(self):
        os.system('rm ' + self.cdr3_fasta + ' ' + self.cdr3_graph)
        os.system('rm -rf ' + self.graph_dir + ' ' + self.lineage_dir)

