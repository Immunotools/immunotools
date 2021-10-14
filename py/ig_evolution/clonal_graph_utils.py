import os
import sys
import operator

import matplotlib as mplt
mplt.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.patches as patches
from matplotlib.patches import Rectangle
import seaborn as sns
sns.set()

import utils
import dataset
import clonal_tree_constructor
import amino_acid_utils
import vj_annotator
import clonal_tree_writer
import mst_algorithms
import clonal_tree_simplification
import clonal_tree_stats
import clonal_graph_algorithms
import shm_visualizer

class ClonalGraph:
    def __init__(self, clonal_tree):
        self.clonal_tree = clonal_tree
        self.full_length_lineage = self.clonal_tree.FullLengthLineage()
        self.vj_annotator = vj_annotator.VJGeneAnnotator(self.full_length_lineage)
        self._InitEdges()
        self._InitRoot()
        self._InitAdjList()
        self._InitSHMs()

    def _InitEdges(self):
        self.aa_dict = amino_acid_utils.AminoAcidDict()
        self.aa_dict.AddClonalTree(self.clonal_tree)
        self.used_aa = []
        self.aa_edges = dict()
        self.vertex_indices = set()
        self.num_edges = 0
        for v in self.clonal_tree.VertexIter():
            v_id = self.clonal_tree.GetSequenceByVertex(v).id
            v_aa = self.aa_dict.GetAAById(v_id)
            v_index = self.aa_dict.GetIndexByAA(v_aa)
            self.vertex_indices.add(v_index)
        for e in self.clonal_tree.EdgeIter():
            self.num_edges += 1
            src_id = self.clonal_tree.GetSequenceByVertex(e[0]).id
            dst_id = self.clonal_tree.GetSequenceByVertex(e[1]).id
            src_aa = self.aa_dict.GetAAById(src_id)
            dst_aa = self.aa_dict.GetAAById(dst_id)
            src_index = self.aa_dict.GetIndexByAA(src_aa)
            dst_index = self.aa_dict.GetIndexByAA(dst_aa)
#            self.vertex_indices.add(src_index)
#            self.vertex_indices.add(dst_index)
            if src_aa == dst_aa:
                continue
            aa_edge = (src_index, dst_index)
            if aa_edge not in self.aa_edges:
                self.aa_edges[aa_edge] = 0
            self.aa_edges[aa_edge] += 1
        print "# amino acid vertices: " + str(len(self.vertex_indices))

    def _InitRoot(self):
        root_aa = self.aa_dict.GetAAById(self.clonal_tree.RootSeq().id)
        self.root = self.aa_dict.GetIndexByAA(root_aa)

    def _InitAdjList(self):
        self.incoming_lists = dict()
        self.outgoing_lists = dict()
        for v in self.vertex_indices:
            self.incoming_lists[v] = []
            self.outgoing_lists[v] = []
        for e in self.aa_edges:
            self.incoming_lists[e[1]].append(e[0])
            self.outgoing_lists[e[0]].append(e[1])

    def _InitSHMs(self):
        self.shms = TreeSHMs(self.clonal_tree.FullLengthLineage(), self.aa_dict)
        for e in self.aa_edges:
            src_aa = self.aa_dict.GetAAByIndex(e[0])
            dst_aa = self.aa_dict.GetAAByIndex(e[1])
            for i in range(len(src_aa)):
                if src_aa[i] == dst_aa[i]:
                    continue
                shm = dataset.SHM(i, -1, src_aa[i], dst_aa[i])
                self.shms.AddSHM(e, shm)

    def FullLengthLineage(self):
        return self.full_length_lineage

    def GetGraphSHMs(self):
        return self.shms

    def VertexIter(self):
        for v in self.vertex_indices:
            yield v

    def GetRootIndex(self):
        return self.root

    def GetVertexMultiplicity(self, v):
        return sum(self.GetVertexMultiplicities(v))

    def GetVertexMultiplicities(self, v):
        aa_id = self.aa_dict.GetAAByIndex(v)
        seq_ids = self.aa_dict.GetIdsByAA(aa_id)
        return [self.full_length_lineage.Dataset().GetSeqMultiplicity(seq_id) for seq_id in seq_ids]

    def GetVertexDiversity(self, v):
        return len(self.aa_dict.GetIdsByAA(self.aa_dict.GetAAByIndex(v)))

    def GetNumVSHMsByVertex(self, v):
        aa_id = self.aa_dict.GetAAByIndex(v)
        seq_ids = self.aa_dict.GetIdsByAA(aa_id)
        return [len(self.full_length_lineage.Dataset().GetVSHMsOutsideCDR3(seq_id)) for seq_id in seq_ids]

    def GetVertexLabels(self, v):
        aa_id = self.aa_dict.GetAAByIndex(v)
        return self._GetAALabels(aa_id)

    def _GetAALabels(self, aa_id):
        aa_labels = []
        seq_ids = self.aa_dict.GetIdsByAA(aa_id)
        for seq_id in seq_ids:
            labels = self.full_length_lineage.Dataset().GetSeqLabels(seq_id)
            aa_labels.extend(labels)
        aa_labels = set(aa_labels)
        return aa_labels

    def EdgeIter(self):
        for e in self.aa_edges:
            yield e

    def GetIncomingVertices(self, v):
        return self.incoming_lists[v]

    def GetOutgoingVertices(self, v):
        return self.outgoing_lists[v]

    def IsLeaf(self, v):
        return len(self.outgoing_lists[v]) == 0

    def NumVertices(self):
        return len(self.vertex_indices)

    def NumEdges(self):
        return self.num_edges

    def GetSHMs(self):
        return self.shms

    def _GetDominantGenes(self):
        v_dict = dict()
        j_dict = dict()
        for v in self.VertexIter():
            aa_id = self.aa_dict.GetAAByIndex(v)
            seq_ids = self.aa_dict.GetIdsByAA(aa_id)
            for s_id in seq_ids:
                v_gene = self.full_length_lineage.Dataset().GetGeneHitBySeqName(s_id, dataset.AnnotatedGene.V)
                j_gene = self.full_length_lineage.Dataset().GetGeneHitBySeqName(s_id, dataset.AnnotatedGene.J)
                if v_gene not in v_dict:
                    v_dict[v_gene] = 0
                v_dict[v_gene] += 1
                if j_gene not in j_dict:
                    j_dict[j_gene] = 0
                j_dict[j_gene] += 1
        return max(v_dict.iteritems(), key=operator.itemgetter(1))[0], max(j_dict.iteritems(), key=operator.itemgetter(1))[0]

    def OutputGraphSHMsToTxt(self, output_fname):
        fh = open(output_fname, 'w')
        v_gene, j_gene = self._GetDominantGenes()
        fh.write('Position\tSrc_AA\tDst_AA\tEdges\tMultiplicity\tRegion\tHas_reverse\tV_gene\tJ_gene\n')
        for shm in self.shms.SHMIter():
            edge_str = ','.join([str(e[0]) + '-' + str(e[1]) for e in self.shms.GetEdgesBySHM(shm)])
            fh.write(str(shm.pos) + '\t' + shm.src_n + '\t' + shm.dst_n + '\t' + edge_str + '\t' + str(self.shms.GetSHMMultiplicity(shm)) + '\t' + self.shms.GetRegionForSHM(shm).name + '\t' + str(self.shms.SHMHasReverse(shm)) + '\t' + v_gene + '\t' + j_gene + '\n')
        fh.close()

    def OutputPutativeAASeqs(self, output_fname):
        v_gene, j_gene = self._GetDominantGenes()
        fh = open(output_fname, 'w')
        fh.write('Index\tAA_seq\tNuclSeqs\tOriginal_mults\tOriginal_headers\tOriginal_labels\tCDR1\tCDR2\tCDR3\tV_gene\tJ_gene\tV_SHMs\tJ_SHMs\n')
        for aa in self.aa_dict:
            aa_ids = self.aa_dict.GetIdsByAA(aa)
            vertex_labels = self._GetAALabels(aa)
            nucl_mults = sorted([self.full_length_lineage.Dataset().GetSeqMultiplicity(seq_id) for seq_id in aa_ids], reverse = True)
            v_shms = sorted([len(self.full_length_lineage.Dataset().GetVSHMsOutsideCDR3(seq_id)) for seq_id in aa_ids], reverse = True)
            j_shms = sorted([len(self.full_length_lineage.Dataset().GetJSHMsOutsideCDR3(seq_id)) for seq_id in aa_ids], reverse = True)
            nucls_seqs = [self.full_length_lineage.Dataset().GetSeqByName(seq_id).seq for seq_id in aa_ids]
            headers = ';'.join(aa_ids)
            cdr1 = aa[self.shms.cdr1_bounds[0] : self.shms.cdr1_bounds[1] + 1]
            cdr2 = aa[self.shms.cdr2_bounds[0] : self.shms.cdr2_bounds[1] + 1]
            cdr3 = aa[self.shms.cdr3_bounds[0] : self.shms.cdr3_bounds[1] + 1]
            fh.write(str(self.aa_dict.GetIndexByAA(aa)) + '\t' + aa + '\t' + ','.join(nucls_seqs) + '\t' + ','.join([str(m) for m in nucl_mults]) + '\t' + headers + '\t' + ','.join(vertex_labels) + '\t' + cdr1 + '\t' + cdr2 + '\t' + cdr3 + '\t' + v_gene + '\t' + j_gene + '\t' + ','.join([str(n) for n in v_shms]) + '\t' + ','.join([str(n) for n in j_shms]) + '\n')
        fh.close()

    def _PositionIsInCDRs(self, aa_pos):
        if aa_pos >= self.shms.cdr1_bounds[0] and aa_pos <= self.shms.cdr1_bounds[1]:
            return True
        if aa_pos >= self.shms.cdr2_bounds[0] and aa_pos <= self.shms.cdr2_bounds[1]:
            return True
        if aa_pos >= self.shms.cdr3_bounds[0] and aa_pos <= self.shms.cdr3_bounds[1]:
            return True
        return False

    def OutputGraphSHMsAsMatrix(self, vertex_orderer, output_base):
        if self.NumVertices() < 2:
            return
        matrix = []
        #vertex_order = vertex_orderer.GetOrder()
        aa_len = len(self.aa_dict.GetAAByIndex(self.GetRootIndex()))
        for e in self.EdgeIter():
            matrix.append([0] * aa_len)
        for i in range(len(matrix)):
            for j in range(len(matrix[i])):
                if self._PositionIsInCDRs(j):
                    matrix[i][j] = 1
        edge_ind = 0
        for e in self.EdgeIter():
            src_aa = self.aa_dict.GetAAByIndex(e[0])
            dst_aa = self.aa_dict.GetAAByIndex(e[1])
            for i in range(len(src_aa)):
                if src_aa[i] != dst_aa[i]:
                    matrix[edge_ind][i] = 2
            edge_ind += 1
#        for i in range(1, len(vertex_order)):
#            cur_aa = self.aa_dict.GetAAByIndex(vertex_order[i])
#            parent_aa = self.aa_dict.GetAAByIndex(self.GetIncomingVertices(vertex_order[i])[0])
#            for j in range(len(cur_aa)):
#                if cur_aa[j] != parent_aa[j]:
#                    matrix[i][j] = 2
#        vertex_levels = clonal_graph_algorithms.GetLevelsByVertexOrder(self, vertex_order)
#        level_colors = []
#        for l in vertex_levels:
#            level_colors.append(utils.GetColorByNormalizedValue('prism', float(l) / max(vertex_levels)))
        sns.heatmap(matrix, cmap = 'coolwarm', xticklabels = [], yticklabels = [], cbar = False)
#        sns.clustermap(matrix, cmap = 'coolwarm', yticklabels = [str(v) for v in vertex_order], row_colors = level_colors, row_cluster = False, col_cluster = False, xticklabels = [])
        plt.savefig(output_base + ".svg")
        utils.OutputPlotToPdf(output_base + '.pdf')

    def OutputSHMsPerPosition(self, output_fname):
        aa_len = len(self.aa_dict.GetAAByIndex(self.GetRootIndex()))
        num_shms = [0] * aa_len
        for shm in self.shms.SHMIter():
            num_shms[shm.pos] += self.shms.GetSHMMultiplicity(shm) * self.shms.GetSHMMultiplicity(shm)
        plt.figure()
        ax = plt.gca()
        ax.add_patch(Rectangle((self.shms.cdr1_bounds[0], 0), self.shms.cdr1_bounds[1] - self.shms.cdr1_bounds[0], max(num_shms), facecolor='#FFB4B6'))
        ax.add_patch(Rectangle((self.shms.cdr2_bounds[0], 0), self.shms.cdr2_bounds[1] - self.shms.cdr2_bounds[0], max(num_shms), facecolor='#FFB4B6'))
        ax.add_patch(Rectangle((self.shms.cdr3_bounds[0], 0), self.shms.cdr3_bounds[1] - self.shms.cdr3_bounds[0], max(num_shms), facecolor='#FFB4B6'))
        plt.bar(range(aa_len), num_shms)
        plt.xlabel('AA position')
        plt.ylabel('Sum of squared SHM multiplicities')
        plt.savefig(output_fname + '.svg')
        utils.OutputPlotToPdf(output_fname + '.pdf')

class TreeSHMs:
    def __init__(self, full_length_lineage, aa_dict):
        self.full_length_lineage = full_length_lineage
        self.aa_seqs = [aa for aa in aa_dict]
        self.aa_dict = aa_dict
        self.shm_edge_map = dict()
        self.edge_shm_map = dict()
        self._InitializeCDRBounds()

    def _InitializeCDRBounds(self):
        first_seq_id = self.aa_dict.GetIdsByAA(self.aa_seqs[0])[0]
        self.cdr1_bounds = self.full_length_lineage.Dataset().GetCDR1BoundsBySeqName(first_seq_id)
        self.cdr1_bounds = (self.cdr1_bounds[0] / 3, self.cdr1_bounds[1] / 3)
        self.cdr2_bounds = self.full_length_lineage.Dataset().GetCDR2BoundsBySeqName(first_seq_id)
        self.cdr2_bounds = (self.cdr2_bounds[0] / 3, self.cdr2_bounds[1] / 3)
        self.cdr3_bounds = self.full_length_lineage.Dataset().GetCDR3BoundsBySeqName(first_seq_id)
        self.cdr3_bounds = (self.cdr3_bounds[0] / 3, self.cdr3_bounds[1] / 3)

    def AddSHM(self, edge, shm):
        if edge not in self.edge_shm_map:
            self.edge_shm_map[edge] = []
        self.edge_shm_map[edge].append(shm)
        if shm not in self.shm_edge_map:
            self.shm_edge_map[shm] = []
        self.shm_edge_map[shm].append(edge)

    def GetRegionForSHM(self, shm):
        if shm.pos < self.cdr1_bounds[0]:
            return utils.StructuralRegion.FR1
        if shm.pos < self.cdr1_bounds[1]:
            return utils.StructuralRegion.CDR1
        if shm.pos < self.cdr2_bounds[0]:
            return utils.StructuralRegion.FR2
        if shm.pos < self.cdr2_bounds[1]:
            return utils.StructuralRegion.CDR2
        if shm.pos < self.cdr3_bounds[0]:
            return utils.StructuralRegion.FR3
        if shm.pos < self.cdr3_bounds[1]:
            return utils.StructuralRegion.CDR3
        return utils.StructuralRegion.FR4

    def SHMIsRepetitive(self, shm):
        return self.GetSHMMultiplicity(shm) > 1

    def EdgeContainsRepetitiveSHMs(self, edge):
        for shm in self.edge_shm_map[edge]:
            if self.SHMIsRepetitive(shm):
                return True
        return False

    def SHMIter(self):
        for shm in sorted(self.shm_edge_map, key = lambda x : x.pos):
            yield shm

    def GetEdgesBySHM(self, shm):
        return self.shm_edge_map[shm] 

    def GetSHMMultiplicity(self, shm):
        return len(set([e[0] for e in self.shm_edge_map[shm]]))

    def GetNumSHMsOnEdge(self, edge):
        return len(self.edge_shm_map[edge])

    def SHMIsNeutral(self, shm):
        reverse_shm = dataset.SHM(shm.pos, -1, shm.dst_n, shm.src_n)
        return reverse_shm in self.shm_edge_map

    def SHMHasReverse(self, shm):
        reverse_shm = dataset.SHM(shm.pos, -1, shm.dst_n, shm.src_n)
        return self.ContainsSHM(reverse_shm)

    def GetReverseSHM(self, shm):
        return dataset.SHM(shm.pos, -1, shm.dst_n, shm.src_n)

    def ContainsSHM(self, shm):
        return shm in self.shm_edge_map

    def Print(self):
        for e in self.edge_shm_map:
            print e, [(str(shm) + ' : ' + str(len(self.shm_edge_map[shm]))) for shm in self.edge_shm_map[e]]

def OutputClonalTree(clonal_tree, full_length_lineage, output_base, vertex_writer):
    if clonal_tree.NumVertices() > 1000:
        return
    #vertex_writer = clonal_tree_writer.LevelMultiplicityVertexWriter(clonal_tree) #clonal_tree_writer.UniqueAAColorWriter(clonal_tree) #.MultiplicityVertexWriter(clonal_tree)
    edge_writer = clonal_tree_writer.SimpleEdgeWriter(clonal_tree) #TypeEdgeWriter(clonal_tree)
    tree_writer = clonal_tree_writer.ClonalTreeWriter(clonal_tree, vertex_writer, edge_writer)
    tree_writer.Output(output_base)

def OutputLineageAminoAcids(full_length_lineage, aa_dict, output_fname):
    fh = open(output_fname, 'w')
    fh.write('Index\tSeq\tDiversity\tNucl_mults\n')
    for aa in aa_dict:
        aa_ids = aa_dict.GetIdsByAA(aa)
        nucl_mults = sorted([full_length_lineage.Dataset().GetSeqMultiplicity(seq_id) for seq_id in aa_ids], reverse = True)
        fh.write(str(aa_dict.GetIndexByAA(aa)) + '\t' + aa + '\t' + str(aa_dict.GetAAMultiplicity(aa)) + '\t' + ','.join([str(m) for m in nucl_mults]) + '\n')
    fh.close()

############################################################################################
def PrintAbundantSequences(full_length_lineage):
    dataset = full_length_lineage.Dataset()
    for seq in full_length_lineage.FullLengthSeqIdIter():
        if dataset.GetSeqMultiplicity(seq.id) > 100:
            print '>' + seq.id
            print seq.seq

def OutputCompressedGraph(clonal_graph, compressed_paths, output_base):
    vertex_path_dict = dict()
    for v in clonal_graph.VertexIter():
        vertex_path_dict[v] = -v - 1
    for ind, path in enumerate(compressed_paths):
        for v in path:
            vertex_path_dict[v] = ind
    fh = open(output_base + '.dot', 'w')
    fh.write('digraph{\n')
    for v in vertex_path_dict:
        if vertex_path_dict[v] >= 0:
            path_ind = vertex_path_dict[v]
            fh.write('p' + str(path_ind) + ' [shape = box, style = rounded, label = \"path (' + str(len(compressed_paths[path_ind])) + ')\"]\n')
        else:
            fh.write('v' + str(v) + ' [label = ' + str(v) + ']\n')
    for e in clonal_graph.EdgeIter():
        if vertex_path_dict[e[0]] == vertex_path_dict[e[1]]:
            continue
        v1_ind = e[0]
        if vertex_path_dict[v1_ind] >= 0:
            v1_ind = 'p' + str(vertex_path_dict[v1_ind])
        else:
            v1_ind = 'v' + str(v1_ind)
        v2_ind = e[1]
        if vertex_path_dict[v2_ind] >= 0:
            v2_ind = 'p' + str(vertex_path_dict[v2_ind])
        else:
            v2_ind = 'v' + str(v2_ind)
        fh.write(v1_ind + ' -> ' + v2_ind + '\n')
    fh.write('}')
    fh.close()
    os.system('dot -Tpdf ' + output_base + '.dot' + ' -o ' + output_base + '.pdf')
    os.system('dot -Tsvg ' + output_base + '.dot' + ' -o ' + output_base + '.svg')

############################################################################################
def DefineClonalGraphName(clonal_graph, lineage_id, comp_frac):
    splits = lineage_id.split('_')
    return splits[0] + '_vertices' + str(clonal_graph.NumVertices()) + '_edges' + str(clonal_graph.NumEdges()) + '_cfrac' + '{:0.2f}'.format(comp_frac)

def OutputAbundantAAGraphs(full_length_lineages, output_dirs, config):
    lineage_index = 1
    for l in sorted(full_length_lineages, key = lambda s : len(s), reverse = True):
        if len(l) < config.min_lineage_size:
            continue
        cdr3_length = l.CDR3Length()
        if cdr3_length < config.min_cdr3_len or cdr3_length > config.max_cdr3_len:
            continue
        # clonal tree construction step
        print "== Processing lineage " + l.id() + ', CDR3 length: ' + str(cdr3_length) + '...'
        custom_filter = clonal_tree_constructor.CustomFilter([clonal_tree_constructor.AbundantLengthFilter(l), clonal_tree_constructor.AbundantVJFilter(l)]) #, clonal_tree_constructor.NaiveSequenceFilter(l)]) 
        seq_iterator = clonal_tree_constructor.AllSequenceIterator(l)
        edge_computer = clonal_tree_constructor.HGToolEdgeComputer(output_dirs['fl_lineages'], 'build/release/bin/./ig_swgraph_construct', config.hg_tau) # TODO: refactor
        tree_computer = mst_algorithms.VertexMultMSTFinder(l) #mst_algorithms.IGraphMSTFinder() # mst_algorithms.VertexMultMSTFinder(l)
        tree_constructor = clonal_tree_constructor.ClonalTreeConstructor(l, seq_iterator, custom_filter, edge_computer, tree_computer, config.min_graph_size, config.min_component_fraction)
        clonal_trees = tree_constructor.GetClonalTrees()
        comp_fractions = tree_constructor.GetComponentFractions()
        if len(clonal_trees) == 0:
            continue
        clonal_tree = clonal_trees[0]
        comp_fraction = comp_fractions[0]
#        PrintAbundantSequences(l)
        # writing clonal tree
        stats_writer = clonal_tree_stats.ClonalTreeStatsWriter(clonal_tree)
        stats_writer.Output(os.path.join(output_dirs['clonal_graphs'], l.id() + '_before_simpl.txt'))
        #OutputClonalTree(clonal_tree, l, os.path.join(output_dirs['clonal_graphs'], l.id() + '_before_simpl'), clonal_tree_writer.LevelMultiplicityVertexWriter(clonal_tree))
        # simplification step
#        print "# vertices before simplification: " + str(clonal_tree.NumVertices())
        leaf_filters = [clonal_tree_simplification.LowFixedAbundanceLeafRemover(clonal_tree, config.min_abs_abundance), clonal_tree_simplification.RelativeAbundanceLeafRemover(clonal_tree, config.max_rel_abundance)]
        leaf_remover = clonal_tree_simplification.IterativeTipRemover(clonal_tree, leaf_filters)
        cleaned_tree = leaf_remover.CleanTips()
        print "# vertices after simplification: " + str(cleaned_tree.NumVertices())
        if cleaned_tree.NumVertices() < config.min_graph_size:
            continue
        # writing clonal tree
        stats_writer = clonal_tree_stats.ClonalTreeStatsWriter(cleaned_tree)
        stats_writer.Output(os.path.join(output_dirs['clonal_graphs'], l.id() + '_after_simpl.txt'))
        #OutputClonalTree(cleaned_tree, l, os.path.join(output_dirs['clonal_graphs'], l.id() + '_after_simpl'), clonal_tree_writer.LevelMultiplicityVertexWriter(cleaned_tree))
        #OutputClonalTree(cleaned_tree, l, os.path.join(output_dirs['clonal_graphs'], l.id() + '_nucl_tree'), clonal_tree_writer.UniqueAAColorWriter(cleaned_tree))
        # aa graph
        clonal_graph = ClonalGraph(clonal_tree)
        if clonal_graph.NumVertices() < config.min_graph_size:
            continue
        graph_name = DefineClonalGraphName(clonal_graph, l.id(), comp_fraction)
        clonal_graph.OutputGraphSHMsToTxt(os.path.join(output_dirs['clonal_graphs'], graph_name + '_shms.txt'))
        clonal_graph.OutputPutativeAASeqs(os.path.join(output_dirs['clonal_graphs'], graph_name + '_seqs.txt'))
        clonal_graph.OutputGraphSHMsAsMatrix(clonal_graph_algorithms.DFSVertexOrderFinder(clonal_graph), os.path.join(output_dirs['shm_matrix'], graph_name))
        clonal_graph.OutputSHMsPerPosition(os.path.join(output_dirs['shm_plot'], graph_name))
        # graph writing
        # multiplicity writing
        edge_writer = clonal_tree_writer.SimpleEdgeWriter(clonal_graph)
        vertex_mult_writer = clonal_tree_writer.LevelMultiplicityVertexWriter(clonal_graph)
        writer = clonal_tree_writer.ClonalGraphVisualizer(clonal_graph, vertex_mult_writer, edge_writer)
        writer.Output(os.path.join(output_dirs['coloring_mult'], graph_name))
        # diversity writing
        vertex_shm_writer = clonal_tree_writer.SHMDepthVertexWriter(clonal_graph)
        writer = clonal_tree_writer.ClonalGraphVisualizer(clonal_graph, vertex_shm_writer, edge_writer)
        writer.Output(os.path.join(output_dirs['coloring_shm'], graph_name))
        # label writing
        vertex_label_writer = clonal_tree_writer.LabelVertexWriter(clonal_graph)
        writer = clonal_tree_writer.ClonalGraphVisualizer(clonal_graph, vertex_label_writer, edge_writer)
        writer.Output(os.path.join(output_dirs['coloring_label'], graph_name))
        # compressing clonal graph
        compressor = clonal_graph_algorithms.GlonalGraphCompressor(clonal_graph)
        compressed_paths = compressor.CompressGraph()
        OutputCompressedGraph(clonal_graph, compressed_paths, os.path.join(output_dirs['compressed'], graph_name))
        # SHM visualization
#        shm_writer = shm_visualizer.SHMVisualizer(clonal_graph)
#        shm_dir = os.path.join(output_dirs['position_shms'], graph_name)
#        shm_writer.OutputSHMsPerPosition(shm_dir)
        lineage_index += 1
        if lineage_index == config.num_lineages:
            break
