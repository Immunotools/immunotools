import os
import sys

class ClonalTreeStatsWriter:
    def __init__(self, directed_clonal_tree):
        self.directed_clonal_tree = directed_clonal_tree

    def Output(self, fname):
        fh = open(fname, 'w')
#        fh.write('Vertices:\n')
#        for v in self.directed_clonal_tree.VertexIter():
#            fh.write(self.directed_clonal_tree.GetSequenceByVertex(v).id + '\n')
#        fh.write('Edges:\n')
#        for e in self.directed_clonal_tree.EdgeIter():
#            src_id = self.directed_clonal_tree.GetSequenceByVertex(e[0]).id
#            dst_id = self.directed_clonal_tree.GetSequenceByVertex(e[1]).id
#            fh.write(src_id + ' -> ' + dst_id + '\n')
        fh.write('Src_ID\tDst_ID\tSrc_mult\tDst_mult\tEdge_type\tNum_reversed_v\tNum_added_v\tNum_reversed_j\tNum_added_j\tnum_CDR3\n')
        for e in self.directed_clonal_tree.EdgeIter():
            edge_stats = self.directed_clonal_tree.GetStatsByEdge(e)
            src_id = self.directed_clonal_tree.GetSequenceByVertex(e[0]).id
            dst_id = self.directed_clonal_tree.GetSequenceByVertex(e[1]).id
            src_mult = self.directed_clonal_tree.Dataset().GetSeqMultiplicity(src_id)
            dst_mult = self.directed_clonal_tree.Dataset().GetSeqMultiplicity(dst_id)
            fh.write(src_id + '\t' + dst_id + '\t' + str(src_mult) + '\t' + str(dst_mult) + '\t' + edge_stats.edge_type.name + '\t' + str(edge_stats.NumReversedV()) + '\t' + str(edge_stats.NumAddedV()) + '\t' + str(edge_stats.NumReversedJ()) + '\t' + str(edge_stats.NumAddedJ()) + '\t' + str(edge_stats.NumCDR3SHMs()) + '\n')
        if self.directed_clonal_tree.NumVertices() == 1:
            v_id = ''
            for v in self.directed_clonal_tree.VertexIter():
                v_id = self.directed_clonal_tree.GetSequenceByVertex(v).id
            v_mult = self.directed_clonal_tree.Dataset().GetSeqMultiplicity(v_id)
            fh.write(v_id + '\tNA\t' + str(v_mult) + '\tNA\tNA\tNA\tNA\tNA\tNA\tNA\n')
        fh.close()
        
