import numpy

class LowFixedAbundanceLeafRemover:
    def __init__(self, clonal_tree, min_abundance):
        self.clonal_tree = clonal_tree
        self.dataset = clonal_tree.Dataset()
        self.min_abundance = min_abundance

    def VertexToBeRemoved(self, v):
        if not self.clonal_tree.IsRoot(v) and not self.clonal_tree.IsLeaf(v):
            return False
        if self.clonal_tree.IsRoot(v) and len(self.clonal_tree.GetVertexNeighs(v)) > 1:
            return False
        vertex_seq_id = self.clonal_tree.GetSequenceByVertex(v).id
        return self.dataset.GetSeqMultiplicity(vertex_seq_id) < self.min_abundance

    def Name(self):
        return 'LowFixedAbundanceLeafRemover'

class RelativeAbundanceLeafRemover:
    def __init__(self, clonal_tree, min_rel_abundance):
        self.clonal_tree = clonal_tree
        self.dataset = clonal_tree.Dataset()
        self.min_rel_abundance = min_rel_abundance

    def _GetNeighbour(self, v):
        if self.clonal_tree.IsRoot(v):
            return list(self.clonal_tree.GetVertexNeighs(v))[0]
        return self.clonal_tree.GetParent(v)

    def _GetVertexMultiplicity(self, v):
        vertex_seq_id = self.clonal_tree.GetSequenceByVertex(v).id
        return self.dataset.GetSeqMultiplicity(vertex_seq_id)

    def VertexToBeRemoved(self, v):
        if not self.clonal_tree.IsRoot(v) and not self.clonal_tree.IsLeaf(v):
            return False
        if self.clonal_tree.IsRoot(v) and len(self.clonal_tree.GetVertexNeighs(v)) != 1:
            return False
        neigh = self._GetNeighbour(v)
        #print str(self._GetVertexMultiplicity(v)) + '\t' + str(self._GetVertexMultiplicity(neigh)) + '\t' + str(float(self._GetVertexMultiplicity(v)) / self._GetVertexMultiplicity(neigh)) + '\t' + str(len(self.clonal_tree.GetVertexNeighs(neigh)))
        rel_mult = float(self._GetVertexMultiplicity(v)) / self._GetVertexMultiplicity(neigh) 
#        fh = open('relative_multiplicity.txt', 'a+')
#        fh.write(str(rel_mult) + '\n')
        return rel_mult <= self.min_rel_abundance

    def Name(self):
        return 'RelativeAbundanceLeafRemover'

class IterativeTipRemover:
    def __init__(self, directed_clonal_tree, vertex_filters):
        self.clonal_tree = directed_clonal_tree
        self.vertex_filters = vertex_filters

    def _CleanTipsInCurrentTree(self, vertex_filter):
        vertices_to_remove = []
        for v in self.clonal_tree.VertexIter():
            if len(vertices_to_remove) == self.clonal_tree.NumVertices() - 1:
                break
            if vertex_filter.VertexToBeRemoved(v):
                vertices_to_remove.append(v)
        for v in vertices_to_remove:
            self.clonal_tree.RemoveVertex(v)
        return len(vertices_to_remove)

    def CleanTips(self):
        for vertex_filter in self.vertex_filters:
#            print "== Applying " + vertex_filter.Name()
            num_removed_vertices = 1
            iteration_step = 1
            while num_removed_vertices != 0: 
                num_removed_vertices = self._CleanTipsInCurrentTree(vertex_filter)
#                print 'Iteration ' + str(iteration_step) + ': ' + str(num_removed_vertices) + ' vertices were removed'
                iteration_step += 1
        return self.clonal_tree

    
