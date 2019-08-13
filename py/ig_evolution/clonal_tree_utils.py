import os
import sys
import Queue
from enum import Enum

import dataset

class UndirectedClonalTree:
    def __init__(self, full_length_lineage, seqs):
        self.full_length_lineage = full_length_lineage
        self.seqs = seqs
        self.adj_list = dict()
        self.edge_weights = dict()

    def _AddVertex(self, v):
        if v not in self.adj_list:
            self.adj_list[v] = set()

    def AddEdge(self, edge, weight):
        self._AddVertex(edge[0])
        self.adj_list[edge[0]].add(edge[1])
        self._AddVertex(edge[1])
        self.adj_list[edge[1]].add(edge[0])
        self.edge_weights[edge] = weight
        self.edge_weights[(edge[1], edge[0])] = weight

    def EdgeIter(self):
        processed_edges = set()
        for v in self.adj_list:
            for w in self.adj_list[v]:
                if (w, v) not in processed_edges:
                    processed_edges.add((v, w))
                    yield (v, w)

    def VertexIter(self):
        for v in self.adj_list:
            yield v

    def SequenceIter(self):
        for seq in self.seqs:
            yield seq

    def GetSequenceByVertex(self, v):
        return self.seqs[v]

    def GetWeightByEdge(self, e):
        return self.edge_weights[e]

    def FullLengthLineage(self):
        return self.full_length_lineage

    def Dataset(self):
        return self.full_length_lineage.Dataset()

    def GetUsedSequences(self):
        return self.seqs

    def NumVertices(self):
        return len(self.adj_list)

    def GetVertexNeighs(self, v):
        return self.adj_list[v]

class SimpleRootComputer:
    def __init__(self, undirected_clonal_tree):
        self.tree = undirected_clonal_tree
        self.seqs = undirected_clonal_tree.GetUsedSequences()
        self.full_length_lineage = undirected_clonal_tree.FullLengthLineage()
        self.dataset = undirected_clonal_tree.FullLengthLineage().Dataset()

    def GetRoot(self):
        min_num_shms = sys.maxint
        root_ind = -1
        for v in self.tree.VertexIter():
            vertex_seq_id = self.tree.GetSequenceByVertex(v).id
            cur_num_shms = self.full_length_lineage.SHMDepth(vertex_seq_id)
            if cur_num_shms < min_num_shms:
                min_num_shms = cur_num_shms
                root_ind = v
        if root_ind == -1:
            print "ERROR: root was not found"
            sys.exit(1)
        return root_ind

class DirectedClonalTree:
    def __init__(self, undirected_clonal_tree, root_index):
        self.full_length_lineage = undirected_clonal_tree.FullLengthLineage()
        self.seqs = undirected_clonal_tree.GetUsedSequences()
        self.root_index = root_index
        self._InitDirectEdges(undirected_clonal_tree)

    def _InitDirectEdges(self, undirected_clonal_tree):
        edge_classifier = DirectedEdgeClassifier(self.full_length_lineage)
        self.edges = dict() # adjacency list
        self.edge_stats = dict() # edge -> edge struct
        self.edge_weights = dict() # edge -> weight
        self.vertex_parent = dict() # vertex -> parent
        queue = Queue.Queue()
        queue.put(self.root_index)
        self.vertex_parent[self.root_index] = -1
        processed_vertices = set()
        while not queue.empty():
            cur_vertex = queue.get()
            processed_vertices.add(cur_vertex)
            self.edges[cur_vertex] = set()
            adj_neighs = undirected_clonal_tree.GetVertexNeighs(cur_vertex)
            for n in adj_neighs:
                if n in processed_vertices:
                    continue
                queue.put(n)
                self.edges[cur_vertex].add(n)
                self.edge_stats[(cur_vertex, n)] = edge_classifier.ClassifyEdge(undirected_clonal_tree.GetSequenceByVertex(cur_vertex).id, undirected_clonal_tree.GetSequenceByVertex(n).id) 
                self.edge_weights[(cur_vertex, n)] = undirected_clonal_tree.GetWeightByEdge((cur_vertex, n))
                self.vertex_parent[n] = cur_vertex
        for edge in self.edge_stats:
            if edge[1] not in self.edges:
                self.edges[edge[1]] = set()

    def RootIndex(self):
        return self.root_index

    def RootSeq(self):
        return self.seqs[self.root_index]

    def IsRoot(self, v):
        return self.RootIndex() == v

    def IsLeaf(self, v):
        return len(self.edges[v]) == 0

    def VertexIter(self):
        for v in self.edges:
            yield v

    def GetVertexNeighs(self, v):
        return self.edges[v]

    def GetParent(self, v):
        return self.vertex_parent[v]

    def EdgeIter(self):
        for v1 in self.edges:
            for v2 in self.edges[v1]: 
                yield (v1, v2)

    def GetStatsByEdge(self, edge):
        return self.edge_stats[edge]

    def GetWeightByEdge(self, e):
        return self.edge_weights[e]

    def GetSequenceByVertex(self, v):
        return self.seqs[v]

    def SequenceIter(self):
        for seq in self.seqs:
            yield seq

    def FullLengthLineage(self):
        return self.full_length_lineage

    def Dataset(self):
        return self.full_length_lineage.Dataset()

    def GetUsedSequences(self):
        return self.seqs

    def NumVertices(self):
        return len(self.edges)

    ################ editing methods
    def RemoveEdge(self, edge):
        self.edge_stats.pop(edge)
        self.edge_weights.pop(edge)

    def RemoveVertex(self, v):
        if not self.IsRoot(v):
            parent = self.GetParent(v)
            self.RemoveEdge((parent, v))
            self.edges[parent].remove(v)
            self.vertex_parent.pop(v)
        else:
            # old root is removed, all its descendants become new roots
            if not self.IsLeaf(v):
                self.root_index = list(self.edges[v])[0] # currently we mark only the first descendant as a new root
        descendants = self.edges[v]
        for w in descendants:
            self.RemoveEdge((v, w))
            self.vertex_parent[w] = -1
#            self.edges[v].remove(w)
        self.edges.pop(v)

class SHMAnalyzer:
    def __init__(self, full_length_lineage):
        self.full_length_lineage = full_length_lineage
        self.dataset = self.full_length_lineage.Dataset()

    def _GetVertexSHMs(self, seq_id):
        v_shms = self.dataset.GetSHMsBySeqName(seq_id, dataset.AnnotatedGene.V)
        j_shms = self.dataset.GetSHMsBySeqName(seq_id, dataset.AnnotatedGene.J)
        cdr3_bounds = self.dataset.GetCDR3BoundsBySeqName(seq_id)
        return [shm for shm in v_shms if shm.read_pos < cdr3_bounds[0]], [shm for shm in j_shms if shm.read_pos > cdr3_bounds[1]]

    def _ComputeSharedSHMs(self, shms_1, shms_2):
        shm_1_set = set(shms_1)
        shared_shms = []
        for shm in shms_2:
            if shm in shm_1_set:
                shared_shms.append(shm)
        return shared_shms 

    def _GetReverseAddedList(self, src_shms, dst_shms):
        src_shm_set = set(src_shms) 
        dst_shm_set = set(dst_shms)
        reverse_shms = []
        added_shms = []
        for shm in src_shms:
            if shm not in dst_shm_set:
                reverse_shms.append(shm)
        for shm in dst_shms:
            if shm not in src_shm_set:
                added_shms.append(shm)
        return reverse_shms, added_shms

    def CompareSHMs(self, src_id, dst_id):
#        print "==== "
#        print src_id
#        print dst_id
        shms_1 = self._GetVertexSHMs(src_id)
        shms_2 = self._GetVertexSHMs(dst_id)
        reverse_v_shms, added_v_shms = self._GetReverseAddedList(shms_1[0], shms_2[0])
#        print "reverse V: " + str(reverse_v_shms)
#        print "added V: " + str(added_v_shms)
        reverse_j_shms, added_j_shms = self._GetReverseAddedList(shms_1[1], shms_2[1])
#        print "reverse J: " + str(reverse_j_shms)
#        print "added J: " + str(added_j_shms)
        return reverse_v_shms, added_v_shms, reverse_j_shms, added_j_shms

class DirectedEdgeType(Enum):
    DIRECTED = 0
    CDR3 = 1
    REVERSE = 2
    SIBLING = 3
    MIXED_ALIGNMENT = 4

class DirectEdge:
    def __init__(self, edge_type, reversed_v, added_v, reversed_j, added_j, cdr3_shms):
        self.edge_type = edge_type
        # V SHMs
        self.reversed_v = reversed_v
        self.added_v = added_v
        # J SHMs
        self.reversed_j = reversed_j
        self.added_j = added_j
        # CDR3
        self.cdr3_shms = cdr3_shms
#        self._SanityCheck()

    def _SanityCheck(self):
        if len(self.reversed_v) == 0 and len(self.added_v) == 0 and len(self.reversed_j) == 0 and len(self.added_j) == 0 and len(self.cdr3_shms) == 0:
            print "ERROR: two identical reads in the sample"
            sys.exit(1)

    def VSHMIter(self):
        num_v_shms = len(self.reversed_v) + len(self.added_v)
        for i in range(num_v_shms):
            if i < len(self.reversed_v):
                yield self.reversed_v[i]
            else:
                yield self.added_v[i - len(self.reversed_v)]

    def JSHMIter(self):
        num_j_shms = len(self.reversed_j) + len(self.added_j)
        for i in range(num_j_shms):
            if i < len(self.reversed_j):
                yield self.reversed_j[i]
            else:
                yield self.added_j[i - len(self.reversed_j)]

    def CDR3SHMIter(self):
        for shm in self.cdr3_shms:
            yield shm

    def NumReversedV(self):
        return len(self.reversed_v)

    def NumAddedV(self):
        return len(self.added_v)

    def NumReversedJ(self):
        return len(self.reversed_j)

    def NumAddedJ(self):
        return len(self.added_j)

    def NumCDR3SHMs(self):
        return len(self.cdr3_shms)

class DirectedEdgeClassifier:
    def __init__(self, full_length_lineage):
        self.full_length_lineage = full_length_lineage
        self.shm_analyzer = SHMAnalyzer(self.full_length_lineage)

    def _SiblingEdgeIsDoubleMutated(self, reverse_v_shms, reverse_j_shms, added_v_shms, added_j_shms):
        added_v_positions = set([shm.pos for shm in added_v_shms])
        added_j_positions = set([shm.pos for shm in added_j_shms])
        for shm in reverse_v_shms:
            if shm.pos not in added_v_positions:
                return False
        for shm in reverse_j_shms:
            if shm.pos not in added_j_positions:
                return False
        return True

    def _GetEdgeType(self, reverse_v_shms, reverse_j_shms, added_v_shms, added_j_shms, cdr3_shms):
        num_reverse = len(reverse_v_shms) + len(reverse_j_shms)
        num_added = len(added_v_shms) + len(added_j_shms)
        if num_reverse + num_added + len(cdr3_shms) == 0:
            return DirectedEdgeType.MIXED_ALIGNMENT
        if num_reverse == 0 and num_added == 0:
            return DirectedEdgeType.CDR3
        if num_reverse == 0:
            return DirectedEdgeType.DIRECTED
        if num_added == 0:
            return DirectedEdgeType.REVERSE
        if self._SiblingEdgeIsDoubleMutated(reverse_v_shms, reverse_j_shms, added_v_shms, added_j_shms):
            return DirectedEdgeType.DIRECTED
        if self._SiblingEdgeIsDoubleMutated(added_v_shms, added_j_shms, reverse_v_shms, reverse_j_shms):
            return DirectedEdgeType.REVERSE
        return DirectedEdgeType.SIBLING

    def _GetCDR3SHMs(self, src_id, dst_id):
        cdr3_src = self.full_length_lineage.Dataset().GetCDR3BySeqName(src_id)
        cdr3_dst = self.full_length_lineage.Dataset().GetCDR3BySeqName(dst_id)
        cdr3_shms = []
        for i in range(len(cdr3_src)):
            if cdr3_src[i] != cdr3_dst[i]:
                cdr3_shms.append(dataset.SHM(i, -1, cdr3_src[i], cdr3_dst[i]))
        return cdr3_shms
    
    def ClassifyEdge(self, src_id, dst_id):
        reverse_v_shms, added_v_shms, reverse_j_shms, added_j_shms = self.shm_analyzer.CompareSHMs(src_id, dst_id)
        cdr3_shms = self._GetCDR3SHMs(src_id, dst_id)
        edge_type = self._GetEdgeType(reverse_v_shms, reverse_j_shms, added_v_shms, added_j_shms, cdr3_shms)
        return DirectEdge(edge_type, reverse_v_shms, added_v_shms, reverse_j_shms, added_j_shms, cdr3_shms)
