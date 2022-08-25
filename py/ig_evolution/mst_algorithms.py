import os
import sys

import queue as Queue

import igraph
from igraph import *

import disjoint_set

def GetVerticesByEdgeDict(edge_dict):
    vertices = set()
    for e in edge_dict:
        vertices.add(e[0])
        vertices.add(e[1])
    return vertices

# contract: ComputeSpanningTree(sequences, edge_dict) -> tree_edges
class IGraphMSTFinder:
    def ComputeSpanningTree(self, sequences, edge_dict):
        graph = Graph()
        num_vertices = max(GetVerticesByEdgeDict(edge_dict)) + 1 #max(len(GetVerticesByEdgeDict(edge_dict)), max(GetVerticesByEdgeDict(edge_dict)))
        graph.add_vertices(num_vertices)
#        for e in edge_dict:
#	    if e[0] >= num_vertices or e[1] >= num_vertices:
#                print(e)
        graph.add_edges(edge_dict.keys())
        spanning_tree = graph.spanning_tree(weights = edge_dict.values(), return_tree = True)
        tree_weights = dict()
        for e in spanning_tree.get_edgelist():
            tree_weights[e] = edge_dict[e]
        return tree_weights

class VertexMultMSTFinder:
    def __init__(self, full_length_lineage):
        self.full_length_lineage = full_length_lineage
        self.dataset = self.full_length_lineage.Dataset()
     
    def ComputeSpanningTree(self, sequences, edge_dict):
#        print "==== "
#        print "Original graph: " + str(edge_dict)
        mult_dict = dict()
        vertices = GetVerticesByEdgeDict(edge_dict)
        for v in vertices:
            seq_id = sequences[v].id
            mult_dict[v] = self.dataset.GetSeqMultiplicity(seq_id)
        tree_edges = self._ComputeTree(edge_dict, mult_dict, vertices)
        if len(tree_edges) + 1 != len(vertices) and len(tree_edges) != 0 and len(vertices) != 0:
            print("ERROR: Tree structure is not correct! # edges: " + str(len(tree_edges)) + ', # vertices: ' + str(len(vertices)))
            sys.exit(1)
#        print 'Tree: ' + str(len(tree_edges)) + ', vertices: ' + str(len(vertices))
        return tree_edges

    def _ComputeTree(self, edge_dict, vertex_mult_dict, vertices):
        weight_queues = self._ComputeWeightDegreeQueue(edge_dict, vertex_mult_dict)
        tree_edges = dict()
        num_processed_edges = 0
        vertex_ds = disjoint_set.DisjointSet(vertices)
        for w in sorted(weight_queues.keys()):
            curr_queue = weight_queues[w]
            while not curr_queue.empty():
                edge_pair = curr_queue.get() 
                weight = edge_pair[0]
                edge = edge_pair[1]
                if vertex_ds.find_index(edge[0]) == vertex_ds.find_index(edge[1]):
                    continue
                num_processed_edges += 1
                tree_edges[edge] = edge_dict[edge]
                vertex_ds.union(edge[0], edge[1])
#        print "# processed: " + str(num_processed_edges) + ' # original edges: ' + str(len(edge_dict))
        return tree_edges

    def _ComputeEdgeDegreeDict(self, edge_dict, vertex_mult_dict):
        edge_degree_dict = dict()
        for e in edge_dict:
            edge_degree_dict[e] = max(vertex_mult_dict[e[0]], vertex_mult_dict[e[1]])
        return edge_degree_dict

    def _ComputeMaxDegreeForEachWeight(self, edge_dict, edge_degree_dict):
        weight_max_degree_dict = dict()
        for e in edge_dict:
            w = edge_dict[e]
            if w not in weight_max_degree_dict:
                weight_max_degree_dict[w] = edge_degree_dict[e]
            weight_max_degree_dict[w] = max(weight_max_degree_dict[w], edge_degree_dict[e])
        return weight_max_degree_dict

    def _ComputeWeightDegreeQueue(self, edge_dict, vertex_mult_dict):
         edge_degree_dict = self._ComputeEdgeDegreeDict(edge_dict, vertex_mult_dict) # edge -> edge_degree
         weight_max_degree_dict = self._ComputeMaxDegreeForEachWeight(edge_dict, edge_degree_dict) # weight -> max edge degree
         weight_degree_queue = dict() # weight -> queue (max_weight_degree - edge_degree, edge) 
         for e in edge_dict:
             w = edge_dict[e]
             if w not in weight_degree_queue:
                 weight_degree_queue[w] = Queue.PriorityQueue()
             edge_priority = weight_max_degree_dict[w] - edge_degree_dict[e]
             weight_degree_queue[w].put((edge_priority, e))
         return weight_degree_queue
