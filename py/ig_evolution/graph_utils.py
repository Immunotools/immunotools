import os
import sys
import shutil
import Queue
from Bio import SeqIO
import operator

class Graph:
    def _GetStartVertex(self):
        for v in self.non_visited:
            return v
        return -1

    def _GetConnectedComponentByStartVertex(self, start_v):
        component = []
        visited = set()
        visited.add(start_v)
        queue = Queue.Queue()
        queue.put(start_v)
        while not queue.empty():
            cur_v = queue.get()
            component.append(cur_v)
            for v in self.edges[cur_v]:
                if v not in visited:
                    queue.put(v)
                    visited.add(v)
        return component

    def _InitEdgeDict(self):
        self.edge_dict = dict() # (v, w) -> weight
        for v in self.edges:
            neighs = self.edges[v]
            weights = self.edge_weights[v]
            for i in range(len(neighs)):
                reverse_edge = (neighs[i], v)
                if reverse_edge in self.edge_dict:
                    continue
                self.edge_dict[(v, neighs[i])] = weights[i]

    def __init__(self, fname):
        self.edges = dict()
        self.edge_weights = dict()
        self.num_edges = 0
        lines = open(fname, 'r').readlines()
        self.num_vertices = len(lines) - 1
        for i in range(1, len(lines)):
            self.edges[i - 1] = []
            self.edge_weights[i - 1] = []
            splits = lines[i].strip().split()
            for j in range(0, len(splits) / 2):
                self.edges[i - 1].append(int(splits[j * 2]) - 1)
                self.edge_weights[i - 1].append(int(splits[j * 2 + 1]))
                self.num_edges += 1
        self.num_edges = self.num_edges / 2
        self._InitEdgeDict()

    def EdgeIterator(self):
        for e in self.edge_dict:
            yield e

    def GetEdgeWeight(self, e):
        return self.edge_dict[e]

    def GetConnectedComponents(self):
        connected_list = []
        self.non_visited = set(range(0, len(self.edges))) #[0] * len(self.edges)
        start_v = self._GetStartVertex()
        while start_v != -1:
            component = self._GetConnectedComponentByStartVertex(start_v)
            connected_list.append(component)
            for c in component:
                self.non_visited.remove(c)
            start_v = self._GetStartVertex()
        return connected_list

    def NumVertices(self):
        return self.num_vertices

    def NumEdges(self):
        return self.num_edges
