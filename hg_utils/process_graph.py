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
#        for i in range(0, len(self.visited)):
#            if self.visited[i] == 0:
#                return i
#        return -1

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
        print "Graph on " + str(len(self.edges)) + " vertices & " + str(self.num_edges) + ' edges was extracted from ' + fname

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
        print str(len(connected_list)) + " connected components were computed for the graph"
        return connected_list

###################################
def GetNumberAlignedCDR3s(component, aligned_cdr3s):
    number_aligned = 0
    for ind in component:
        if ind in aligned_cdr3s:
            number_aligned += 1
    return number_aligned

def OutputGraphComponent(graph, component, aligned_cdr3s, output_fname):
    fh = open(output_fname, 'w')
    fh.write('digraph{\n')
    for v in component:
        color = 'white'
        if v in aligned_cdr3s:
            color = 'red'
        fh.write(str(v) + ' [fillcolor = ' + color + ', style = filled];\n')
    processed_edges = set()
    for v in component:
        for i in range(len(graph.edges[v])):
            v2 = graph.edges[v][i]
            weight = graph.edge_weights[v][i]
            if (v2, v) in processed_edges:
                continue
            fh.write(str(v) + ' -> ' + str(v2) + ' [label = ' + str(weight) + '];\n')
            processed_edges.add((v, v2))
    fh.write('}')
    fh.close()

###################################
graph_fname = sys.argv[1]
#aligned_cdr3s = sys.argv[2]
output_dir = sys.argv[2]

graph = Graph(graph_fname)
components = graph.GetConnectedComponents()
max_size = 0
for c in components:
    max_size = max(max_size, len(c))
print "Size of the maximal component: " + str(max_size)

aligned_indices = set() #set(range(graph.num_vertices))
#if aligned_cdr3s != '-':
#    aligned_indices = set([(int(l.strip().split()[0]) - 1) for l in open(aligned_cdr3s).readlines() if l != ''])

num_aligned_comps = 0
num_large = 0
num_components = 0
output_stats = open(os.path.join(output_dir, "graph_stats.txt"), 'w')
for c in components:
    num_aligned_cdr3 = GetNumberAlignedCDR3s(c, aligned_indices)
    output_stats.write(str(len(c)) + '\t' + str(num_aligned_cdr3) + '\n')
    if num_aligned_cdr3 == 0:
        num_aligned_comps += 1
    if len(c) >= 10:
        num_large += 1
    num_components += 1
    OutputGraphComponent(graph, c, aligned_indices, os.path.join(output_dir, "component" + str(num_components) + '_nodes' + str(len(c)) + '.dot'))
output_stats.close()
print "Graph contains " + str(num_large) + ' large (>= 10) components'
print str(num_aligned_comps) + ' out of ' + str(num_components) + ' contain at least one aligned CDR3s'
