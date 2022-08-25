import os
import sys
import queue as Queue

class DFSVertexOrderFinder: 
    def __init__(self, clonal_graph):
        self.clonal_graph = clonal_graph

    def GetOrder(self):
        root = self.clonal_graph.GetRootIndex()
        queue = Queue.PriorityQueue()
        queue.put((0, root))
        order = []
        processed_vertices = set()
        while not queue.empty():
            cur_item = queue.get()
#            print cur_item
            order.append(cur_item[1])
            processed_vertices.add(cur_item[1])
            child_vertices = self.clonal_graph.GetOutgoingVertices(cur_item[1])
            if len(child_vertices) == 1 and child_vertices[0] not in processed_vertices:
                queue.put((cur_item[0], child_vertices[0]))
                processed_vertices.add(child_vertices[0])
                continue
            child_prior = cur_item[0]
            for v in child_vertices:
                if v in processed_vertices:
                    continue
                child_prior -= 1
                queue.put((child_prior, v))
                processed_vertices.add(v)
        return order

def GetLevelsByVertexOrder(clonal_graph, vertex_order):
    levels = []
    cur_level = 1
    levels.append(cur_level)
    for i in range(0, len(vertex_order) - 1):
        child_vertices = clonal_graph.GetOutgoingVertices(vertex_order[i])
        if not(len(child_vertices) == 1 and child_vertices[0] == vertex_order[i + 1]):
            cur_level += 1
        levels.append(cur_level)
    return levels

class GlonalGraphCompressor:
    def __init__(self, clonal_graph):
        self.clonal_graph = clonal_graph

    def _GetNextVertex(self, non_visited_vertices):
        for v in non_visited_vertices:
            return v
        return -1

    def CompressGraph(self):
        non_visited_vertices = set()
        for v in self.clonal_graph.VertexIter():
            non_visited_vertices.add(v)
        start_vertex = self._GetNextVertex(non_visited_vertices)
        non_branching_paths = []
        while start_vertex != -1:
            path = self._GetNonBranchingPath(start_vertex)
            if len(path) == 0:
                if start_vertex in non_visited_vertices:
                    non_visited_vertices.remove(start_vertex)
            else:
                #print 'Non-branching path: ' + str(path)
                for v in path:
                    non_visited_vertices.remove(v)
                if len(path) > 1:
                    non_branching_paths.append(path)
            start_vertex = self._GetNextVertex(non_visited_vertices)
        print(str(len(non_branching_paths)) + ' non-trivial non-branching paths were extracted from clonal graph')
        return non_branching_paths

    def _VertexIsGood(self, v, visited_vertices): 
        return len(self.clonal_graph.GetIncomingVertices(v)) <= 1 and len(self.clonal_graph.GetOutgoingVertices(v)) <= 1 and v not in visited_vertices

    def _GetNonBranchingPath(self, start_vertex):
        path = []
        cur_vertex = start_vertex
        visited = set()
        while self._VertexIsGood(cur_vertex, visited):
            path.insert(0, cur_vertex)
            visited.add(cur_vertex)
            if len(self.clonal_graph.GetIncomingVertices(cur_vertex)) == 0:
                break
            cur_vertex = self.clonal_graph.GetIncomingVertices(cur_vertex)[0]
        if len(path) != 0:
            path = path[ : len(path) - 1]
            visited.remove(start_vertex)
        cur_vertex = start_vertex
        while self._VertexIsGood(cur_vertex, visited):
            path.append(cur_vertex)
            visited.add(cur_vertex)
            if len(self.clonal_graph.GetOutgoingVertices(cur_vertex)) == 0:
                break
            cur_vertex = self.clonal_graph.GetOutgoingVertices(cur_vertex)[0]
        return path
