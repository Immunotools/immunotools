import os

class SHMVisualizer:
    def __init__(self, clonal_graph):
        self.clonal_graph = clonal_graph
        self.shms = clonal_graph.GetSHMs()

    def OutputSHMsPerPosition(self, output_dir):
        if not os.path.exists(output_dir):
            os.mkdir(output_dir)
        self._InitSHMPositionMap()
        for pos in self.pos_shm_dict:
            if len(self.pos_shm_dict[pos]) == 1:
                continue
            output_base = os.path.join(output_dir, self._GetFnameForSHMPosition(pos))
            self._OutputSHMGraph(pos, output_base) 

    def _SHMIsToStopCodon(self, shm):
        return shm.dst_n == '*' or shm.src_n == '*'

    def _InitSHMPositionMap(self):
        self.pos_shm_dict = dict()
        for shm in self.shms.SHMIter():
            if self._SHMIsToStopCodon(shm):
                continue
            if shm.pos not in self.pos_shm_dict:
                self.pos_shm_dict[shm.pos] = []
            self.pos_shm_dict[shm.pos].append(shm)

    def _GetFnameForSHMPosition(self, pos):
        return 'pos' + str(pos) + '_' + self.shms.GetRegionForSHM(self.pos_shm_dict[pos][0]).name
    
    def _OutputSHMGraph(self, pos, output_base):
        output_dot = output_base + '.dot'
        fh = open(output_dot, 'w')
        fh.write('digraph{\n')
        for shm in self.pos_shm_dict[pos]:
            fh.write(shm.src_n + ' -> ' + shm.dst_n + ' [label = \"' + str(self.shms.GetSHMMultiplicity(shm)) + '\"]\n')
        fh.write('}')
        fh.close() 
        os.system('dot -Tpdf ' + output_dot + ' -o ' + output_base + '.pdf')
#        os.system('dot -Tsvg ' + output_dot + ' -o ' + output_base + '.svg')
