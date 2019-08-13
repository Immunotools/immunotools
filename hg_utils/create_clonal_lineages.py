import os
import sys
import shutil
from Bio import SeqIO

def ReadFasta(fasta_fname):
    records = []
    for r in SeqIO.parse(fasta_fname, 'fasta'):
        r.seq = str(r.seq).upper()
        records.append(r)
    return records

class AlignedCDR3s:
    def __init__(self, file_line):
        splits = file_line.strip().split()
        self.cdr3_index = int(splits[0]) - 1
        self.cdr3_name = splits[1]
        self.d_name = splits[2]

def ReadAlignedCDR3s(cdr3_fname):
    aligned_cdr3s = []
    lines = open(cdr3_fname).readlines()
    for l in lines:
        aligned_cdr3s.append(AlignedCDR3s(l))
    return aligned_cdr3s

def ReadGraphVertices(graph_fname):
    lines = open(graph_fname).readlines()
    vertices = []
    for i in range(1, len(lines)):
        l = lines[i].strip()
        if l.find('->') != -1 or l == '}':
            break
        splits = l.split()
        vertices.append(int(splits[0]))
    return vertices

###########################
cdr3_fasta = sys.argv[1]
graph_dir = sys.argv[2]
output_dir = sys.argv[3]

if os.path.exists(output_dir):
    shutil.rmtree(output_dir)
os.mkdir(output_dir)

cdr3s = ReadFasta(cdr3_fasta)
print str(len(cdr3s)) + ' CDR3s were extracted from ' + cdr3_fasta

graph_files = os.listdir(graph_dir)
print str(len(graph_files)) + ' graph files were found in ' + graph_dir
lineage_index = 1
cdr3_lineage_map = dict()
for g in graph_files:
    if g.find('dot') == -1:
        continue
    full_graph_fname = os.path.join(graph_dir, g)
#    print "Reading " + g + '...'
    vertices = ReadGraphVertices(full_graph_fname)
#    print vertices
    max_len = max([len(cdr3s[v].seq) for v in vertices])
    lineage_name =  g.split('.')[0] #"lineage_" + str(lineage_index) + '_size_' + str(len(vertices)) + '_len_' + str(max_len)
    output_fname = os.path.join(output_dir, lineage_name + '.fasta')
#    print "Writing " + output_fname
    output_fh = open(output_fname, 'w')
    for v in vertices:
        output_fh.write('>' + cdr3s[v].id + '\n')
        output_fh.write(cdr3s[v].seq + '\n')
        cdr3_lineage_map[cdr3s[v].id] = lineage_name
    output_fh.close()
    lineage_index += 1

output_map = os.path.join(output_dir, "cdr3_lineage_map.txt")
fh = open(output_map, 'w')
for cdr3 in cdr3_lineage_map:
    fh.write(cdr3 + '\t' + cdr3_lineage_map[cdr3] + '\n')
fh.close()
