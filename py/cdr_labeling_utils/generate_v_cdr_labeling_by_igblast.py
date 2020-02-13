import os
import sys

class Block:
    def __init__(self, query):
        self.query = query
        self.region_dict = dict()
        self.regions = ['FR1', 'CDR1', 'FR2', 'CDR2', 'FR3']

    def AddRegion(self, region, start, end):
        self.region_dict[region] = (start, end)

    def Empty(self):
        return len(self.region_dict) == 0

    def __iter__(self):
        for r in self.regions:
            yield r

    def GetBoundsByRegion(self, r):
        return self.region_dict[r]

    def __repr__(self):
        return self.query + ': ' + str(self.region_dict)

def StringHasPrefix(str, prefix):
    return str[ : len(prefix)] == prefix

def LineIsQuery(line):
    return StringHasPrefix(line, '# Query:')

input_txt = sys.argv[1]
output_fname = sys.argv[2]

blocks = []
region_starts = ['FR1-IMGT', 'CDR1-IMGT', 'FR2-IMGT', 'CDR2-IMGT', 'FR3-IMGT']
cdr3_start = 'CDR3-IMGT (germline)'
cur_block = Block('')
for l in open(input_txt).readlines():
    l = l.strip()
    if LineIsQuery(l):
        if not cur_block.Empty():
            blocks.append(cur_block)
        gene_name = l.split()[2]
        cur_block = Block(gene_name)
        continue
    for r in region_starts:
        if StringHasPrefix(l, r):
            splits = l.split()
            cur_block.AddRegion(r.split('-')[0], int(splits[1]), int(splits[2]))
blocks.append(cur_block)

for b in blocks:
    print b
print len(blocks)

output_fh = open(output_fname, 'w')
for b in blocks:
    output_fh.write(b.query + '\t')
    for r in b:
        bounds = b.GetBoundsByRegion(r)
        output_fh.write(str(bounds[0]) + '\t' + str(bounds[1]) + '\t')
    output_fh.write('VH\t0\n')
output_fh.close()
