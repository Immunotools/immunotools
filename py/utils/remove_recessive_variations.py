import os
import sys
from Bio import SeqIO

input_fasta = sys.argv[1]
output_fasta = sys.argv[2]

output_fh = open(output_fasta, 'w')
for r in SeqIO.parse(input_fasta, 'fasta'):
    if int(r.id.split('*')[1]) != 1:
        continue
    output_fh.write('>' + r.id + '\n')
    output_fh.write(str(r.seq) + '\n')
output_fh.close()
