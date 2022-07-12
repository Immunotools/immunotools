import os
import sys
from Bio import SeqIO

gene_fasta = sys.argv[1]
output_fname = sys.argv[2]
cons_seq = 'tggggc'
if len(sys.argv) == 4:
    cons_seq = sys.argv[3]

aligned_genes = []
for r in SeqIO.parse(gene_fasta, 'fasta'):
    r.seq = str(r.seq) #.upper()
    aligned_genes.append(r)

output_fh = open(output_fname, 'w')
for gene in aligned_genes:
    w_pos = gene.seq.lower().find(cons_seq)
    output_fh.write(gene.id + '\t' + str(w_pos) + '\t' + str(w_pos + 1) + '\n')
output_fh.close()
