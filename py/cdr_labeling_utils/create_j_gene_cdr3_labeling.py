import os
import sys
from Bio import SeqIO

gene_fasta = sys.argv[1]
output_fname = sys.argv[2]
chain = sys.argv[3]

chain_dict = {'IGH' : ['TGGGGC'], 'IGK' : ['TTAGG', 'TTAGG', 'TTCGG', 'TTTGG'], 'IGL' : ['TTAGG', 'TTCGG', 'TTCAT'], 'TRA' : ['TTTGG', 'TTCGG', 'TTCAG'], 'TRB' : ['TTTGG', 'TTCGG', 'TTCAG'], 'TRG' : ['TTTGG', 'TTCGG', 'TTGGA']}
chain_id_dict = {'IGH' : 'HJ', 'IGK' : 'KJ', 'IGL' : 'LJ', 'TRA' : 'AJ', 'TRB' : 'BJ', 'TRG' : 'GJ', 'TRD' : 'DJ'}

aligned_genes = []
for r in SeqIO.parse(gene_fasta, 'fasta'):
    r.seq = str(r.seq) #.upper()
    aligned_genes.append(r)

output_fh = open(output_fname, 'w')
for gene in aligned_genes:
    w_pos = -1
    for cons_seq in chain_dict[chain]:
        w_pos = gene.seq.upper().find(cons_seq)
        if w_pos != -1:
            break
    if w_pos == -1:
        print('WARN: ' + gene.id + ' motif position was not found')
        continue
    output_fh.write(gene.id + '\t' + str(w_pos) + '\t' + str(w_pos + 1) + '\t' + chain_id_dict[chain] + '\n')
output_fh.close()
