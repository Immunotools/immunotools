import os
import sys
import shutil
from Bio import SeqIO

def ReadFasta(fasta_fname):
    records = []
    if not os.path.exists(fasta_fname) or fasta_fname == '':
        print "ERROR: file " + fasta_fname + ' was not found'
        sys.exit(1)
    for r in SeqIO.parse(fasta_fname, 'fasta'):
        r.seq = str(r.seq).upper()
        records.append(r)
    return records

def PrepareOutputDir(output_dir):
    if os.path.exists(output_dir):
        shutil.rmtree(output_dir)
    os.mkdir(output_dir)

def CollapseIdenticalSequences(seqs):
    added_seqs = set()
    collapsed_seqs = []
    for s in seqs:
        if s.seq in added_seqs:
            continue
        collapsed_seqs.append(s)
        added_seqs.add(s.seq)
    return collapsed_seqs
