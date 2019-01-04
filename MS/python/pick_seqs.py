#!/usr/bin/env python
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import os, sys

fst = sys.argv[1]
lst = sys.argv[2]

gene_lst = []
lines = open( lst, "r" ).readlines()
for l in lines:
    gene_lst.append(l.split()[0])

handle = open( fst, "rU" )
picked_handle = open( "sample_all-match.fasta", "w" )
for record in SeqIO.parse(handle, "fasta"):
    if record.id in gene_lst:
        SeqIO.write( record, picked_handle, "fasta" )

handle.close()
picked_handle.close()
