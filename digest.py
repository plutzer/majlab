#!/usr/bin/env python
import sys
from Bio import SeqIO
from collections import defaultdict
from pyopenms import *

input_fasta = sys.argv[1]

dig = ProteaseDigestion()
dig.setMissedCleavages(0) # I would set it to 0 first because missed cleavages are rarer

pep2prots = defaultdict(set)
prot2peps = defaultdict(set)

with open(input_fasta, "r") as handle:
    for record in SeqIO.parse(input_fasta, "fasta"):
        peptides = []
        seq = AASequence.fromString(str(record.seq))
        dig.digest(seq, peptides)
        for pep in peptides:
            s = pep.toString().decode("utf-8")
            if "X" not in s and len(s) >= 6 and pep.getMonoWeight() <= 4600:
                pep2prots[s].add(record.id)
                prot2peps[acc].add(s)
                
                


        