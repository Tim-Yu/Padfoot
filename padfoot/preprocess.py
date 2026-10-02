#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 13 17:39:44 2025

@author: keskusa2
"""

from collections import defaultdict
import pandas as pd

def generate_gff(new_gff, gff_name):
    fopen = open(new_gff, 'r')
    gene_ls = []
    exonls = []
    old_gene_name = ''
    old_strand = ''
    old_ref_id = ''
    genes = defaultdict(list)

    def flush():
        # write the gene collected so far under ITS OWN chromosome (not the chromosome of the row that ends it)
        if gene_ls:
            gene_ls[-1].append(exonls)
            lns = [x[0][2] - x[0][1] for x in gene_ls[1:]]
            ind = lns.index(max(lns))
            genes[(old_ref_id, old_gene_name, old_strand)] = [gene_ls[0]] + gene_ls[ind+1]
    for line in fopen:
        if not 'gene_type=protein_coding' in line:
            continue
        gene_name = [l.split('=')[1] for l in line.split()[-1].split(';') if 'gene_name' in l][0]
        l = line.split()
        ref_id, typ, start, end, strand = l[0], l[2], int(l[3]), int(l[4]), l[6]
        if typ == 'CDS':
            continue
        if typ == 'gene':
            flush()
            old_gene_name = gene_name
            old_strand =  strand
            old_ref_id = ref_id
            gene_ls = [(typ,start, end)]
            exonls = []
        elif typ == 'transcript':
            if exonls:
                gene_ls[-1].append(exonls)
            exonls = []
            tr_id  = [l.split('=')[1] for l in line.split()[-1].split(';') if 'transcript_name' in l][0]
            gene_ls.append([(tr_id, start, end)])
        else:
            if typ == 'exon':
                typ = 'exon' + [l.split('=')[1] for l in line.split()[-1].split(';') if 'exon_number' in l][0]
            exonls.append((typ,start, end))
    flush()   # the last gene of the file

    out_file3 = gff_name
    with open(out_file3, "w") as fout3:
        for (ref_id, gene_name, strand), exons in genes.items():
            line = '\t'.join([ref_id, gene_name, strand, exons[0][0], str(exons[0][1]), str(exons[0][2])])
            fout3.write(line)
            fout3.write('\n')
            line = '\t'.join([ref_id, gene_name, strand, exons[1][0], str(exons[1][1]), str(exons[1][2])])
            fout3.write(line)
            fout3.write('\n')
            if len(exons) > 2:
                for exon in exons[2]:
                    line = '\t'.join([ref_id, gene_name, strand, exon[0], str(exon[1]), str(exon[2])])
                    fout3.write(line)
                    fout3.write('\n')
    fout3.close()

def generate_rm(rm_file, rm_path):
    df = pd.read_csv(rm_file, sep='\t', header=None)
    df_small = df[[4, 5, 6, 10]]
    df_small.to_csv(rm_path, sep='\t', header=False, index=False)
    