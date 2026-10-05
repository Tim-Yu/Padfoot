#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 13 17:39:44 2025

@author: keskusa2
"""

from collections import defaultdict
import gzip

def open_text(path):
    """Open a plain or gzip/bgzip-compressed text file for reading (told apart by the gzip magic bytes, not the name)."""
    with open(path, 'rb') as fh:
        magic = fh.read(2)
    return gzip.open(path, 'rt') if magic == b'\x1f\x8b' else open(path, 'r')

def generate_gff(new_gff, gff_name):
    """GENCODE-style GFF3 (plain or gzipped) -> Padfoot's 6-column gene table (chr gene strand type start end), the format
    of the bundled beds/<genome>.gff3.gz that get_genes() reads; gzip-compressed when gff_name ends with .gz. Returns the
    number of protein-coding genes written."""
    fopen = open_text(new_gff)
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
    fopen.close()

    out_file3 = gff_name
    with (gzip.open if out_file3.endswith('.gz') else open)(out_file3, "wt") as fout3:
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
    return len(genes)

def generate_rm(rm_file, rm_path):
    """Repeat annotation (plain or gzipped) -> the 4-column BED (chr start end class/family) of the bundled
    beds/<genome>_rm.bed.gz, which annot_bp_repeat() intersects the breakpoints with. Accepts
    - RepeatMasker .out / .fa.out: whitespace-aligned, 3 header lines, an optional 16th column '*'; query = column 5,
      begin/end = columns 6/7, class/family = column 11. begin is kept as written, as in the bundled files;
    - BED with at least 4 columns: chr, start, end, class/family (further columns are dropped).
    Returns the number of repeats written."""
    n = 0
    with open_text(rm_file) as fin, open(rm_path, 'w') as fout:
        for i, line in enumerate(fin, 1):
            l = line.split()
            if not l or l[0].startswith('#') or l[0] in ('SW', 'score', 'track', 'browser'):
                continue   # blank line, comment, RepeatMasker column header, BED track/browser line
            if len(l) >= 11 and l[5].isdigit() and l[6].isdigit():
                row = [l[4], l[5], l[6], l[10]]
            elif len(l) >= 4 and l[1].isdigit() and l[2].isdigit():
                row = l[:4]
            else:
                raise ValueError(f"{rm_file} line {i} is neither a RepeatMasker .out row nor a BED row "
                                 f"(chr start end class/family): {line.strip()}")
            fout.write('\t'.join(row) + '\n')
            n += 1
    return n
    