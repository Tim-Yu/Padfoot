#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 17 14:30:38 2025

@author: keskusa2
"""

import pysam
from collections import defaultdict, Counter
from typing import Dict, List
import bisect
import subprocess
from Bio import Align
import os
import copy
import re
import numpy as np
import pandas as pd
import logging

from padfoot.preprocess import open_text

logger = logging.getLogger()

BND_MATE_RE = re.compile(r"[\[\]]([^:\[\]]+):(\d+)[\[\]]")

class SV(object):
    __slots__ = ("bp_1", "direction_1", "bp_2", "direction_2", "supp",'supp_read_ids', 'has_ins', 'sv_type','vaf','loh', 'prec',\
                 'vcf_id', 'is_single', 'vntr', 'cluster_id', 'detailed_type', 'ins_seq', 'genes', 'repeat', 'microh', 'ins_align',\
                     'telomere', 'repeat_bp', 'cn_assigned', 'cn_altering', 'hp1', 'hp2', 'score', 'cancer', 'impact')
    def __init__(self, bp_1, direction_1, bp_2, direction_2, supp, has_ins, sv_type,vaf, vcf_id, is_single, vntr, ins_seq, cluster_id, detailed_type, hp1, hp2):
        self.bp_1 = bp_1
        self.direction_1 = direction_1
        self.bp_2 = bp_2
        self.direction_2 = direction_2
        self.supp = supp
        self.supp_read_ids = ''
        self.ins_seq = ins_seq
        self.has_ins = has_ins
        self.hp1 = hp1
        self.hp2 = hp2
        self.vaf = vaf
        self.loh = ''
        self.prec = ''
        self.vcf_id = vcf_id
        self.sv_type = sv_type
        self.is_single = is_single
        self.vntr = vntr
        self.cluster_id = cluster_id
        self.detailed_type = detailed_type
        self.genes = []
        self.repeat = []
        self.microh = []
        self.ins_align = []
        self.telomere = ''
        self.repeat_bp = ['', '']
        self.cn_assigned  = ''
        self.cn_altering = False
        self.score = 0
        self.cancer = [' ',' ']
        self.impact = ' '
        
    def get_jun(self):
         f1 = 'T' if self.direction_1 == '-' else 'H'
         f2 = 'T' if self.direction_2 == '-' else 'H'
         return f1+f2
     
    def get_name(self):
        label_1 = "{0}:{1}".format(self.bp_1[0], self.bp_1[1])
        label_2 = "{0}:{1}".format(self.bp_2[0], self.bp_2[1])
        return label_1 + "|" + label_2
    
    def to_str(self):
        gg = self.genes[0] if self.genes[0] else [' ', ' ', ' ']
        gene1 = '\t'.join(gg)
        
        gg = self.genes[1] if self.genes[1] else [' ', ' ', ' ']
        gene2 = '\t'.join(gg)
        
        typ = self.genes[2] if len(self.genes)==3 else ' '
        microh = self.microh[0] if self.microh else ' '
        tel = str(self.telomere) if self.telomere else ' '
        ins_al = ';'.join(set(self.ins_align)) if self.ins_align else ' '
        rep = ','.join([':'.join([s[0], str(s[1])]) for s in self.repeat])

        st = self.direction_1 + self.direction_2
        return '\t'.join([self.vcf_id, self.bp_1[0], str(self.bp_1[1]), self.bp_2[0], str(self.bp_2[1]), str(self.supp),\
         str(round(self.vaf, 3)), st, str(self.score), self.impact, self.cancer[0], gene1,  self.repeat_bp[0], self.cancer[1], gene2, self.repeat_bp[1], typ, microh, str(self.vntr), tel, ins_al, rep])

    def to_sum(self):
        if self.bp_1[0] == self.bp_2[0]:
            pos = self.bp_1[0] + ':' + str(self.bp_1[1]) + '-' + str(self.bp_2[1])
        else:
            pos = self.bp_1[0] + ':' + str(self.bp_1[1]) + '-' + self.bp_2[0] + ':' + str(self.bp_2[1])
        fus = ''
        if self.genes[2] in ['possible_fusion', 'oncogenic_fusion']:
            fus = ',' + self.genes[0][0]+ '::'+ self.genes[1][0]
        return self.impact +':' +'[' + self.vcf_id + ',' + pos + fus + ','+ self.genes[2] +']'
            
class CNA(object):
    __slots__ = ("ref_id", "pos_1", "pos_2", "cn", "haplotype",'LOH', 'sv1', 'sv2', 'genes', 'dir1')
    def __init__(self, ref_id, pos_1, pos_2, cn, haplotype, LOH, sv1, sv2):
        self.ref_id = ref_id
        self.pos_1 = pos_1
        self.pos_2 = pos_2
        self.cn = cn
        self.haplotype = haplotype
        self.LOH = LOH
        self.sv1 = sv1
        self.sv2 = sv2
        self.genes = []
        self.dir1 = 'NEUT'
    def get_name(self):
        return 'HP' + str(self.haplotype) +':' + self.ref_id + ':' + str(self.pos_1) + '-' + str(self.pos_2)

class Gene(object):
    __slots__ = ("ref_id", "pos_1", "pos_2", "Gene_symbol", "CN", 'SV', 'score', 'cancer', 'impact', 'CN_impact')
    def __init__(self, ref_id, pos_1, pos_2, Gene_symbol):
        self.ref_id = ref_id
        self.pos_1 = pos_1
        self.pos_2 = pos_2
        self.Gene_symbol = Gene_symbol
        self.CN = ['NA', 'NA']   # per-haplotype copy number at the gene start; 'NA' = unknown
        self.SV = []
        self.score = [0, 0, 0]
        self.cancer = ''
        self.impact = ['','','']
        self.CN_impact = ['','']
    def to_str(self):
        pos = self.ref_id + ':' + str(self.pos_1) + '-' + str(self.pos_2)
        sv_inf = [[],[],[]]
        sv_info = []
        for sv, bp in self.SV:
            hp = sv.hp1 if bp == 1 else sv.hp2
            sv_inf[hp].append(sv.to_sum())
        for svs in sv_inf:
            if svs:
                svs= list(set(svs))
                svs.sort()
                sv_info.append(','.join(svs))
            else:
                sv_info.append('')
        return '\t'.join([self.Gene_symbol, pos, str(sum(self.score)), self.cancer, self.CN_impact[0], fmt_cn(self.CN[0]), self.impact[1], sv_info[1],  self.CN_impact[1], fmt_cn(self.CN[1]), self.impact[2], sv_info[2], self.impact[0], sv_info[0]])

def get_bps(vcf_file):
    svs = defaultdict(list)
    vcf = pysam.VariantFile(vcf_file)
    for var in vcf:
        if "_2" in var.id:
            continue
        bp_1 = (var.chrom, var.pos)
        dir_ls= var.info['STRANDS'] if 'STRANDS' in var.info.keys() else ('','')
        dir1, dir2 = (dir_ls[0], dir_ls[-1])
        
        ins_seq = '' if not var.info['SVTYPE'] == 'INS'  else var.alts[0]
        sample_id = var.samples.keys()[0]
        supp = var.samples[sample_id]['DV']
        vaf = var.samples[sample_id]['VAF']
        vntr = True if 'INSIDE_VNTR' in var.info.keys() else False
        has_ins = '' if not 'INSSEQ' in var.info.keys() else var.info['INSSEQ']
        is_single = False if not var.info['SVTYPE'] == 'sBND' else True
        sv_type = var.info['SVTYPE']
        vcf_id = var.id
        hp1 = hp2 = 0
        bp_2 = bp_1
        if 'HP' in var.info.keys():
            hp1, hp2 = [int(h) for h in var.info['HP'].split('|')]
        
        if var.info['SVTYPE'] == 'sBND':
            chr2 = var.chrom
            
        elif var.info['SVTYPE'] == 'BND':
            if '_2' in var.id:
                if 'HP' in var.info.keys():
                    hp2 = var.info['HP']
                    svs[var.id.replace('_2', '_1')].hp2 = hp2
                continue
            dir1, dir2 = var.info['STRANDS']
            chr2, pos2 = var.alts[0].replace('[','').replace(']','').replace('N', '').split(':')
            pos2 = int(pos2)
            bp_2 = (chr2, pos2)
            
        elif var.info['SVTYPE'] in ['DEL', 'DUP', 'INV']:
            end_pos = var.stop if var.stop else var.pos + var.info['SVLEN']
            bp_2 = (var.chrom, end_pos)
            
        elif var.info['SVTYPE'] == 'INS':
            bp_2 = bp_1
            
        cluster_id = ''
        if 'CLUSTERID' in var.info.keys():
            cluster_id = var.info['CLUSTERID']
            
        detailed_type = ''
        if 'DETAILED_TYPE' in var.info.keys():
            detailed_type = var.info['DETAILED_TYPE']
            
        svs[vcf_id]=(SV(bp_1, dir1, bp_2, dir2, supp, has_ins, sv_type,vaf, vcf_id, is_single, vntr,ins_seq, cluster_id, detailed_type, hp1, hp2))
    return list(svs.values())

def get_savana_bps(vcf_file):
    svs = defaultdict(list)
    vcf = pysam.VariantFile(vcf_file)
    seen_mates = set()
    for var in vcf:
        if var.id in seen_mates:
            continue
        bp_1 = (var.chrom, var.pos)
        sv_type = var.info.get('SVTYPE', '')
        dir_ls = var.info.get('BP_NOTATION', '')
        dir1, dir2 = (dir_ls[0], dir_ls[1]) if len(dir_ls) == 2 else ('','')

        supp = int(var.info.get('TUMOUR_READ_SUPPORT', 0))
        vaf_values = var.info.get('TUMOUR_AF', (0,))
        if not isinstance(vaf_values, (tuple, list)):
            vaf_values = (vaf_values,)
        vaf = max([float(x) for x in vaf_values]) if vaf_values else 0
        hp_counts = var.info.get('TUMOUR_ALT_HP', (0, 0, 0))
        hp1 = hp2 = savana_haplotype(hp_counts)
        vcf_id = re.sub(r'_[12]$', '', var.id)
        bp_2 = bp_1
        ins_seq = ''

        if sv_type == 'BND':
            mate = parse_bnd_alt(var.alts[0])
            if not mate:
                logger.warning("Skipping Savana BND with unparseable ALT: %s %s", var.id, var.alts[0])
                continue
            chr2, pos2 = mate
            bp_2 = (chr2, pos2)
            if 'MATEID' in var.info.keys():
                seen_mates.add(var.info['MATEID'])

        elif sv_type in ['DEL', 'DUP', 'INV']:
            svlen = abs(int(var.info.get('SVLEN', 0)))
            end_pos = var.stop if var.stop else var.pos + svlen
            bp_2 = (var.chrom, end_pos)

        elif sv_type == 'INS':
            bp_2 = bp_1
            if var.alts and not var.alts[0].startswith('<'):
                ins_seq = var.alts[0]

        cluster_id = var.info.get('MATEID', '')
        detailed_type = var.info.get('CLASS', '')
        svs[vcf_id]=(SV(bp_1, dir1, bp_2, dir2, supp, '', sv_type, vaf, vcf_id, False, False, ins_seq, cluster_id, detailed_type, hp1, hp2))
    return list(svs.values())

def savana_haplotype(hp_counts):
    if not isinstance(hp_counts, (tuple, list)) or len(hp_counts) < 2:
        return 0
    hp1, hp2 = int(hp_counts[0]), int(hp_counts[1])
    if hp1 == hp2:
        return 0
    return 1 if hp1 > hp2 else 2

def parse_bnd_alt(alt):
    match = BND_MATE_RE.search(alt)
    if not match:
        return None
    return match.group(1), int(match.group(2))

def get_SVs(vcf_file, caller='severus'):
    if caller == 'savana':
        return get_savana_bps(vcf_file)
    return get_bps(vcf_file)

def fill_missing_segments_cn1(
    chr_lens: Dict[str, int],
    hp1ls: Dict[str, List[List[int]]],
    default_cn: float = 1.0,
    assume_1based_inclusive: bool = True,
):
    hp1ls_filled: Dict[str, List[List[int]]] = {}
    for chrom , L in chr_lens.items():
        L = int(L)
        if assume_1based_inclusive:
            span_start, span_end = 1, L
        else:
            span_start, span_end = 0, L 
        if chrom not in hp1ls or len(hp1ls[chrom]) != 3 or len(hp1ls[chrom][0]) == 0:
            if assume_1based_inclusive:
                hp1ls_filled[chrom] = [[default_cn], [span_start], [span_end]]
            else:
                hp1ls_filled[chrom] = [[default_cn], [span_start], [span_end]]
            continue
        cn_list, s_list, e_list = hp1ls[chrom]
        segs = sorted(
            [(float(cn_list[i]), int(s_list[i]), int(e_list[i])) for i in range(len(cn_list))],
            key=lambda x: (x[1], x[2]),
        )
        clamped = []
        for cn, s, e in segs:
            if assume_1based_inclusive:
                s = max(s, span_start)
                e = min(e, span_end)
                if s <= e:
                    clamped.append((cn, s, e))
            else:
                s = max(s, span_start)
                e = min(e, span_end)
                if s < e:
                    clamped.append((cn, s, e))
        merged = []
        for cn, s, e in clamped:
            if not merged:
                merged.append([cn, s, e])
                continue
            pcn, ps, pe = merged[-1]
            if assume_1based_inclusive:
                overlap_or_touch = s <= pe + 1
            else:
                overlap_or_touch = s <= pe  # [ps,pe) touches when s==pe
            if overlap_or_touch:
                merged[-1][2] = max(pe, e)
            else:
                merged.append([cn, s, e])
        out_cn, out_s, out_e = [], [], []
        cur = span_start
        for cn, s, e in merged:
            if assume_1based_inclusive:
                if cur < s:
                    out_cn.append(default_cn)
                    out_s.append(cur)
                    out_e.append(s - 1)
                out_cn.append(cn); out_s.append(s); out_e.append(e)
                cur = e + 1
            else:
                if cur < s:
                    out_cn.append(default_cn)
                    out_s.append(cur)
                    out_e.append(s)
                out_cn.append(cn); out_s.append(s); out_e.append(e)
                cur = e
        if assume_1based_inclusive:
            if cur <= span_end:
                out_cn.append(default_cn)
                out_s.append(cur)
                out_e.append(span_end)
        else:
            if cur < span_end:
                out_cn.append(default_cn)
                out_s.append(cur)
                out_e.append(span_end)
        hp1ls_filled[chrom] = [out_cn, out_s, out_e]
    return hp1ls_filled

def pad_unlisted_segments(chr_lens, hpls, default_cn=1.0):
    """Return a copy of ``hpls`` ({chrom: [cn_list, start_list, end_list]}, 1-based inclusive) in which every
    chromosome of ``chr_lens`` is covered end to end: bases before the first listed segment, between non-adjacent
    segments and after the last one get a ``default_cn`` segment, and chromosomes without any listed segment get a
    single ``default_cn`` segment. Listed segments are kept exactly as they are (never merged, even when adjacent).
    Contigs absent from ``chr_lens`` are passed through unchanged."""
    out = {}
    for chrom, length in chr_lens.items():
        length = int(length)
        if chrom not in hpls or len(hpls[chrom]) != 3 or not hpls[chrom][0]:
            out[chrom] = [[default_cn], [1], [length]]
            continue
        cns, starts, ends = hpls[chrom]
        order = sorted(range(len(cns)), key=lambda i: (int(starts[i]), int(ends[i])))
        new_cn, new_s, new_e = [], [], []
        cur = 1
        for i in order:
            s, e = int(starts[i]), int(ends[i])
            if s > cur:
                new_cn.append(default_cn)
                new_s.append(cur)
                new_e.append(s - 1)
            new_cn.append(cns[i])
            new_s.append(s)
            new_e.append(e)
            cur = max(cur, e + 1)
        if cur <= length:
            new_cn.append(default_cn)
            new_s.append(cur)
            new_e.append(length)
        out[chrom] = [new_cn, new_s, new_e]
    for chrom, value in hpls.items():
        if chrom not in out:
            out[chrom] = value
    return out

def fmt_cn(value):
    """Copy number for the tables: integral values without a trailing .0, unknown as NA."""
    if value is None or value == 'NA' or (isinstance(value, float) and np.isnan(value)):
        return 'NA'
    try:
        f = float(value)
    except (TypeError, ValueError):
        return str(value)
    return str(int(f)) if f.is_integer() else str(value)

def is_unknown(cn):
    return cn is None or cn == 'NA' or (isinstance(cn, float) and np.isnan(cn))

def round_half_up(x):
    """Nearest integer, halves rounded up (Python's round() is banker's rounding: round(2.5) == 2)."""
    return int(np.floor(float(x) + 0.5))

def haplotype_baseline(ploidy):
    """Expected copy number of one haplotype (Wakhan HP1/HP2) or one allele (SAVANA major/minor) in a copy-neutral
    region of a tumour with the given ploidy: round(ploidy / 2)."""
    return max(0, round_half_up(float(ploidy) / 2.0))

def cn_state(cn, baseline):
    """'AMP' / 'DEL' / 'NEUT' for a copy number against the per-haplotype baseline; '' when the copy number is unknown."""
    if is_unknown(cn):
        return ''
    if cn > baseline:
        return 'AMP'
    if cn < baseline:
        return 'DEL'
    return 'NEUT'

def read_ploidy_file(path):
    """Tumour ploidy from a caller's fit table: SAVANA ``*_fitted_purity_ploidy.tsv`` (purity, ploidy, distance, rank) or
    Wakhan ``solutions_ranks.tsv`` (..., ploidy, confidence, solution_rank). Any TSV with a ``ploidy`` column works; the row
    with the lowest ``rank`` / ``solution_rank`` is used when such a column exists, otherwise the first row."""
    df = pd.read_csv(path, sep='\t')
    df.columns = [str(c).strip() for c in df.columns]
    if 'ploidy' not in df.columns:
        raise ValueError(f"{path}: no 'ploidy' column (columns: {', '.join(df.columns)})")
    df = df.dropna(subset=['ploidy'])
    if df.empty:
        raise ValueError(f"{path}: no row with a ploidy value")
    rank_col = next((c for c in ('rank', 'solution_rank') if c in df.columns), None)
    row = df.sort_values(rank_col).iloc[0] if rank_col else df.iloc[0]
    ploidy = float(row['ploidy'])
    if not np.isfinite(ploidy) or ploidy <= 0:
        raise ValueError(f"{path}: ploidy is not a positive number: {row['ploidy']!r}")
    return ploidy

def mean_cn(segments):
    """Length-weighted mean copy number over (cn, start, end) segments, skipping unknown copy numbers."""
    num = den = 0.0
    for cn, start, end in segments:
        if is_unknown(cn):
            continue
        length = max(1, int(end) - int(start) + 1)
        num += float(cn) * length
        den += length
    return num / den if den else float('nan')

def resolve_ploidy(ploidy, estimate, caller):
    """Return (ploidy, per-haplotype baseline); fall back to the profile-derived estimate with a warning."""
    if ploidy is None:
        ploidy = estimate if np.isfinite(estimate) else 2.0
        logger.warning("No tumour ploidy supplied (--ploidy / --ploidy-file): estimated %.2f from the %s copy-number "
                       "profile as the length-weighted mean total copy number. Prefer the caller's fitted ploidy.", ploidy, caller)
    base = haplotype_baseline(ploidy)
    logger.info("Tumour ploidy %.2f -> per-haplotype baseline %d (copy number > %d = AMP, < %d = DEL, = %d = NEUT)",
                ploidy, base, base, base, base)
    return ploidy, base

def round_savana_cn(total_cn, minor_cn):
    """SAVANA's absolute copy numbers are purity-scaled segment means and therefore fractional. Round them to integers so
    the AMP/DEL/NEUT call works like Wakhan's integer calls: total first, minor capped at floor(total/2) so that
    minor <= major, major = total - minor. Returns (major, minor); (None, None) when the total or the minor-allele copy
    number is missing (segments without het SNPs): the split is unknown and must not be mistaken for a loss."""
    if is_unknown(total_cn) or is_unknown(minor_cn):
        return None, None
    total_int = max(0, round_half_up(total_cn))
    minor_int = min(max(0, round_half_up(minor_cn)), total_int // 2)
    return total_int - minor_int, minor_int

def get_CNA(cna_vcf, svs, ploidy=None):
    vcf = pysam.VariantFile(cna_vcf)
    hp1ls = defaultdict(list)
    hp2ls = defaultdict(list)
    LOH  = defaultdict(list)
    THR=1
    CNAs = defaultdict(list)
    cov1 = []
    chrs = vcf.header.contigs.keys()
    chr_lens = [x.length for x in vcf.header.contigs.values() ]
    chr_lens = dict(zip(chrs, chr_lens))
    # Segment building below assumes position-sorted records; sort explicitly rather than rely on the writer.
    for var in sorted(vcf, key=lambda v: (v.chrom, v.pos, v.stop)):
        ref_id, pos_1, pos_2 = var.chrom, var.pos, var.stop
        hp1, hp2 = var.samples['Sample']['CN1'], var.samples['Sample']['CN2']
        if hp1:
            cov1.append( var.samples['Sample']['COV1']/hp1)
        if ref_id in hp1ls:
            if not hp1ls[ref_id][0][-1] == hp1:
                if abs(hp1ls[ref_id][2][-1] - pos_1) > THR:
                    hp1ls[ref_id][0].append(1.0)
                    hp1ls[ref_id][1].append(hp1ls[ref_id][2][-1]+1)
                    hp1ls[ref_id][2].append(pos_1-1)
                pos_1 = pos_1 if pos_1 > hp1ls[ref_id][2][-1] else pos_1 + 1
                hp1ls[ref_id][0].append(hp1)
                hp1ls[ref_id][1].append(pos_1)
                hp1ls[ref_id][2].append(pos_2)
            elif pos_1 <= hp1ls[ref_id][2][-1] + THR + 1:
                # same copy number and adjacent (or overlapping): extend the open segment
                hp1ls[ref_id][2][-1] = max(pos_2, hp1ls[ref_id][2][-1])
            else:
                # same copy number but separated by unlisted (neutral) bases: a new segment, the gap is padded later
                hp1ls[ref_id][0].append(hp1)
                hp1ls[ref_id][1].append(pos_1)
                hp1ls[ref_id][2].append(pos_2)
            if not hp2ls[ref_id][0][-1] == hp2:
                if abs(hp2ls[ref_id][2][-1] - pos_1) > THR:
                    hp2ls[ref_id][0].append(1.0)
                    hp2ls[ref_id][1].append(hp2ls[ref_id][2][-1]+1)
                    hp2ls[ref_id][2].append(pos_1-1)
                pos_1 = pos_1 if pos_1 > hp2ls[ref_id][2][-1] else pos_1 + 1
                hp2ls[ref_id][0].append(hp2)
                hp2ls[ref_id][1].append(pos_1)
                hp2ls[ref_id][2].append(pos_2)
            elif pos_1 <= hp2ls[ref_id][2][-1] + THR + 1:
                # same copy number and adjacent (or overlapping): extend the open segment
                hp2ls[ref_id][2][-1] = max(pos_2, hp2ls[ref_id][2][-1])
            else:
                # same copy number but separated by unlisted (neutral) bases: a new segment, the gap is padded later
                hp2ls[ref_id][0].append(hp2)
                hp2ls[ref_id][1].append(pos_1)
                hp2ls[ref_id][2].append(pos_2)
        else:
            hp1ls[ref_id]= [[hp1],[pos_1],[pos_2]]
            hp2ls[ref_id] = [[hp2],[pos_1], [pos_2]]
        if hp1 == 0 or hp2 == 0:
            if ref_id in LOH.keys():
                LOH[ref_id][0].append(pos_1)
                LOH[ref_id][1].append(pos_2)
            else:
                LOH[ref_id]= [[pos_1],[pos_2]]

    cn1_cov = int(np.median(cov1)) *0.75
    # Wakhan's integer VCF lists altered segments only. Pad every primary chromosome (and any contig that has a
    # listed segment) with neutral one-copy-per-haplotype segments. Without this, annot_CNAs() looks genes up by
    # bisect in a list that does not cover them: a gene before the first listed segment lands on index -1 (the
    # chromosome's LAST segment), a gene after the last one on that last segment, and chromosomes without any
    # listed segment are skipped altogether.
    primary = re.compile(r'^(chr)?([0-9]{1,2}|X|Y)$')
    pad_lens = {c: L for c, L in chr_lens.items() if L and (c in hp1ls or primary.match(c))}
    hp1ls = pad_unlisted_segments(pad_lens, hp1ls)
    hp2ls = pad_unlisted_segments(pad_lens, hp2ls)
    for hp, hpls in enumerate([hp1ls, hp2ls]):
        for ref_id, (cnls, pos1ls, pos2ls) in hpls.items():
            loh = False
            for i, (pos_1, pos_2) in enumerate(zip(pos1ls, pos2ls)):
                cn = cnls[i]
                ls1 = [(ref_id, pos_1),(ref_id, pos_1+1),(ref_id, pos_1-1)]
                sv1 = [sv for sv in svs if sv.bp_1 in ls1 and sv.direction_1 == '-' or sv.bp_2 in ls1 and sv.direction_2 == '-']
                sv1 = sv1[0] if sv1 else ''
                ls2 = [(ref_id, pos_2),(ref_id, pos_2+1),(ref_id, pos_2-1)]
                sv2 = [sv for sv in svs if sv.bp_1 in ls2 and sv.direction_1 == '+' or sv.bp_2 in ls2 and sv.direction_2 == '+' ]
                sv2 = sv2[0] if sv2 else ''
                if ref_id in LOH.keys() and  pos_1 in LOH[ref_id][0]:
                    loh = True
                CNAs[(ref_id, hp+1)].append(CNA(ref_id, pos_1, pos_2, cn, hp+1, loh, sv1, sv2))

    # Baseline = the tumour ploidy (fitted by the caller, passed in by the user), not a mean of this profile. Both padded
    # haplotype profiles cover the same bases, so the fallback estimate of the ploidy is mean(HP1) + mean(HP2).
    estimate = mean_cn([(c, s, e) for cnls, s1, e1 in hp1ls.values() for c, s, e in zip(cnls, s1, e1)]) + \
               mean_cn([(c, s, e) for cnls, s1, e1 in hp2ls.values() for c, s, e in zip(cnls, s1, e1)])
    ploidy, base = resolve_ploidy(ploidy, estimate, 'Wakhan')
    check_hp_svs(CNAs)
    for cnas in CNAs.values():
        for cna in cnas:
            cna.dir1 = cn_state(cna.cn, base)
    check_cn_altering_svs(svs, cn1_cov)
    return (CNAs, [base, base])

def get_savana_CNA(cna_tsv, svs, ploidy=None):
    """SAVANA ``*_segmented_absolute_copy_number.tsv``: HP1 = major allele (copyNumber - minorAlleleCopyNumber), HP2 =
    minor allele, both rounded to integers (see round_savana_cn). SAVANA has no haplotype phasing, so 'HP1'/'HP2' are
    the larger and the smaller allele of each segment; both are compared with the same baseline, round(ploidy / 2)."""
    df = pd.read_csv(cna_tsv, sep='\t')
    hp1ls = defaultdict(list)
    hp2ls = defaultdict(list)
    LOH = defaultdict(list)
    CNAs = defaultdict(list)
    totals = []
    n_unknown = 0
    for _, row in df.sort_values(['chromosome', 'start', 'end']).iterrows():
        ref_id = row['chromosome']
        pos_1, pos_2 = int(row['start']), int(row['end'])
        total_cn = pd.to_numeric(row['copyNumber'], errors='coerce')
        minor_cn = pd.to_numeric(row['minorAlleleCopyNumber'], errors='coerce')
        major_int, minor_int = round_savana_cn(total_cn, minor_cn)
        if minor_int is None:
            n_unknown += 1
        hp1ls[ref_id].append((major_int, pos_1, pos_2))
        hp2ls[ref_id].append((minor_int, pos_1, pos_2))
        totals.append((total_cn, pos_1, pos_2))
        if minor_int == 0:
            LOH[ref_id].append((pos_1, pos_2))
    if n_unknown:
        logger.warning("%d SAVANA segment(s) have no minor-allele copy number (no het SNPs): their allele copy numbers "
                       "are left unknown (NA), not treated as a loss", n_unknown)

    ploidy, base = resolve_ploidy(ploidy, mean_cn(totals), 'SAVANA')
    for hp, hpls in enumerate([hp1ls, hp2ls]):
        for ref_id, segments in hpls.items():
            for cn, pos_1, pos_2 in segments:
                sv1 = find_segment_boundary_sv(svs, ref_id, pos_1, '-')
                sv2 = find_segment_boundary_sv(svs, ref_id, pos_2, '+')
                loh = any(loh_start == pos_1 and loh_end == pos_2 for loh_start, loh_end in LOH.get(ref_id, []))
                CNAs[(ref_id, hp+1)].append(CNA(ref_id, pos_1, pos_2, cn, hp+1, loh, sv1, sv2))

    check_hp_svs(CNAs)
    for cnas in CNAs.values():
        for cna in cnas:
            cna.dir1 = cn_state(cna.cn, base)
    return (CNAs, [base, base])

def find_segment_boundary_sv(svs, ref_id, pos, direction):
    boundary = [(ref_id, pos),(ref_id, pos+1),(ref_id, pos-1)]
    if direction == '-':
        candidates = [sv for sv in svs if (sv.bp_1 in boundary and sv.direction_1 == '-') or (sv.bp_2 in boundary and sv.direction_2 == '-')]
    else:
        candidates = [sv for sv in svs if (sv.bp_1 in boundary and sv.direction_1 == '+') or (sv.bp_2 in boundary and sv.direction_2 == '+')]
    return candidates[0] if candidates else ''

def get_CNAs(cna_file, svs, caller='wakhan', ploidy=None):
    """Returns (CNAs, baselines): baselines = [per-haplotype baseline] * 2, derived from the tumour ploidy."""
    if caller == 'savana':
        return get_savana_CNA(cna_file, svs, ploidy)
    return get_CNA(cna_file, svs, ploidy)

def check_cn_altering_svs(svs, cn1_cov):
    for sv in svs:
        sv.cn_altering = round(sv.supp/cn1_cov)

def check_hp_svs(cnas):
    multsvs = defaultdict(list)
    for cns in cnas.values():
        for cn in cns:
            if cn.sv1:
                multsvs[cn.sv1].append(cn)
            if cn.sv2:
                multsvs[cn.sv2].append(cn)
    for sv, cns in multsvs.items():
        xx = Counter([(cn.ref_id, cn.pos_1) for cn in cns])
        xx = [x for x in xx.values() if x > 1]
        
        yy = Counter([(cn.ref_id, cn.pos_2) for cn in cns])
        yy = [x for x in yy.values() if x > 1]
        
        if not xx and not yy:
            continue
        hp2 = [cn for cn in cns if cn.haplotype == 2]
        new_sv = copy.copy(sv)
        for cn in hp2:
            if cn.sv1 == sv:
                cn.sv1 = new_sv
            else:
                cn.sv2 = new_sv
    
GENE_PAD = 10000   # bp added on both sides of every gene in get_genes(); hits inside the pad are 'promoter/utr'

def get_genes(gff_file):
    THR = GENE_PAD
    genes = defaultdict(list)
    exon_pos = defaultdict(list)
    fopen = open_text(gff_file)
    for line in fopen:
        ref_id, gene_name, strand, typ, start, end = line.split()
        start = int(start)
        end = int(end)
        THR1, THR2 = (THR, 0) if strand == '+' else (0, THR)
        if typ == 'gene':
            if ref_id in genes.keys():
                genes[ref_id][0].append(gene_name)
                genes[ref_id][1].append(start-THR)
                genes[ref_id][2].append(end+THR)
            else:
                genes[ref_id] = [[],[],[]]
                genes[ref_id][0].append(gene_name)
                genes[ref_id][1].append(start-THR)
                genes[ref_id][2].append(end+THR)
        elif not re.fullmatch(r'exon\d+', typ):
            # transcript rows (<gene>-NNN, or ENST... for unnamed ENSG genes), codons and UTRs are not exons
            continue
        else:
            if gene_name in exon_pos.keys():
                exon_pos[gene_name][0].append(start)
                exon_pos[gene_name][1].append(end)
                exon_pos[gene_name][2].append(typ)
            else:
                exon_pos[gene_name] = [[],[],[],[strand]]
                exon_pos[gene_name][0].append(start)
                exon_pos[gene_name][1].append(end)
                exon_pos[gene_name][2].append(typ)
    for gene_name,exons in exon_pos.items():
        exonls = [ind for ind, typ in enumerate(exons[2]) if 'exon' in typ]
        utrs5 = [ind for ind, typ in enumerate(exons[2]) if typ == 'five_prime_UTR']
        utrs3 = [ind for ind, typ in enumerate(exons[2]) if typ == 'three_prime_UTR']
        if utrs5:
            st5 = min([exons[0][i] for i in utrs5])
            end5 = max([exons[1][i] for i in utrs5])
        if utrs3:
            st3 = min([exons[0][i] for i in utrs3])
            end3 = max([exons[1][i] for i in utrs3])
        if utrs5 and utrs3:
            exons[0] = [st5] + [exons[0][i] for i in exonls] + [st3]
            exons[1] = [end5] + [exons[1][i] for i in exonls] + [end3]
            exons[2] = ['five_prime_UTR'] + [exons[2][i] for i in exonls] + ['three_prime_UTR']
        elif utrs5:
            exons[0] = [st5] + [exons[0][i] for i in exonls]
            exons[1] = [end5] + [exons[1][i] for i in exonls]
            exons[2] = ['five_prime_UTR'] + [exons[2][i] for i in exonls]
        elif utrs3:
            exons[0] = [exons[0][i] for i in exonls] + [st3]
            exons[1] = [exons[1][i] for i in exonls] + [end3]
            exons[2] = [exons[2][i] for i in exonls] + ['three_prime_UTR']
        exons[2] = [x for _, x in sorted(zip(exons[0], exons[2]))]
        exons[1] = [x for _, x in sorted(zip(exons[0], exons[1]))]
        exons[0] = sorted(exons[0])
    for ref_id, genels in genes.items():
        genes[ref_id] = sort_gene_lists(genels)
    return (genes, exon_pos)

def sort_gene_lists(genels):
    """[names, padded starts, padded ends] -> the same lists ordered by padded start (ties by end, then name), plus a fourth
    list with the running maximum of the padded ends. The lookup in annotBPs() bisects the start list, so it must be
    sorted; the file order is not (the strand-dependent promoter padding reorders neighbours, and the bundled files used to
    carry each chromosome's last gene at the front of the next chromosome)."""
    names, starts, ends = genels[0], genels[1], genels[2]
    order = sorted(range(len(names)), key=lambda i: (starts[i], ends[i], names[i]))
    names = [names[i] for i in order]
    starts = [starts[i] for i in order]
    ends = [ends[i] for i in order]
    max_end, running = [], float('-inf')
    for e in ends:
        running = max(running, e)
        max_end.append(running)
    return [names, starts, ends, max_end]

def find_gene_at(genels, pos):
    """Index of the gene to annotate a breakpoint at ``pos`` with, or None. Candidates are the genes whose padded span
    (gene +/- GENE_PAD) contains ``pos``; a gene whose body (the unpadded gene) contains it beats genes that only reach it
    with their flank, and among equals the one with the largest start (the innermost of nested genes) wins. The walk back
    from the bisect position stops as soon as no gene further left can reach ``pos`` (running maximum of the ends)."""
    names, starts, ends = genels[0], genels[1], genels[2]
    max_end = genels[3] if len(genels) > 3 else None
    idx = bisect.bisect_right(starts, pos) - 1
    best_pad = None
    while idx >= 0:
        if max_end is not None and max_end[idx] < pos:
            break
        if ends[idx] >= pos:
            if starts[idx] + GENE_PAD <= pos <= ends[idx] - GENE_PAD:
                return idx                      # inside the gene body: the innermost such gene
            if best_pad is None:
                best_pad = idx                  # flank hit; keep looking for a gene body further left
        idx -= 1
    return best_pad

def exon_label(exons, pos, strand):
    """Where ``pos`` falls relative to a gene's exons (exons = [starts, ends, names, [strand]], sorted by start):
    'exonK' inside exon K; 'intronK' between exon K and exon K+1 in transcript order (the lower number of the two
    flanking exons, so the label is right on both strands and for K >= 10); in the padded flank before the first
    exon of the transcript 'promoter/utr', after its last exon 'downstream'. Distinct labels for the two flanks
    keep an SV spanning the whole gene a 'between_exons' event in annot_SVS()."""
    starts, ends, names = exons[0], exons[1], exons[2]
    n = len(starts)
    after = bisect.bisect_right(starts, pos)     # exons starting at or before pos
    ended = bisect.bisect_left(ends, pos)        # exons ending before pos
    if after == 0 or ended == n:
        five_prime = (after == 0) == (strand != '-')
        return 'promoter/utr' if five_prime else 'downstream'
    if after == ended:
        nums = [int(m.group(1)) for m in (re.fullmatch(r'exon(\d+)', x) for x in (names[after - 1], names[after])) if m]
        return 'intron' + str(min(nums)) if nums else 'intron'
    return names[after - 1]

def annotBPs(sv, bp_1, genes, exon_pos, by_gene, bp):
    if not bp_1[0] in list(genes.keys()) and not 'chr' + bp_1[0] in list(genes.keys()):
        sv.genes.append(())
    else:
        if bp_1[0] in list(genes.keys()):
            genels = genes[bp_1[0]]
        else:
            genels = genes['chr' + bp_1[0]]
        # The old test `bisect_right(starts) == bisect_left(ends) + 1` assumed both padded lists are sorted; the end list
        # is not wherever genes overlap, so the answer depended on the search path (and on the file order of the genes).
        ind2 = find_gene_at(genels, bp_1[1])
        if ind2 is not None:
            ty, strand = 'promoter/utr', '.'
            if genels[0][ind2] in exon_pos:
                exons = exon_pos[genels[0][ind2]]
                strand = exons[3][0]
                ty = exon_label(exons, bp_1[1], strand)
            sv.genes.append((genels[0][ind2], ty, strand))
            add_gene(by_gene, genels[0][ind2], bp_1[0], genels[1][ind2], genels[2][ind2])
            by_gene[genels[0][ind2]].SV.append((sv, bp))
        else:
            sv.genes.append(())
            
def annot_SVS(genes, exon_pos, svs, by_gene):
    for sv in svs:
        annotBPs(sv, sv.bp_1, genes, exon_pos, by_gene,1)
        annotBPs(sv, sv.bp_2, genes, exon_pos, by_gene,2)
    for sv in svs:
        if not sv.genes[1] and not sv.genes[0]:
            sv.genes.append('')
            continue
        if not sv.genes[1] or not sv.genes[0]:
            sv.genes.append('Gene-NonCoding')
            sv.score += 2
            sv.impact = 'LOW'
            continue
        if sv.genes[0] == sv.genes[1]:
            if 'exon'  in sv.genes[0][1]:
                svlen = sv.bp_2[1] - sv.bp_1[1]
                if svlen % 3:
                    sv.genes.append('frameshift')
                    sv.score += 3
                    sv.impact = 'HIGH'
                else:
                    sv.genes.append('exon_altering')
                    sv.score += 3
                    sv.impact = 'HIGH'
            elif 'intron' in sv.genes[0][1]:
                sv.genes.append('intronic')
                sv.score += 1
                sv.impact = 'LOW'
            else:
                sv.genes.append('promoter/UTR')
                sv.score += 1
                sv.impact = 'LOW'
        else:
            if sv.genes[0][0] == sv.genes[1][0]:
                if not sv.genes[0][1] == sv.genes[1][1]:
                    sv.genes.append('between_exons')
                    sv.score += 3
                    sv.impact = 'HIGH'
            elif not sv.genes[0][0] == sv.genes[1][0]:
                if sv.direction_1 == sv.direction_2 and not sv.genes[0][2] == sv.genes[1][2]:
                    sv.genes.append('possible_fusion')
                    sv.score += 3
                    sv.impact = 'HIGH'
                elif not sv.direction_1 == sv.direction_2 and sv.genes[0][2] == sv.genes[1][2]:
                    sv.genes.append('possible_fusion')
                    sv.score += 3
                    sv.impact = 'HIGH'
                else:
                    sv.genes.append('between_genes')
                    sv.score += 2
                    sv.impact = 'LOW'
                    
def add_gene(by_gene, gene, ref_id, pos1, pos2):
    if not gene in by_gene.keys():
        by_gene[gene] = Gene(ref_id, pos1, pos2, gene)

def annot_CNAs(genes, cnas, baselines, by_gene):
    hps = [1,2]
    for ref_id, genels in genes.items():
        if not (ref_id,1) in cnas.keys():
            continue
        for hp in hps:
            cn_prof = cnas[ref_id, hp]
            cn_start = [cn.pos_1 for cn in cn_prof]
            cn_end = [cn.pos_2 for cn in cn_prof]
            cn = [cn.cn for cn in cn_prof]
            for i, gene in enumerate(genels[0]):
                ind1 = bisect.bisect_right(cn_start, genels[1][i])
                ind2 = bisect.bisect_left(cn_end, genels[2][i])
                add_gene(by_gene, gene, ref_id, genels[1][i], genels[2][i])
                if by_gene[gene].ref_id != ref_id:
                    # by_gene is keyed by symbol: a symbol annotated on two contigs (the 18 PAR genes on chrX and chrY)
                    # keeps the copy number of the contig seen first instead of being overwritten by the second.
                    continue
                if genels[1][i] > cn_end[-1]:
                    # gene start beyond the last segment (i.e. beyond the contig length): the annotation row is
                    # mislabelled (hg38.gff3.gz lists the last gene of each chromosome under the next chromosome);
                    # leave its copy number unset rather than assign the trailing segment's.
                    if hp == 1:
                        logger.warning("Gene %s at %s:%d-%d lies beyond the end of %s; copy number left unset",
                                       gene, ref_id, genels[1][i], genels[2][i], ref_id)
                    continue
                state = cn_state(cn[ind1-1], baselines[hp-1])
                if not ind1 - ind2 == 1:
                    cn_prof[ind1-1].genes.append((gene, 'disturbed'))
                elif state in ('AMP', 'DEL'):
                    cn_prof[ind1-1].genes.append(gene)
                    cn_prof[ind1-1].dir1 = state
                by_gene[gene].CN_impact[hp-1] = cn_prof[ind1-1].dir1
                by_gene[gene].CN[hp-1] = 'NA' if state == '' else cn[ind1-1]

def check_complexSV(cnas, svs):
    svls = defaultdict(list)
    for sv in svs:
         svls[sv.vcf_id] = sv

def run_command(cmd):
    result = subprocess.run(
        ["bash", "-o", "pipefail", "-c", cmd],
        capture_output=True,
        text=True,
    )
    if result.returncode != 0:
        logger.error("External command failed (%d): %s", result.returncode, cmd)
        if result.stderr:
            logger.error(result.stderr.rstrip())
        raise RuntimeError(f"External command failed: {cmd}")

def write_ins(svs, ref, t, specie, run_repeatmasker):
    fa_out = open('temp_ins.fa', 'w')
    wrote_sequence = False
    for sv in svs:
        if sv.ins_seq:
            fa_out.write('>' + sv.vcf_id + '\n')
            fa_out.write(sv.ins_seq)
            fa_out.write('\n')
            wrote_sequence = True
        elif sv.has_ins:
            fa_out.write('>' + sv.vcf_id + '\n')
            fa_out.write(sv.has_ins)
            fa_out.write('\n')
            wrote_sequence = True
    fa_out.close()
    if not wrote_sequence:
        return False
    if run_repeatmasker:
        run_command(f"RepeatMasker -species {specie} temp_ins.fa")
    run_command(f"minimap2 -ax map-ont {ref} temp_ins.fa -k 17 -y -K 5G -t {t} --eqx | samtools sort -@ {t} -m 4G > temp_ins.bam")
    run_command(f"samtools index -@ {t} temp_ins.bam")
    return True
    

def get_repeat(svls):
    if os.path.exists('temp_ins.fa.out'):
        rep_file = open('temp_ins.fa.out', 'r')
    else:
        logger.warning("RepeatMasker output file not found: temp_ins.fa.out")
        return ''
    # RepeatMasker does not group the .out lines by query (a query's hits can appear in several blocks), and the
    # old streaming parser summarised a query every time its id changed and never after the last line. Accumulate
    # per query over the whole file instead, then summarise each query once.
    masked = defaultdict(lambda: defaultdict(int))   # query id -> repeat family -> masked bases
    length = {}                                        # query id -> inserted-sequence length
    ll = 0
    for line in rep_file:
        ll += 1
        if ll < 4:
            continue
        l = line.split()
        if len(l) < 11:
            # blank line, or RepeatMasker's "There were no repetitive sequences detected" notice
            continue
        masked[l[4]][l[10]] += int(l[6]) - int(l[5])
        length[l[4]] = int(l[7][1:-1]) + int(l[6])
    rep_file.close()
    for sv_id, reps in masked.items():
        # Keep the dominant repeat family when it masks more than 80% of the inserted sequence.
        max_rep = max(reps, key=reps.get)
        replen = reps[max_rep]
        if replen > length[sv_id] * 0.8:
            if sv_id in svls:
                svls[sv_id].repeat.append((max_rep, replen))
            else:
                logger.warning("RepeatMasker hit for an insertion id not in the SV set: %s", sv_id)

def check_homology(seq, sv):
    seq1, seq2 = seq
    THR = 25
    MINSCORE = 3
    aligner = Align.PairwiseAligner()
    aligner.mode = 'local'
    aligner.gap_score = -100
    aligner.mismatch_score = -100
    ss = []
    for i in [THR-10,THR]:
        for j in [THR-10,THR]:
            a = aligner.align(seq1[i:i+11], seq2[j:j+11])
            if not a:
                continue
            else:
                aa = a[0]
                if aa.score < MINSCORE:
                    continue
                pos1 = abs(THR-i-aa.aligned[0][0][0]) + abs(THR-j - aa.aligned[1][0][0])
                if aa.score-pos1 > 0: 
                    ss.append((aa.score-pos1, aa.sequences[0][aa.aligned[0][0][0]:aa.aligned[0][0][1]]))
    if ss:
        ss.sort(key=lambda s:s[0])
        sv.microh.append(ss[-1][1])
    
def get_microhomology(svs, ref):
    THR = 25
    bpls = []
    f = open("temp_bpseq.tsv", "w")
    for sv in svs:
        if sv.sv_type == 'INS':
            continue
        bpls.append(sv)
        run_command(f'samtools faidx {ref} \'{sv.bp_1[0]}:{max(0, sv.bp_1[1]-THR)}-{sv.bp_1[1]+THR}\' >> temp_bpseq.tsv')
        run_command(f'samtools faidx {ref} \'{sv.bp_2[0]}:{max(0, sv.bp_2[1]-THR)}-{sv.bp_2[1]+THR}\' >> temp_bpseq.tsv')
    
    f.close()
  
    f = open('temp_bpseq.tsv', 'r')
    t = 0
    seq = []
    tt = 0
    for line in f:
        if line.startswith('>'):
            continue
        else:
            seq.append(line.strip())
            t+=1
        if t == 2:
            check_homology(seq, bpls[tt])
            t = 0
            seq = []
            tt+=1

def get_align(svls):
    MAPQ_THR = 45
    aln_file = pysam.AlignmentFile('temp_ins.bam', "rb")
    for aln in aln_file.fetch():
        if not aln.is_secondary and not aln.is_unmapped and aln.mapq > MAPQ_THR:
            pos = aln.reference_name +':' + str(aln.reference_start) +'-' + str(aln.reference_end)
            svls[aln.query_name].ins_align.append(pos)
    
def get_tel(svs):
    telseq = ['CCCTAA', 'TTAGGG']
    MIN_THR = 4
    for sv in svs:
        if sv.ins_seq or sv.has_ins:
            seq = sv.ins_seq if sv.ins_seq else sv.has_ins
            c1 = seq.count(telseq[0])
            if c1 < MIN_THR:
                c1 = seq.count(telseq[1])
            if c1 >= MIN_THR:
                sv.telomere = c1

def annot_bp_repeat(svls, rm_bed):
    f = open('temp_bps.bed', 'w')
    for sv in svls.values():
        f.write('\t' .join([sv.bp_1[0], str(max(sv.bp_1[1]-5,0)), str(sv.bp_1[1]+5), sv.vcf_id, 'BP1']))
        f.write('\n')
        if not sv.sv_type == 'INS':
            f.write('\t' .join([sv.bp_2[0], str(max(0, sv.bp_2[1]-5)), str(sv.bp_2[1]+5), sv.vcf_id, 'BP2']))
            f.write('\n')
            
    run_command(f'bedtools intersect -a temp_bps.bed -b {rm_bed} -wb > temp_int.bed')
    
    f1 = open('temp_int.bed')
    for line in f1:
        l = line.split()
        ind = 0 if l[4] == 'BP1' else 1
        svls[l[3]].repeat_bp[ind] = l[-1]
        
def annot_ins(svs, ref,t, rm_bed, specie, run_repeatmasker):
    svls = defaultdict(list)
    for sv in svs:
         svls[sv.vcf_id] = sv
    has_insertions = write_ins(svs, ref, t, specie, run_repeatmasker)
    if has_insertions:
        get_repeat(svls)
        get_align(svls)
    else:
        logger.info("No insertion sequences found; skipping insertion repeat/alignment annotation")
    get_tel(svs)
    annot_bp_repeat(svls, rm_bed)
    
    
def get_cancer_anno(sv, gene, cancer_genes):
    if gene in ['between_exons','exon_altering', 'frameshift']:
        sv.score += 3
        sv.impact = 'HIGH_ONCO'
    elif gene in ['between_genes', 'Gene-NonCoding']:
        sv.score += 2
        sv.impact = 'HIGH_ONCO'
    elif gene in ['intronic', 'promoter/UTR']:
        sv.score += 1
        sv.impact = 'MODERATE'
    

def cancer_annot_svs(svs, cancer_genes, fusion):
    for sv in svs:
        if not sv.genes[2]:
            continue
        if sv.genes[2] == 'possible_fusion':
            g1, g2 = sv.genes[0][0], sv.genes[1][0]
            if g2 in fusion.get(g1, ()) or g1 in fusion.get(g2, ()):
                sv.score += 3
                sv.impact = 'HIGH_ONCO'
                sv.genes[2] = 'oncogenic_fusion'
        else:
            if sv.genes[0] and sv.genes[0][0] in cancer_genes.keys():
                sv.cancer[0] = cancer_genes[sv.genes[0][0]]
                get_cancer_anno(sv, sv.genes[2], cancer_genes)
            if sv.genes[1] and sv.genes[1][0] in cancer_genes.keys():
                if sv.genes[0] and sv.genes[0][0] == sv.genes[1][0]:
                    sv.cancer[1] = sv.cancer[0]
                else:
                    sv.cancer[1] = cancer_genes[sv.genes[1][0]]
                    get_cancer_anno(sv, sv.genes[2], cancer_genes)

def cancer_annot_genes(by_gene, cancer_genes):
    for gene, gs in by_gene.items():
        if gene in cancer_genes.keys():
            gs.cancer = cancer_genes[gene]
        if 'TSG' in gs.cancer:
            for hp, cn in enumerate(gs.CN_impact):
                if cn == 'DEL':
                    gs.score[hp+1] += 3
                    gs.impact[hp+1] = 'TSG_DEL'
        elif 'oncogene' in gs.cancer:
            for hp, cn in enumerate(gs.CN_impact):
                if cn == 'AMP':
                    gs.score[hp+1] += 3
                    gs.impact[hp+1] = 'ONC_AMP' 
        for (sv, bp) in gs.SV:
            hp = sv.hp1 if bp == 1 else sv.hp2
            if sv.genes[2] in ['between_exons','exon_altering', 'frameshift']:
                gs.score[hp] += 3
                if not gs.impact[hp]:
                    gs.impact[hp] = sv.impact
            elif sv.genes[2] in ['between_genes', 'Gene-NonCoding']:
                gs.score[hp] += 2
                if not gs.impact[hp]:
                    gs.impact[hp] = sv.impact
            elif sv.genes[2] == 'oncogenic_fusion':
                gs.score[hp] += 2
                gs.impact[hp] = 'Oncogenic_fusion'
                             
BEDS_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'beds'))
# Bundled copy of https://raw.githubusercontent.com/aysegokce/PadfootFiles/refs/heads/main/cancer_genes.tsv
# (PadfootFiles commit 92c1842e, 2025-07-11; retrieved 2026-10-02; see beds/README.md). Padfoot used to download it on
# every human run, which made the tool fail on nodes without internet access.
CANCER_GENES_PATH = os.path.join(BEDS_DIR, 'cancer_genes.tsv')

def load_cancer_genes(path=None):
    """Read the cancer-gene table (columns Gene_symbol, Role_in_cancer, fusion, ...) into the two lookups
    cancer_annot_svs() / cancer_annot_genes() use: {symbol: role} and {symbol: [fusion partner, ...]}."""
    path = path or CANCER_GENES_PATH
    df = pd.read_csv(path, sep='\t')
    missing = [c for c in ('Gene_symbol', 'Role_in_cancer', 'fusion') if c not in df.columns]
    if missing:
        raise ValueError(f"{path}: missing column(s) {', '.join(missing)}")
    cancer_genes = defaultdict(list)
    fusion = defaultdict(list)
    for index, row in df.iterrows():
        if not pd.isna(row['fusion']):
            fusion[row['Gene_symbol']] = row['fusion'].split(',')
        if not pd.isna(row['Role_in_cancer']):
            cancer_genes[row['Gene_symbol']] = row['Role_in_cancer']
    return cancer_genes, fusion

def cancer_annot(svs, by_gene, path=None):
    logger.info("Cancer gene table: %s", path or CANCER_GENES_PATH)
    cancer_genes, fusion = load_cancer_genes(path)
    cancer_annot_svs(svs, cancer_genes, fusion)
    cancer_annot_genes(by_gene, cancer_genes)

def output_svs(svs, out_dir):
    fopen = open(out_dir + '/annotated_svs.tsv', 'w')
    header = '\t'.join(['SV ID', 'chr_1', 'pos_1', 'chr_2', 'pos_2', 'N_supp', 'Vaf', 'Strand','Score', 'Impact',\
                         'Role in Cancer1', 'Gene Name1',  'Pos_1', 'Strand_1', 'Repeat_Anno1', 'Role in Cancer2', 'Gene Name2',  'Pos_2', \
                             'Strand_2', 'Repeat_Anno2', 'Type', 'Microhomology', 'VNTR', \
                                 'Telomere_repeat', 'Aligned_pos(INS)', 'Repeat_annot'])
    fopen.write(header)
    fopen.write('\n')
    for sv in svs:
        fopen.write(sv.to_str())
        fopen.write('\n')
    fopen.close()
    
def output_genes(by_gene, out_dir):
    f = open(out_dir + '/by_gene.tsv', 'w')
    header = '\t'.join(['Gene Name','Gene Pos','Score','Role in Cancer', 'HP1', 'HP1 CN','HP1 Impact','HP1 SV', 'HP2', 'HP2 CN','HP2 Impact', 'HP2 SV', 'Unphased annot', 'Unphased SV'])
    f.write(header)
    f.write('\n')
    for gene, val in by_gene.items():
        f.write(val.to_str())
        f.write('\n')
    f.close()


def annotate_things(args):
    cna_vcf, vcf_file, t, ref, out_dir = args.cna_vcf, args.vcf_file, args.threads, args.ref, args.out_dir
    by_gene = defaultdict(list)
    (genes, exon_pos) = get_genes(args.gff_file)
    svs = get_SVs(vcf_file, args.sv_caller)
    cnas = []
    if cna_vcf:
        cnas, baselines = get_CNAs(cna_vcf, svs, args.cna_caller, getattr(args, 'ploidy', None))
        annot_CNAs(genes, cnas, baselines, by_gene)
    annot_SVS(genes, exon_pos, svs, by_gene)
    if args.specie == 'human':
        cancer_annot(svs, by_gene, getattr(args, 'cancer_genes', None))
    annot_ins(svs, ref,t, args.rm_file, args.specie, args.run_repeatmasker)
    get_microhomology(svs, ref)
    output_svs(svs, out_dir)
    output_genes(by_gene, out_dir)
    return (genes, cnas, exon_pos, svs, by_gene)
    

