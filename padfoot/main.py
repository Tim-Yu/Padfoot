#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun  6 14:05:48 2023

@author: keskusa2
"""

#!/usr/bin/env python3

import sys
import shutil
import argparse
import os
import importlib
import logging

from padfoot.annot import annotate_things, read_ploidy_file
from padfoot.preprocess import generate_gff, generate_rm
#from padfoot_cluster_cns import cluster_cn
from padfoot.__version__ import __version__


logger = logging.getLogger()

def _version():
    return __version__

def _check_python_dependencies():
    required_modules = {
        "pysam": "pysam",
        "pandas": "pandas",
        "numpy": "numpy",
        "Bio": "biopython",
    }

    missing = []
    for module_name, package_name in required_modules.items():
        try:
            importlib.import_module(module_name)
        except ImportError:
            missing.append(package_name)

    if missing:
        logger.error(
            "Missing required Python packages: %s",
            ", ".join(sorted(set(missing)))
        )
        logger.error(
            "Install them with conda, for example:\n"
            "  conda env create -f environment.yml\n"
            "or\n"
            "  conda install -c bioconda -c conda-forge %s",
            " ".join(sorted(set(missing)))
        )
        sys.exit(1)



def _check_external_dependencies(args, ref_path=None):
    required_bins = ["samtools", "minimap2"]
    missing_bins = [tool for tool in required_bins if shutil.which(tool) is None]

    if missing_bins:
        logger.error(
            "Missing required external tools in PATH: %s",
            ", ".join(missing_bins)
        )
        logger.error(
            "Please install them in your conda environment, e.g.:\n"
            "  conda install -c bioconda -c conda-forge %s",
            " ".join(missing_bins)
        )
        sys.exit(1)
    
    if shutil.which("RepeatMasker") is None:
        if args.run_repeatmasker:
            logger.error(
                "RepeatMasker is required for insertion repeat annotation but was not found in PATH. "
                "Install and configure RepeatMasker/FamDB, or use --skip_RepeatMasker."
            )
            sys.exit(1)
        logger.info("RepeatMasker is not found; insertion repeat annotation was skipped")

    if ref_path is not None:
        if not os.path.exists(ref_path):
            logger.error("Reference FASTA not found: %s", ref_path)
            sys.exit(1)

        fai_path = ref_path + ".fai"
        if not os.path.exists(fai_path):
            logger.error(
                "Reference FASTA index not found: %s\n"
                "Create it with:\n"
                "  samtools faidx %s",
                fai_path,
                ref_path,
            )
            sys.exit(1)
            
def _enable_logging(log_file, debug, overwrite):
    log_formatter = logging.Formatter(
        "[%(asctime)s] %(name)s: %(levelname)s: %(message)s",
        "%Y-%m-%d %H:%M:%S"
    )
    console_formatter = logging.Formatter(
        "[%(asctime)s] %(levelname)s: %(message)s",
        "%Y-%m-%d %H:%M:%S"
    )

    logger.handlers.clear()

    console_log = logging.StreamHandler()
    console_log.setFormatter(console_formatter)
    console_log.setLevel(logging.DEBUG if debug else logging.INFO)

    if overwrite:
        open(log_file, "w").close()

    file_handler = logging.FileHandler(log_file, mode="a")
    file_handler.setFormatter(log_formatter)

    logger.setLevel(logging.DEBUG)
    logger.addHandler(console_log)
    logger.addHandler(file_handler) 
    
    
def main():
    parser = argparse.ArgumentParser(
        description="Annotate structural variants and copy-number alterations"
    )
    parser.add_argument("-v", "--version", action="version", version=_version())
    parser.add_argument("--severus-vcf", "--sv-vcf", dest="vcf_file", metavar="path", required=True,
                        help="path to SV VCF file")
    parser.add_argument("--sv-caller", dest="sv_caller", choices=["severus", "savana"], default="severus",
                        help="SV caller input format [severus]")
    parser.add_argument("--wakhan-vcf", "--cna-file", dest="cna_vcf", metavar="path", required=True,
                        help="path to CNA VCF/TSV file")
    parser.add_argument("--cna-caller", dest="cna_caller", choices=["wakhan", "savana"], default="wakhan",
                        help="CNA caller input format [wakhan]")
    parser.add_argument("--ref", "-r", dest="ref", metavar="path", required=True,
                        help="path to reference file")
    parser.add_argument("--out-dir", dest="out_dir", required=True, metavar="path",
                        help="Output directory")
    parser.add_argument("--gff", dest="new_gff", default=None, metavar="path",
                        help="GFF annotation file")
    parser.add_argument("--genome", dest="genome", default="hg38", metavar="name",
                        help="Either hg38, chm13, mm10")
    parser.add_argument("--rm", dest="rm_file", default=None, metavar="path",
                        help="Repeat masker annotation file")
    parser.add_argument("-t", "--threads", dest="threads", default=8, type=int, metavar="int",
                        help="number of parallel threads [8]")
    parser.add_argument("--specie", dest="specie", default="human", help="Specie")
    parser.add_argument("--skip_RepeatMasker", dest="run_repeatmasker", action = "store_false", help="Skip RepeatMasker [True]")
    parser.add_argument("--ploidy", dest="ploidy", type=float, default=None, metavar="float",
                        help="tumour ploidy (mean total copy number) used as the copy-number baseline: a haplotype or allele "
                             "is AMP above round(ploidy/2) and DEL below it [estimated from the CN profile when absent]")
    parser.add_argument("--ploidy-file", dest="ploidy_file", default=None, metavar="path",
                        help="TSV with a 'ploidy' column to read the tumour ploidy from: SAVANA *_fitted_purity_ploidy.tsv or "
                             "Wakhan solutions_ranks.tsv (the rank-1 row is used); ignored when --ploidy is given")
    args = parser.parse_args()

    beds = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "beds"))

    args.out_dir = os.path.abspath(args.out_dir)
    args.vcf_file = os.path.abspath(args.vcf_file)
    args.cna_vcf = os.path.abspath(args.cna_vcf)
    args.ref = os.path.abspath(args.ref)
    if args.ploidy_file:
        args.ploidy_file = os.path.abspath(args.ploidy_file)

    os.makedirs(args.out_dir, exist_ok=True)

    temp_dir = os.path.join(args.out_dir, "temp")
    os.makedirs(temp_dir, exist_ok=True)
    os.chdir(temp_dir)

    log_file = os.path.join(args.out_dir, "padfoot.log")
    _enable_logging(log_file, debug=True, overwrite=True)
    _check_python_dependencies()
    _check_external_dependencies(args, args.ref)

    logger.info("Starting Padfoot %s", _version())
    logger.debug("Command: %s", " ".join(sys.argv))
    logger.debug("Python version: %s", sys.version)
    logger.info("Working directory: %s", os.getcwd())
    logger.info("Temporary directory: %s", temp_dir)
    logger.info("Output directory: %s", args.out_dir)

    logger.info("Input SV VCF (%s): %s", args.sv_caller, args.vcf_file)
    logger.info("Input CNA file (%s): %s", args.cna_caller, args.cna_vcf)
    logger.info("Reference FASTA: %s", args.ref)
    logger.info("Threads: %s", args.threads)
    logger.info("Requested species: %s", args.specie)
    logger.info("Requested genome preset: %s", args.genome)

    if not os.path.exists(args.vcf_file):
        logger.error("SV VCF does not exist: %s", args.vcf_file)
        sys.exit(1)

    if not os.path.exists(args.cna_vcf):
        logger.error("CNA file does not exist: %s", args.cna_vcf)
        sys.exit(1)

    if not os.path.exists(args.ref):
        logger.error("Reference FASTA does not exist: %s", args.ref)
        sys.exit(1)

    if not os.path.exists(args.ref + ".fai"):
        logger.error("Reference FASTA index is missing: %s.fai", args.ref)
        logger.error("Create it with: samtools faidx %s", args.ref)
        sys.exit(1)

    logger.info("Reference FASTA index found: %s.fai", args.ref)

    if args.ploidy is not None:
        if not args.ploidy > 0:
            logger.error("--ploidy must be a positive number: %s", args.ploidy)
            sys.exit(1)
        logger.info("Tumour ploidy from --ploidy: %.3f", args.ploidy)
    elif args.ploidy_file:
        if not os.path.exists(args.ploidy_file):
            logger.error("Ploidy file does not exist: %s", args.ploidy_file)
            sys.exit(1)
        try:
            args.ploidy = read_ploidy_file(args.ploidy_file)
        except (ValueError, KeyError, OSError) as exc:
            logger.error("Cannot read the tumour ploidy from %s: %s", args.ploidy_file, exc)
            sys.exit(1)
        logger.info("Tumour ploidy from %s: %.3f", args.ploidy_file, args.ploidy)
    else:
        logger.warning("No --ploidy / --ploidy-file given: the tumour ploidy will be estimated from the copy-number profile")

    if args.new_gff:
        gff_path = os.path.join(args.out_dir, "gff_file.gff3")
        logger.info("Using custom GFF input: %s", args.new_gff)
        logger.debug("Generating normalized GFF file: %s", gff_path)
        generate_gff(args.new_gff, gff_path)
        args.gff_file = gff_path
        logger.info("Prepared GFF annotations: %s", args.gff_file)
    else:
        if not args.genome:
            logger.error("Please provide --gff or select --genome [hg38, mm10, chm13]")
            sys.exit(1)
        args.gff_file = os.path.join(beds, f"{args.genome}.gff3.gz")
        logger.info("Using bundled GFF annotations: %s", args.gff_file)

    if args.rm_file:
        rm_path = os.path.join(args.out_dir, "rm.bed")
        logger.info("Using custom RepeatMasker input: %s", args.rm_file)
        logger.debug("Generating normalized RepeatMasker BED: %s", rm_path)
        generate_rm(args.rm_file, rm_path)
        args.rm_file = rm_path
        logger.info("Prepared RepeatMasker annotations: %s", args.rm_file)
    else:
        if not args.genome:
            logger.error("Please provide --rm or select --genome [hg38, mm10, chm13]")
            sys.exit(1)
        args.rm_file = os.path.join(beds, f"{args.genome}_rm.bed.gz")
        logger.info("Using bundled RepeatMasker annotations: %s", args.rm_file)

    if not args.specie:
        logger.error("Please provide --specie for repeat annotation")
        sys.exit(1)

    if args.genome == "mm10" and args.specie != "mouse":
        logger.warning("Genome preset mm10 usually implies species 'mouse'; current value: %s", args.specie)

    logger.info("Annotation inputs resolved successfully")
    logger.info("GFF file: %s", args.gff_file)
    logger.info("RepeatMasker file: %s", args.rm_file)

    logger.info("Starting annotation pipeline")
    genes, cnas, exon_pos, svs, by_gene = annotate_things(args)

    logger.info("Annotation pipeline finished")
    logger.info("Summary: genes=%d, cnas=%d, exon_positions=%d, svs=%d, by_gene=%d",
                len(genes), len(cnas), len(exon_pos), len(svs), len(by_gene))
    logger.info("Log file written to: %s", log_file)
    logger.info("Padfoot completed successfully")