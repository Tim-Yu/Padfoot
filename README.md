# Padfoot
<p>
<img src="docs/logo.png" alt="Padfoot logo" align="left" style="width:100px;"/>
<b>Padfoot</b> is a structural variant (SV) and copy number alteration (CNA) annotation tool designed to work seamlessly with output from <b>Severus</b>, <b>Wakhan</b>, and <b>Savana</b>. It provides biological and functional annotations for SVs and CNAs, enabling downstream interpretation of genomic alterations.
</p>

<br/>


## Features

- Annotates somatic SVs called by [Severus](https://github.com/KolmogorovLab/Severus)
- Annotates CNAs reported by [Wakhan](https://github.com/KolmogorovLab/Wakhan)
- Annotates Savana classified somatic SV VCFs and segmented absolute copy-number TSVs
- Outputs functional annotations for prioritization
- Designed for cancer genome analysis using long-read data



## Contents

* [Installation](#installation)
* [Quick Usage](#quick-usage)
* [Input and Parameters](#inputs-and-parameters)


## Installation

The easiest way to install is through conda:

```
conda create -n padfoot_env padfoot
conda activate padfoot_env
padfoot --help
```

Or alternatively, you can clone the repository and run without installation,
but you'll still need to install the dependencies via conda:

```
git clone https://github.com/KolmogorovLab/Padfoot.git
cd Padfoot
conda env create --name padfoot_env --file environment.yml
conda activate padfoot_env
./padfoot.py
```

## Quick Usage

```
padfoot --severus-vcf severus_somatic.vcf --wakhan-vcf wakhan_cna.vcf --ploidy-file solutions_ranks.tsv --ref ref.fa --out-dir padfoot_out -t 16 --specie human
```

Savana input can be run with caller-specific format flags:

```
padfoot --sv-vcf sample.classified.somatic.vcf --sv-caller savana --cna-file sample_segmented_absolute_copy_number.tsv --cna-caller savana --ploidy-file sample_fitted_purity_ploidy.tsv --ref ref.fa --out-dir padfoot_savana_out -t 16 --specie human
```

## Inputs and Parameters

### Required

```
--severus-vcf   path to Severus vcf
--wakhan-vcf    path to Wakhan vcf
--sv-vcf        generic alias for SV VCF input
--sv-caller     SV input format: severus or savana [severus]
--cna-file      generic alias for CNA VCF/TSV input
--cna-caller    CNA input format: wakhan or savana [wakhan]
--ref           path to reference fasta file (needs to be indexed)
--out-dir       path to output directory
```

### Optional parameters

```
--threads               number of threads [8]
--ploidy                tumour ploidy (mean total copy number); genes are AMP above round(ploidy/2) copies per
                        haplotype/allele and DEL below it [estimated from the CN profile when absent, with a warning]
--ploidy-file           read the ploidy from the caller's fit table instead: SAVANA *_fitted_purity_ploidy.tsv or
                        Wakhan solutions_ranks.tsv (rank-1 row); ignored when --ploidy is given
--specie                human or mouse or user defined(GFF and rm files need to be provided) [human]
--genome                Either hg38, chm13 or mm10 [hg38]
--gff                   If user want to use a alternative gff
--rm                    Repeat masker file (.fa.out)
```

### Genome Annotations

We provide default **GFF** and **RepeatMasker** annotation files in the `bed/` directory.

If you would like to use a different genome assembly or custom annotations, you can specify your own files using the following parameters:

- `--gff`: Path to a custom GFF annotation file  
- `--rm`: Path to a custom RepeatMasker BED file  
- `--specie`: Specify the species or genome build (e.g., `hg38`, `chm13`)

These options allow flexibility for applying **Padfoot** to a variety of genome references and annotation sources.

### 📁 Output

The output is a tab-delimited file containing:

- Annotated SVs and CNAs
- Gene overlaps and exon-level effects  
- Functional impact scores  
- Complex SV plots


License
-------

Padfoot is distributed under a BSD license. See the [LICENSE file](LICENSE) for details.


Credits
-------

Padfoot was developed by the Kolmogorov Lab at the National Cancer Institute as part of ongoing efforts to improve interpretation of long-read somatic variant data in cancer genomics.

Key contributors:

* Ayse Keskus
* Mikhail Kolmogorov

---
### Contact
For advising, bug reporting and requiring help, please submit an [issue](https://github.com/KolmogorovLab/Padfoot/issues).
You can also contact the developer: aysegokce.keskus@nih.gov
