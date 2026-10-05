# Bundled annotation files

| file | content | source |
| --- | --- | --- |
| `hg38.gff3.gz` | protein-coding genes, transcripts, exons and UTRs of GENCODE v47 basic (GRCh38), in Padfoot's 6-column format (`chr gene strand type start end`), as written by `padfoot/preprocess.py::generate_gff` | GENCODE v47 basic annotation |
| `mm10.gff3.gz` | the same for mouse (GRCm38/mm10) | GENCODE mouse basic annotation |
| `hg38_rm.bed.gz`, `mm10_rm.bed.gz` | RepeatMasker annotations (`chr start end class/family`) | UCSC RepeatMasker tracks |
| `cancer_genes.tsv` | cancer gene roles and fusion partners (`Gene_symbol`, `Role_in_cancer`, `fusion`, `Supp`, `nSupp`), used with `--specie human`; override with `--cancer-genes` | copy of <https://raw.githubusercontent.com/aysegokce/PadfootFiles/refs/heads/main/cancer_genes.tsv>, PadfootFiles commit `92c1842e` (2025-07-11), retrieved 2026-10-02, sha256 `1fa94960d5a4fe27` (first 16 hex digits of the file hash recorded at import) |

Notes

- Before 2026-10-02 `cancer_genes.tsv` was downloaded from GitHub on every human run; it is now read from this directory so
  Padfoot runs without network access. Padfoot `--cancer-genes <path>` points it at a different table.
- `generate_gff()` used to write the last gene of each chromosome under the next chromosome's name (and dropped the last gene
  of the file). The two bundled `.gff3.gz` files were corrected for the first defect on 2026-10-02: every chromosome block's
  first gene whose coordinates lie beyond the chromosome or after the block's second gene was moved back to the previous
  chromosome (hg38: 24 blocks / 348 rows -- PGBD2, FAM240C, LMLN, DUX4, OR4F3, PDCD2, VIPR2, C8orf33, CACNA1B, FRG2B, B3GAT1,
  ANHX, CHAMP1, TMEM121, OR4F4, PRDM7, METRNL, PARD6G, MZF1, PCMTD2, PRMT2, RABL2B, WASH6P on chrX and on chrY; mm10: 21 blocks /
  216 rows -- Cr2 ... Gm21748). Only the chromosome field changed. The second defect cannot be repaired without regenerating
  the files: the last mitochondrial gene (hg38 `MT-CYB`, mm10 `mt-Cytb`) is absent from both. Regenerating with the fixed
  `generate_gff()` from the GENCODE GFF3 restores it.
- `padfoot/annot.py::get_genes` sorts each chromosome's gene list by (padded) start and `annotBPs` finds the gene containing a
  breakpoint by a containment search, so the result no longer depends on the file order of the genes.
