# MVA Hackathon 2026 — Track 1 pipeline

Phenotype-guided, selection-corrected variant prioritisation for one proband
(*Rare Disease, Real Kid: MVA Hackathon 2026*, Track 1). Author: **yasld**.

Input: the proband's pre-called VCF and eight HPO-coded clinical features. Output: a ranked prediction CSV in the
challenge format. Designed for N = 1: no model is trained; ~4.7 M variants are triaged with explicit, pre-specified
filters, genes are scored by phenotype fit weighted by term specificity, and the result is tested against a permutation
null that corrects for choosing the top gene after the fact.

> **Status.** Code consolidated from an interactive analysis and covered by synthetic tests (`tests/`). The numbers in the
> report come from the interactive run on the real data. Run the full pipeline and compare with `check_results.py`
> before relying on them.

## Stages

| # | Module | Purpose | Output |
|---|---|---|---|
| 1 | `stage01_03_data_qc_coding` | Download VCF/index/phenotype doc (HF), Ensembl GTF, HPO table; sizes + hashes | `provenance.json` |
| 2 | same | QC: PASS, non-reference, DP>=10, GQ>=30; **every exclusion counted** | `variants_pass.parquet`, `qc_stats.json` |
| 3 | same | Coding restriction (CDS ±2 bp, merged across isoforms) | `variants_coding.parquet` |
| 4 | `stage04_05_annotate_scan` | MyVariant.info annotation (gnomAD, ClinVar, effect); raw JSON cached | `myvariant_raw.json` |
| 5 | same | Local gene assignment from GTF; **pre-specified recessive scan** (A homozygous, B >=2 het in a gene); ClinVar check | `scan_*.csv` |
| 6 | `stage06_phenotype` | Specificity-weighted ranking, leave-one-symptom-out, **permutation null**, rule-based selection | `phenotype_ranking.csv`, `permutation_results.json`, `selected_candidates.json` |
| 7 | `stage07_verify` | Calls read from the VCF (allele balance, phasing tags), markers between pair members, ClinVar, predictors, UniProt | `verification.json` |
| 8 | `stage08_locus` | **Supp.** Unfiltered review of the whole candidate locus | `locus_summary.json` |
| 9 | `stage09_10_nmd_gnomad` | **Supp.** Exon number + NMD rule for truncating candidates | `nmd_check.json` |
| 10 | same | **Supp.** gnomAD v4 lookup (`unavailable` is NOT absence) | `gnomad_check.json` |
| 11 | `stage11_baf` | **Supp.** Allele-fraction + depth imbalance screen, with detection limit | `baf_summary.json` |
| 12 | `stage12_submission` | Write + validate the submission CSV (a pair is ONE row) | `submission.csv` |

Candidate selection is rule-based and contains **no gene names**.

## Run

Kaggle: accept the dataset terms on Hugging Face, add secret `HF_TOKEN`, turn Internet on, open
`notebook/mva_track1_pipeline.ipynb` (self-contained) and Run all.

Local:
```bash
pip install -r requirements.txt
export HF_TOKEN=hf_...
python run_pipeline.py                 # stages 1-12 (outputs in ./work or $MVA_WORKDIR)
python run_pipeline.py --stages 8-11   # supplementary analyses only
python check_results.py                # compare with the reported numbers
python -m unittest discover -s tests -v
```
Pin releases for a frozen run: `MVA_GTF_URL`, `MVA_HPO_URL`. `provenance.json` and `environment.json` record versions.

## Design decisions
* Pre-specified filters: impact high/moderate; gnomAD AF < 0.1 % **or no value** (kept, flagged); no benign ClinVar entry.
* Genes assigned locally from the GTF (the web service covered <50 % of coding variants).
* dbSNP membership is not used as a rarity or pathogenicity proxy.
* Damage weights are **not blind** (fixed after the candidate was known); the permutation null mitigates, not removes, this.
* The depth/allele-fraction screen cannot exclude mosaic variegated aneuploidy; it reports what a bulk test could detect.
* Undeclared FORMAT tags (PGT/PID in an unphased callset) are handled explicitly.

## Limitations
Coding sequence only in the main scan; no SV/CNV calling; single sample (phase unresolved); gnomAD API may be
unavailable; annotation depends on MyVariant.info.

## Privacy
The proband is a real child; the challenge forbids publishing anything that could re-identify the family. Do not commit
the VCF, FASTQ, parquet files, raw annotation JSON, locus table or phenotype document (`.gitignore` excludes them).
The proband's HPO terms are not in the code: stage 1 reads them from the challenge phenotype document into the git-ignored working directory. Check the challenge rules before publishing variant-level output. This repository contains code and synthetic test data only.

## Licence and credits
Add a licence before making the repository public. Data: Sage Bionetworks (CC BY 4.0). Resources: Ensembl, MyVariant.info,
UniProt, Human Phenotype Ontology.
