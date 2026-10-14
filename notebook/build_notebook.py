#!/usr/bin/env python3
"""Generate notebook/mva_track1_pipeline.ipynb (embeds the package via %%writefile cells).
Regenerate after any code change:  python notebook/build_notebook.py"""
from pathlib import Path
import nbformat as nbf

ROOT = Path(__file__).resolve().parents[1]
PKG = ["__init__.py", "config.py", "utils.py", "stage01_03_data_qc_coding.py", "stage04_05_annotate_scan.py",
       "stage06_phenotype.py", "stage07_verify.py", "stage08_locus.py", "stage09_10_nmd_gnomad.py",
       "stage11_baf.py", "stage12_submission.py"]
cells = []
md = lambda t: cells.append(nbf.v4.new_markdown_cell(t.strip()))
code = lambda t: cells.append(nbf.v4.new_code_cell(t.strip()))

md("""
# MVA Hackathon 2026 — Track 1 · full pipeline (yasld)
**Phenotype-guided, selection-corrected variant prioritisation for one proband.** Self-contained: the first cells write the
`mvapipe` package to disk, then each stage runs in order (stages 8–11 are supplementary analyses).

**Antes de ejecutar:** (1) acepta las condiciones del dataset `SageBio/mva-hackathon-2026-data` en Hugging Face;
(2) Kaggle → *Add-ons → Secrets*: `HF_TOKEN` (token de lectura) y **Internet ON**; (3) *Run all*. Las etapas 2 y 4 tardan varios
minutos; cada etapa guarda su salida en disco, así que se puede reanudar. Si una etapa falla, relanza solo esa celda.

**Datos protegidos.** Es un niño real: no publiques el VCF, los parquet, el JSON de anotación ni la tabla del locus.
""")
md("## 0. Setup")
code('!pip install -q cyvcf2 python-docx huggingface_hub pyarrow\nimport os, sys\nos.makedirs("mvapipe", exist_ok=True)\nsys.path.insert(0, os.getcwd())')
for name in PKG:
    code(f"%%writefile mvapipe/{name}\n" + (ROOT / "mvapipe" / name).read_text())
code("%%writefile expected_results.json\n" + (ROOT / "expected_results.json").read_text())
code("%%writefile check_results.py\n" + (ROOT / "check_results.py").read_text())
code('''
import json
import pandas as pd
pd.set_option("display.width", 200); pd.set_option("display.max_columns", 30); pd.set_option("display.max_colwidth", 80)
from mvapipe import config as cfg
from mvapipe.utils import hf_token, read_json
import mvapipe.stage01_03_data_qc_coding as s1
import mvapipe.stage04_05_annotate_scan as s4
import mvapipe.stage06_phenotype as s6
import mvapipe.stage07_verify as s7
import mvapipe.stage08_locus as s8
import mvapipe.stage09_10_nmd_gnomad as s9
import mvapipe.stage11_baf as s11
import mvapipe.stage12_submission as s12
assert hf_token(), "Add HF_TOKEN as a Kaggle secret (and accept the dataset terms) before continuing"
print("working directory:", cfg.workdir())
''')
md("## 1. Data and provenance\nVCF, index and phenotype document from Hugging Face; Ensembl GTF (GRCh38, r112); HPO gene–phenotype table. Records sizes and SHA-256 and checks the HPO terms against the document.")
code('prov = s1.stage01_data()\nprint("HPO terms match the document:", prov.get("hpo_terms_match_config"))\npd.DataFrame(prov["files"]).T')
md("## 2. Quality control\nPASS, non-reference, DP ≥ 10, GQ ≥ 30. **Every exclusion is counted** (the interactive run recorded only the final count).")
code("pd.Series(s1.stage02_qc())")
md("## 3. Coding restriction\nCDS ±2 bp (canonical splice sites), merged across isoforms.")
code("pd.Series(s1.stage03_coding())")
md("## 4. Annotation\nBatched MyVariant.info queries (hg38). Multi-allelic records are not annotated and are counted. Raw JSON is cached.")
code("pd.Series(s4.stage04_annotate())")
md("## 5. Phenotype-agnostic recessive scan (pre-specified filters)\nImpact high/moderate; gnomAD AF < 0.1 % **or no value** (kept, flagged); no benign ClinVar. Genes assigned locally from the GTF. A = homozygous; B = ≥2 heterozygous in one gene. Check X-linked 'homozygous' calls (hemizygous in XY individuals) and genes whose pair is an indel cluster (possible artefact).")
code('''scan = s4.stage05_scan()
print(json.dumps(scan, indent=2))
print("\\n--- A: homozygous ---"); display(pd.read_csv(cfg.path("scan_A")))
print("--- B: >=2 heterozygous in one gene ---"); display(pd.read_csv(cfg.path("scan_B")))''')
md("## 6. Phenotype ranking, robustness and permutation null\nHPO terms weighted by IC = −log₂(genes with term / all genes). Gene score = phenotype fit × damage score. Leave-one-symptom-out and a permutation null (random terms of similar specificity): (i) chance of the top gene's fit, (ii) chance that **any** gene reaches the best score (corrects for post-hoc selection). Selection is rule-based (no gene names). *Caveat:* damage weights are not blind.")
code('''ph = s6.stage06_phenotype()
print(json.dumps({k: v for k, v in ph.items() if k != "primary"}, indent=2, default=str))
print("\\nTop of ranking:"); display(pd.read_csv(cfg.path("ranking")).head(10))
print("Leave-one-symptom-out:"); display(pd.read_csv(cfg.path("loo")))
print("Selected:"); print(json.dumps(read_json(cfg.path("selected")), indent=2))''')
md("## 7. Verification of the candidates\nCall read from the VCF (allele balance, depth, GQ, FILTER, phasing tags); heterozygous markers between the pair members and the largest gap (can short reads phase them?); ClinVar, predictors, UniProt features.")
code('''ver = s7.stage07_verify()
p = ver["primary"]
if p:
    print("Primary gene:", p["gene"])
    display(pd.DataFrame(p["vcf_calls"])[["chrom","pos","ref","alt","filter","genotype","phased","ad_ref","ad_alt","alt_fraction","dp","gq","has_PGT","has_PID"]])
    print(json.dumps(p["pair_geometry_and_phase"], indent=2))
    for a in p["annotation"]:
        print("\\n", a.get("hgvs"), "| gnomAD genome/exome AF:", a.get("gnomad_genome_af"), a.get("gnomad_exome_af"))
        print("  ClinVar:", json.dumps(a.get("clinvar"))[:600]); print("  predictors:", a.get("predictors"))
    for u in p["uniprot"]:
        print("\\nUniProt", u.get("accession"), "residue", u.get("residue")); display(pd.DataFrame(u.get("features_near_residue", [])))''')
md("## 8. Supplementary — unfiltered locus review\nEvery record in the gene body ±50 kb with **no quality filter**, classified (CDS, UTR exon, intron, flank). Is there a non-coding/splice-region second hit? Did a quality filter remove a coding call?")
code('''loc = s8.stage08_locus()
print(json.dumps({k: v for k, v in loc.items() if k not in ("coding_non_reference","rare_noncoding_candidates")}, indent=2, default=str))
print("\\nCoding records, unfiltered:"); display(pd.DataFrame(loc["coding_non_reference"]))
print("Rare non-coding candidates:"); display(pd.DataFrame(loc["rare_noncoding_candidates"]))''')
md("## 9. Supplementary — exon position and NMD rule\nMANE Select transcript. A premature stop >~50 nt upstream of the last exon–exon junction is expected to trigger NMD; a stop in the last exon is expected to escape it.")
code('''s9.stage09_nmd()
for r in read_json(cfg.path("nmd"))["results"]:
    print(r["variant"].get("gene"), r["variant"].get("hgvs_p"), "->", r.get("transcript"), "| exon", r.get("exon_number"), "of", r.get("total_exons"), "|", r.get("interpretation"))''')
md("## 10. Supplementary — gnomAD v4 lookup\n**`unavailable` means the query failed, not that the variant is absent.** Check such entries in the gnomAD browser and note the result by hand.")
code('print(s9.stage10_gnomad())\ndisplay(pd.json_normalize(read_json(cfg.path("gnomad"))["results"]))')
md("## 11. Supplementary — allele-fraction and depth screen\nHet-SNV allele fractions and depth per chromosome and 5 Mb window (robust z-scores) and the smallest *uniform* mosaic fraction detectable in depth. In MVA the aneuploidies differ cell to cell: a null result bounds a bulk test, it does **not** exclude the condition.")
code('bs = s11.stage11_baf()\nprint(json.dumps(bs, indent=2, default=str))\ndisplay(pd.read_csv(cfg.path("baf_chrom"), index_col=0).round(4))')
md("## 12. Submission file\nOne row per prediction; the compound-heterozygous pair is **one row**; GRCh38 with `chr` prefix. The validator checks the schema and that each variant exists in the VCF.")
code('sub = s12.stage12_submission()\nprint(sub)\ndisplay(pd.read_csv(cfg.path("submission")))')
md("## 13. Compare with the reported numbers")
code("!python check_results.py")
md("## 14. Files\nDownload `submission.csv`. **Do not make variant-level files public** (`variants_*.parquet`, `myvariant_raw.json`, `locus_review.csv`, `annotated.parquet`).")
code('for f in sorted(cfg.workdir().glob("*.json")) + sorted(cfg.workdir().glob("*.csv")):\n    print(f"{f.stat().st_size:>10,d}  {f.name}")')

nb = nbf.v4.new_notebook(cells=cells)
nb.metadata = {"kernelspec": {"display_name": "Python 3", "language": "python", "name": "python3"},
               "language_info": {"name": "python"}}
out = ROOT / "notebook" / "mva_track1_pipeline.ipynb"
nbf.validate(nb); nbf.write(nb, out)
print("written", out, "-", len(cells), "cells")
