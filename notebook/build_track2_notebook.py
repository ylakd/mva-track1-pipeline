#!/usr/bin/env python3
"""Generate notebook/mva_track2_framework.ipynb (embeds the Track 2 package and knowledge files via %%writefile)."""
from pathlib import Path
import nbformat as nbf

ROOT = Path(__file__).resolve().parents[1]
PKG = ["__init__.py", "faulttree.py", "fmea.py", "readthrough.py", "scenarios.py", "checkpoint_sim.py", "db_screen.py", "run.py"]
KB = ["bub1b_mva1_kb.json", "TEMPLATE_new_genotype.json"]
cells = []
md = lambda t: cells.append(nbf.v4.new_markdown_cell(t.strip()))
code = lambda t: cells.append(nbf.v4.new_code_cell(t.strip()))

md("""
# MVA Hackathon 2026 — Track 2 · fault-tree / FMEA triage of repurposing hypotheses (yasld)
Runs without patient data: it reads a curated knowledge file (`knowledge/*.json`). Optional last step: if the Track 1
output files `nmd_check.json` and `locus_summary.json` are present in `TRACK1_WORK`, the nonsense allele is checked for
translational-readthrough context (needs Internet for two reference bases from the UCSC API).
Outputs are **hypotheses for follow-up, not treatment advice**. Severity and detection scores are expert judgements; see the rank-stability table.
""")
code('!pip install -q pandas matplotlib numpy\nimport os, sys, json\nsys.path.insert(0, os.getcwd())\nos.makedirs("mvapipe/track2", exist_ok=True); os.makedirs("knowledge", exist_ok=True)\nopen("mvapipe/__init__.py", "a").close()')
for n in PKG:
    code(f"%%writefile mvapipe/track2/{n}\n" + (ROOT / "mvapipe" / "track2" / n).read_text())
for n in KB:
    code(f"%%writefile knowledge/{n}\n" + (ROOT / "knowledge" / n).read_text())
md("## 1. Fault trees and minimal cut sets\nThe root tree is AND(allele 1 fails, allele 2 fails), each allele an OR of mechanisms. The consequence tree is an OR of downstream harms.")
code('''import pandas as pd
from mvapipe.track2 import run as t2, faulttree as ft
kb = t2.load_kb("knowledge/bub1b_mva1_kb.json")
res = t2.analyse(kb)
print("root cut sets:", res["minimal_cut_sets"]["root"])
print("smallest sets that break every root cut set:", res["fully_mitigating_sets_root"])''')
md("## 2. FMEA ranking (lower RPN = better) and hard stops")
code('pd.set_option("display.width", 200); pd.set_option("display.max_colwidth", 70)\ndisplay(pd.DataFrame(res["ranking"]))\ndisplay(pd.DataFrame(res["excluded"])[["id", "label", "stop"]])')
md("## 3. Rank stability\nEvery S, O and D is perturbed by -1, 0 or +1 (10,000 draws). Reports how often the order or the first place survives.")
code('print(json.dumps(res["rank_stability"], indent=1))')
md("## 4. Figures and files")
code('''TRACK1_WORK = os.environ.get("MVA_WORKDIR", "work_track2")   # point to the Track 1 work directory to enable the readthrough step
res = t2.run(outdir=TRACK1_WORK, network=True)
print(json.dumps(res.get("readthrough_check", "no Track 1 outputs found: readthrough step skipped"), indent=1))
print(sorted(os.listdir(TRACK1_WORK)))''')
md("## 5. Genotype-uncertainty scenarios\nr = probability that the compound-heterozygous genotype is the causal one. Root-tier candidates lose value as r falls; consequence-tier ones do not. A second tumour-predisposition variant raises severity for tagged interventions.")
code('print(json.dumps(res["scenarios"], indent=1, default=str))')
md("## 6. Spindle-checkpoint simulation (toy model, direction of effect only)\nNot fitted to data. Reports the share of 5,000 random parameter sets in which each mitosis-acting intervention gives fewer, equal or more uncorrected chromosome errors in a checkpoint-weakened cell.")
code('print(json.dumps(res["checkpoint_simulation"], indent=1, default=str)[:3500])')
md("## 7. Database and pathway screen (needs Internet; only gene symbols are sent)\nDGIdb, Open Targets and Enrichr for the mechanism's gene sets. Failures are listed, not hidden. Paste the printed table back for review.")
code('import subprocess; print(subprocess.run([sys.executable, "-m", "mvapipe.track2.db_screen", "--out", "work_track2"], capture_output=True, text=True).stdout[-6000:])')
md("## 8. Same engine, a different genotype class\nThe template knowledge file is synthetic; it shows that a new disease needs a new JSON file, not new code.")
code('kb2 = t2.load_kb("knowledge/TEMPLATE_new_genotype.json")\ndisplay(pd.DataFrame(t2.analyse(kb2)["ranking"]))')

nb = nbf.v4.new_notebook(cells=cells)
nb.metadata = {"kernelspec": {"display_name": "Python 3", "language": "python", "name": "python3"},
               "language_info": {"name": "python"}}
out = ROOT / "notebook" / "mva_track2_framework.ipynb"
nbf.validate(nb); nbf.write(nb, out)
print("written", out, "-", len(cells), "cells")
