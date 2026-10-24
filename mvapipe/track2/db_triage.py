"""Triage of the database screen with explicit, auditable rules (a review list for a human, not a recommendation).

Input: db_drug_gene.csv from db_screen.py. Per drug, interactions are grouped and rules applied in this order:

  evidence   an interaction counts only if it is *typed* (DGIdb gives an interaction type) or is a clinical-stage record
             in Open Targets. Untyped DGIdb rows are largely text-mined and are ignored (reported as a count).
  STOP  mitotic machinery   typed/clinical interaction with a spindle-checkpoint/mitotic gene (wrong direction in a
                            checkpoint-weakened cell). A STRICT variant also counts untyped rows.
  STOP  not approved        no approval flag (DGIdb) and no Phase 4/approval record (Open Targets)
  STOP  opposite direction  axis-specific: immune/cGAS-STING axis hit by an agonist/activator; NAD+ axis hit by an inhibitor
  STOP  no primary target   no typed/clinical interaction with a primary target of the mechanism axes (readout genes such
                            as TP53, CDKN1A/2A, IL6, MMP3 are not primary targets)
  REVIEW direction unknown  passes the stops but no typed direction (e.g. biologics recorded without a type)
  PASS                      typed interaction in the wanted direction, approved
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd

from .db_screen import GENE_SETS

MITOTIC = set(GENE_SETS["spindle_checkpoint"]) - {"BUB1B"}
PRIMARY = {"interferon/JAK": {"JAK1", "JAK2", "TYK2", "IFNAR1", "IFNAR2", "STAT1", "STAT2"},
           "cGAS-STING": {"CGAS", "STING1", "TBK1", "IRF3"},
           "senescence": {"BCL2", "BCL2L1", "SRC"},
           "NAD+/sirtuin": {"SIRT1", "SIRT2", "NAMPT", "NMNAT1", "NAPRT", "PARP1"}}
WANT_INHIBIT = {"interferon/JAK", "cGAS-STING", "senescence"}
RAISE_TYPES = ("agonist", "activator", "positive modulator", "vaccine")
LOWER_TYPES = ("inhibitor", "antagonist", "blocker", "negative modulator")


def _approved(row) -> bool:
    a = row.get("approved")
    if a is True or str(a).lower() == "true":
        return True
    return str(row.get("stage", "")).upper() in {"APPROVAL", "4", "4.0", "PHASE_4", "PHASE IV"}


def _has(t, words) -> bool:
    return any(w in str(t).lower() for w in words)


def triage(df: pd.DataFrame, kb: dict, strict_mitotic: bool = False) -> pd.DataFrame:
    df = df.copy()
    df["drug"] = df["drug"].astype(str).str.upper()
    df["typed"] = df["type"].notna() | df["stage"].notna()
    df["is_approved"] = df.apply(_approved, axis=1)
    kb_names = {d.upper(): i["id"] for i in kb["interventions"] for d in i.get("drugs", [])}
    kb["screen_review"] = {k.upper(): v for k, v in kb.get("screen_review", {}).items()}
    rows = []
    for drug, g in df.groupby("drug"):
        ev = g if strict_mitotic else g[g["typed"]]
        reasons = []
        mit = set(ev["gene"]) & MITOTIC
        if mit:
            reasons.append("acts on mitotic machinery: " + ",".join(sorted(mit)))
        if not g["is_approved"].any():
            reasons.append("not approved")
        axes, wanted, opposite = set(), False, False
        for axis, genes in PRIMARY.items():
            sub = g[g["typed"] & g["gene"].isin(genes)]
            if sub.empty:
                continue
            axes.add(axis)
            ty = sub["type"].dropna()
            if axis in WANT_INHIBIT:
                opposite |= any(_has(t, RAISE_TYPES) for t in ty)
                wanted |= any(_has(t, LOWER_TYPES) for t in ty)
            else:
                opposite |= any(_has(t, LOWER_TYPES) for t in ty)
                wanted |= any(_has(t, RAISE_TYPES) for t in ty)
        if opposite and not wanted:
            reasons.append("acts opposite to the wanted direction on " + "/".join(sorted(axes)))
        if not axes:
            reasons.append("no typed interaction with a primary target")
        status = "STOP" if reasons else ("PASS" if wanted else "REVIEW: direction unknown")
        rows.append({"drug": drug, "status": status, "stop": "; ".join(reasons) or None, "axes": ";".join(sorted(axes)),
                     "genes": ",".join(sorted(set(ev["gene"]))), "approved": bool(g["is_approved"].any()),
                     "n_untyped_rows_ignored": int((~g["typed"]).sum()),
                     "manual_review": kb.get("screen_review", {}).get(drug), "in_knowledge_file": next((v for k, v in kb_names.items() if drug == k or drug.startswith(k + " ")), None)})
    return pd.DataFrame(rows).sort_values(["status", "axes", "drug"])


def main(argv=None) -> int:
    from .run import load_kb
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--csv", required=True); ap.add_argument("--out", default="work_track2")
    a = ap.parse_args(argv)
    df, kb = pd.read_csv(a.csv), load_kb()
    out = Path(a.out); out.mkdir(parents=True, exist_ok=True)
    res, strict = triage(df, kb), triage(df, kb, True)
    res.to_csv(out / "db_triage.csv", index=False)
    summary = {"rows": int(len(df)), "rows_untyped_ignored": int((~(df["type"].notna() | df["stage"].notna())).sum()),
               "drugs": int(len(res)), "status_counts": res["status"].value_counts().to_dict(),
               "status_counts_strict_mitotic_rule": strict["status"].value_counts().to_dict(),
               "stopped_mitotic": int(res["stop"].fillna("").str.contains("mitotic").sum()),
               "stopped_opposite": int(res["stop"].fillna("").str.contains("opposite").sum()),
               "knowledge_file_drugs": {r.drug: r.status for r in res.itertuples() if pd.notna(r.in_knowledge_file)}}
    (out / "db_triage_summary.json").write_text(json.dumps(summary, indent=1))
    print(json.dumps(summary, indent=1))
    print(res[res["status"] != "STOP"].to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
