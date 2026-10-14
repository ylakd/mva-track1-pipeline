"""Stage 12: write the Track 1 submission CSV and validate it.

Schema (from the challenge Space): one row per prediction; a compound-heterozygous PAIR is ONE row
(chrom_1..alt_1 and chrom_2..alt_2). Splitting a pair into two rows earns only partial credit.
Coordinates are GRCh38 with the 'chr' prefix; `epcr` in (0, 1] orders rows (highest first).
"""
from __future__ import annotations

import csv
import re

from . import config as cfg
from .utils import read_json

COLUMNS = ["proband_id", "chrom_1", "pos_1", "ref_1", "alt_1", "chrom_2", "pos_2", "ref_2", "alt_2",
           "epcr", "finding_type", "notes"]


def _chr(c: str) -> str:
    c = str(c)
    return c if c.startswith("chr") else "chr" + c


def build_rows(sel: dict) -> list[dict]:
    rows = []
    p = sel.get("primary")
    if p:
        v1, v2 = p["variants"][0], p["variants"][1]
        rows.append({"proband_id": cfg.PROBAND_ID,
                     "chrom_1": _chr(v1["chrom"]), "pos_1": v1["pos"], "ref_1": v1["ref"], "alt_1": v1["alt"],
                     "chrom_2": _chr(v2["chrom"]), "pos_2": v2["pos"], "ref_2": v2["ref"], "alt_2": v2["alt"],
                     "epcr": cfg.EPCR_PRIMARY, "finding_type": "primary",
                     "notes": f"{p['gene']} compound-heterozygous candidate: {v1.get('hgvs_p')} + "
                              f"{v2.get('hgvs_p')}. Phase unresolved."})
    for i, s in enumerate(sel.get("secondary", [])[: cfg.MAX_ROWS - len(rows)]):
        rows.append({"proband_id": cfg.PROBAND_ID,
                     "chrom_1": _chr(s["chrom"]), "pos_1": s["pos"], "ref_1": s["ref"], "alt_1": s["alt"],
                     "chrom_2": "", "pos_2": "", "ref_2": "", "alt_2": "",
                     "epcr": round(max(cfg.EPCR_SECONDARY - 0.01 * i, 0.01), 2), "finding_type": "secondary",
                     "notes": f"{s['gene']} {s.get('hgvs_p')} heterozygous; possible incidental finding."})
    return rows


def write_csv(rows: list[dict], path) -> None:
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=COLUMNS)
        w.writeheader()
        w.writerows(rows)


def validate_submission(path, vcf_path=None) -> list[str]:
    """Return a list of problems (empty = valid)."""
    problems: list[str] = []
    with open(path, newline="") as f:
        rd = csv.DictReader(f)
        missing = [c for c in COLUMNS[:10] if c not in (rd.fieldnames or [])]
        if missing:
            return [f"missing columns: {missing}"]
        rows = list(rd)
    if not rows:
        problems.append("no rows")
    if len(rows) > cfg.MAX_ROWS:
        problems.append(f"more than {cfg.MAX_ROWS} rows")
    prev, seen = 2.0, set()
    for i, r in enumerate(rows, 1):
        tag = f"row {i}"
        if r["proband_id"].strip() != cfg.PROBAND_ID:
            problems.append(f"{tag}: proband_id must be {cfg.PROBAND_ID}")
        try:
            e = float(r["epcr"])
            if not 0 < e <= 1:
                problems.append(f"{tag}: epcr must be in (0, 1]")
            if e > prev:
                problems.append(f"{tag}: rows are not ordered by decreasing epcr")
            prev = e
        except ValueError:
            problems.append(f"{tag}: epcr is not a number")
        if r.get("finding_type", "") not in ("", "primary", "secondary"):
            problems.append(f"{tag}: finding_type must be primary or secondary")
        variants = []
        for k in ("1", "2"):
            c, p, ref, alt = (r[f"chrom_{k}"].strip(), r[f"pos_{k}"].strip(),
                              r[f"ref_{k}"].strip().upper(), r[f"alt_{k}"].strip().upper())
            if not any((c, p, ref, alt)):
                continue
            if not all((c, p, ref, alt)):
                problems.append(f"{tag}: variant {k} is incomplete")
                continue
            if not re.fullmatch(r"chr([1-9]|1[0-9]|2[0-2]|X|Y|M)", c):
                problems.append(f"{tag}: chrom_{k} '{c}' must look like chr15 (GRCh38, 'chr' prefix)")
            if not p.isdigit():
                problems.append(f"{tag}: pos_{k} is not an integer")
            if not re.fullmatch(r"[ACGT]+", ref) or not re.fullmatch(r"[ACGT]+", alt):
                problems.append(f"{tag}: ref_{k}/alt_{k} must contain only A, C, G, T")
            variants.append((c, p, ref, alt))
        if not variants:
            problems.append(f"{tag}: no variant")
        key = frozenset(variants)
        if key in seen:
            problems.append(f"{tag}: duplicate prediction")
        seen.add(key)
        if vcf_path:
            from cyvcf2 import VCF
            for c, p, ref, alt in variants:
                if not p.isdigit():
                    continue
                found = any(v.POS == int(p) and v.REF == ref and alt in v.ALT
                            for v in VCF(str(vcf_path))(f"{c.replace('chr', '')}:{p}-{p}"))
                if not found:
                    problems.append(f"{tag}: {c}:{p} {ref}>{alt} is not in the VCF")
    return problems


def stage12_submission(out=None, verify_against_vcf: bool = True) -> dict:
    sel = read_json(cfg.path("selected"))
    rows = build_rows(sel)
    out = out or cfg.path("submission")
    write_csv(rows, out)
    problems = validate_submission(out, cfg.path("vcf") if verify_against_vcf else None)
    return {"rows": len(rows), "path": str(out), "problems": problems, "valid": not problems}
