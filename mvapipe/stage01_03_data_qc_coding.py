"""Stages 1-3: data download (with provenance), VCF quality control, coding-region restriction."""
from __future__ import annotations

import datetime as dt
import re
from collections import Counter

import numpy as np
import pandas as pd
import requests

from . import config as cfg
from .utils import (fmt, gtf_iter, hf_token, in_intervals, merge_intervals,
                    sha256_file, write_json)


def _download(url: str, dest) -> None:
    r = requests.get(url, stream=True, timeout=300)
    r.raise_for_status()
    with open(dest, "wb") as f:
        for chunk in r.iter_content(1 << 20):
            f.write(chunk)


def extract_hpo_terms(docx_path) -> dict[str, str]:
    """Read HPO identifiers and labels from the challenge phenotype document."""
    from docx import Document
    terms: dict[str, str] = {}
    doc = Document(docx_path)
    for table in doc.tables:
        for row in table.rows:
            cells = [c.text.strip() for c in row.cells]
            for i, c in enumerate(cells):
                if re.fullmatch(r"HP:\d{7}", c):
                    terms[c] = cells[i - 1] if i > 0 else ""
    return terms


def stage01_data(force: bool = False) -> dict:
    from huggingface_hub import hf_hub_download

    token = hf_token()
    wd = cfg.workdir()
    for name in (cfg.VCF_NAME, cfg.TBI_NAME, cfg.PHENO_NAME):
        dest = wd / "data" / name
        if dest.exists() and not force:
            continue
        hf_hub_download(cfg.REPO_ID, filename=name, repo_type="dataset",
                        local_dir=str(wd / "data"), token=token)
    for key, url in (("gtf", cfg.GTF_URL), ("hpo", cfg.HPO_URL)):
        dest = cfg.path(key)
        if not dest.exists() or force:
            _download(url, dest)

    prov = {"date_utc": dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"),
            "gtf_url": cfg.GTF_URL, "hpo_url": cfg.HPO_URL, "files": {}}
    for key in ("tbi", "pheno", "gtf", "hpo"):
        p = cfg.path(key)
        prov["files"][p.name] = {"bytes": p.stat().st_size, "sha256": sha256_file(p)}
    vcf = cfg.path("vcf")
    prov["files"][vcf.name] = {"bytes": vcf.stat().st_size}
    try:
        found = extract_hpo_terms(cfg.path("pheno"))
        write_json(found, cfg.path("hpo_terms"))          # private: stays in the working directory
        cfg.HPO_TERMS.clear(); cfg.HPO_TERMS.update(found)
        prov["n_hpo_terms_in_document"] = len(found)
    except Exception as e:
        prov["hpo_terms_in_document"] = f"could not parse: {e}"
    write_json(prov, cfg.path("provenance"))
    return prov


def stage02_qc(vcf_path=None, out=None, limit: int | None = None) -> dict:
    """Keep PASS, non-reference, DP >= 10, GQ >= 30 calls; record every exclusion count."""
    from cyvcf2 import VCF

    vcf_path = vcf_path or cfg.path("vcf")
    out = out or cfg.path("qc")
    counts: Counter = Counter()
    rows = []
    for v in VCF(str(vcf_path)):
        counts["total_records"] += 1
        if v.FILTER is not None:
            counts["excluded_filter_not_pass"] += 1
        else:
            gt = int(v.gt_types[0])
            if gt in (0, 2):
                counts["excluded_reference_or_unknown_gt"] += 1
            else:
                dp, gq = int(v.gt_depths[0]), int(v.gt_quals[0])
                if dp < cfg.QC_MIN_DP or gq < cfg.QC_MIN_GQ:
                    counts["excluded_low_dp_or_gq"] += 1
                else:
                    ad = fmt(v, "AD")
                    ad_ref, ad_alt = (int(ad[0][0]), int(ad[0][1])) if ad is not None else (-1, -1)
                    rows.append((v.CHROM, v.POS, v.REF, ",".join(v.ALT), v.QUAL, gt, dp, gq,
                                 bool(v.INFO.get("DB", False)), ad_ref, ad_alt))
                    counts["kept"] += 1
        if limit and counts["total_records"] >= limit:
            break

    df = pd.DataFrame(rows, columns=["chrom", "pos", "ref", "alt", "qual", "gt", "dp", "gq",
                                     "in_dbsnp", "ad_ref", "ad_alt"])
    df["chrom"] = df["chrom"].astype(str)
    df.to_parquet(out)
    stats = dict(counts)
    stats["het"] = int((df["gt"] == 1).sum())
    stats["hom_alt"] = int((df["gt"] == 3).sum())
    stats["in_dbsnp"] = int(df["in_dbsnp"].sum())
    stats["not_in_dbsnp"] = int((~df["in_dbsnp"]).sum())
    write_json(stats, cfg.path("qc_stats"))
    return stats


def cds_intervals(gtf_path, pad: int = cfg.CDS_PAD) -> dict[str, tuple[np.ndarray, np.ndarray]]:
    raw: dict[str, list[tuple[int, int]]] = {}
    for chrom, _feat, s, e, _strand, _attrs in gtf_iter(gtf_path, {"CDS"}):
        raw.setdefault(chrom, []).append((s - pad, e + pad))
    return {ch: merge_intervals(iv) for ch, iv in raw.items()}


def stage03_coding(qc_path=None, gtf_path=None, out=None) -> dict:
    qc_path = qc_path or cfg.path("qc")
    gtf_path = gtf_path or cfg.path("gtf")
    out = out or cfg.path("coding")
    df = pd.read_parquet(qc_path)
    merged = cds_intervals(gtf_path)

    mask = np.zeros(len(df), dtype=bool)
    pos_all = df["pos"].to_numpy()
    for chrom, idx in df.groupby("chrom").indices.items():
        if chrom in merged:
            s, e = merged[chrom]
            mask[idx] = in_intervals(s, e, pos_all[idx])

    coding = df[mask].reset_index(drop=True)
    coding.to_parquet(out)
    stats = {"input_variants": int(len(df)), "coding_variants": int(len(coding)),
             "coding_het": int((coding["gt"] == 1).sum()),
             "coding_hom_alt": int((coding["gt"] == 3).sum())}
    write_json(stats, cfg.workdir() / "coding_stats.json")
    return stats
