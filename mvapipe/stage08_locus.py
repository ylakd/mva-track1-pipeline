"""Stage 8 (supplementary): unfiltered review of the primary candidate's whole locus.

All records in the gene body plus flanks are read from the VCF with NO quality filter and each is
classified by position (CDS, UTR exon, intron, flank). This asks whether anything outside coding
sequence, or anything our quality filters removed, could complete or contradict the coding pair.
"""
from __future__ import annotations

from collections import Counter

import numpy as np
import pandas as pd
import requests

from . import config as cfg
from .utils import (flat_nums, fmt, gtf_iter, hgvs_genomic, in_intervals, maxnum, merge_intervals,
                    read_json, write_json)


def gene_structure(gtf_path, gene: str) -> dict | None:
    span = None
    exons, cds = [], []
    for chrom, feat, s, e, strand, attrs in gtf_iter(gtf_path, {"gene", "exon", "CDS"}):
        if f'gene_name "{gene}"' not in attrs:
            continue
        if feat == "gene":
            span = (chrom, s, e, strand)
        elif feat == "exon":
            exons.append((s, e))
        else:
            cds.append((s, e))
    if span is None:
        return None
    return {"chrom": span[0], "start": span[1], "end": span[2], "strand": span[3],
            "exons": merge_intervals(exons), "cds": merge_intervals(cds)}


def classify(pos: int, g: dict) -> tuple[str, int | None]:
    p = np.array([pos])
    if in_intervals(*g["cds"], p)[0]:
        cls = "CDS"
    elif in_intervals(*g["exons"], p)[0]:
        cls = "UTR_exon"
    elif g["start"] <= pos <= g["end"]:
        cls = "intron"
    else:
        up = pos < g["start"] if g["strand"] == "+" else pos > g["end"]
        cls = "upstream_flank" if up else "downstream_flank"
    s, e = g["exons"]
    edges = np.concatenate([s, e]) if len(s) else np.array([])
    dist = int(np.min(np.abs(edges - pos))) if len(edges) else None
    return cls, dist


def read_locus(vcf_path, chrom: str, start: int, end: int) -> pd.DataFrame:
    from cyvcf2 import VCF
    rows = []
    for v in VCF(str(vcf_path))(f"{chrom}:{start}-{end}"):
        ad = fmt(v, "AD")
        rows.append({"chrom": str(v.CHROM), "pos": v.POS, "ref": v.REF, "alt": ",".join(v.ALT),
                     "qual": v.QUAL, "filter": v.FILTER or "PASS", "gt": int(v.gt_types[0]),
                     "dp": int(v.gt_depths[0]), "gq": int(v.gt_quals[0]),
                     "ad_ref": int(ad[0][0]) if ad is not None else -1,
                     "ad_alt": int(ad[0][1]) if ad is not None else -1,
                     "in_dbsnp": bool(v.INFO.get("DB", False))})
    return pd.DataFrame(rows)


def annotate_rare(df: pd.DataFrame, cap: int = 300, session=None) -> pd.DataFrame:
    """Add gnomAD AF and CADD for non-dbSNP records (best effort; failures are recorded)."""
    http = session or requests
    df = df.copy()
    df["gnomad_af"] = np.nan
    df["cadd"] = np.nan
    df["annotation_status"] = "not_requested"
    sub = df[(~df["in_dbsnp"]) & (df["gt"].isin([1, 3]))].head(cap)
    ids = {i: hgvs_genomic(r.chrom, r.pos, r.ref, r.alt) for i, r in sub.iterrows()}
    ids = {i: h for i, h in ids.items() if h}
    if not ids:
        return df
    try:
        r = http.post("https://myvariant.info/v1/variant",
                      data={"ids": ",".join(ids.values()), "assembly": "hg38",
                            "fields": "gnomad_genome.af.af,gnomad_exome.af.af,dbnsfp.cadd.phred,cadd.phred"},
                      timeout=120)
        r.raise_for_status()
        by_q = {x.get("query"): x for x in r.json()}
    except Exception as e:
        df.loc[list(ids), "annotation_status"] = f"failed: {type(e).__name__}"
        return df
    for i, h in ids.items():
        x = by_q.get(h, {})
        if x.get("notfound"):
            df.at[i, "annotation_status"] = "not_found"
            continue
        afs = []
        for k in ("gnomad_genome", "gnomad_exome"):
            afs += flat_nums(((x.get(k) or {}).get("af") or {}).get("af"))
        df.at[i, "gnomad_af"] = max(afs) if afs else np.nan
        cadd = (x.get("dbnsfp") or {}).get("cadd")
        df.at[i, "cadd"] = maxnum(cadd.get("phred")) if isinstance(cadd, dict) \
            else maxnum((x.get("cadd") or {}).get("phred"))
        df.at[i, "annotation_status"] = "annotated"
    return df


def stage08_locus(gene: str | None = None, vcf_path=None, gtf_path=None,
                  flank: int = cfg.LOCUS_FLANK, session=None, annotate: bool = True) -> dict:
    gene = gene or (read_json(cfg.path("selected")).get("primary") or {}).get("gene")
    if not gene:
        raise RuntimeError("No primary gene selected (run stage 6 first).")
    g = gene_structure(gtf_path or cfg.path("gtf"), gene)
    if g is None:
        raise RuntimeError(f"Gene {gene} not found in the GTF")
    df = read_locus(vcf_path or cfg.path("vcf"), g["chrom"], max(1, g["start"] - flank), g["end"] + flank)
    cls = [classify(int(p), g) for p in df["pos"]]
    df["region"] = [c[0] for c in cls]
    df["dist_to_exon_edge"] = [c[1] for c in cls]
    df["splice_region_pm8"] = (df["region"] == "intron") & (df["dist_to_exon_edge"].fillna(10**9) <= 8)
    if annotate:
        df = annotate_rare(df, session=session)
    else:
        df["gnomad_af"], df["cadd"], df["annotation_status"] = np.nan, np.nan, "not_requested"
    df.to_csv(cfg.path("locus"), index=False)

    nonref = df[df["gt"].isin([1, 3])]
    brief = ["pos", "ref", "alt", "gt", "filter", "dp", "gq", "in_dbsnp"]
    summary = {
        "gene": gene, "chrom": g["chrom"], "gene_start": g["start"], "gene_end": g["end"],
        "strand": g["strand"], "flank_bp": flank,
        "records_total": int(len(df)), "records_non_reference": int(len(nonref)),
        "by_region_non_reference": dict(Counter(nonref["region"])),
        "by_filter_non_reference": dict(Counter(nonref["filter"])),
        "non_reference_not_in_dbsnp": int((~nonref["in_dbsnp"]).sum()),
        "coding_non_reference": nonref[nonref["region"] == "CDS"][brief].to_dict("records"),
        "splice_region_non_reference": nonref[nonref["splice_region_pm8"]][brief].to_dict("records"),
        "non_pass_non_reference": nonref[nonref["filter"] != "PASS"][brief].to_dict("records"),
        "rare_noncoding_candidates": nonref[(nonref["region"] != "CDS") & (~nonref["in_dbsnp"])
                                            & (nonref["gnomad_af"].isna() | (nonref["gnomad_af"] < cfg.RARE_AF))]
        [["pos", "ref", "alt", "gt", "region", "dist_to_exon_edge", "dp", "gq", "gnomad_af", "cadd",
          "annotation_status"]].to_dict("records"),
    }
    write_json(summary, cfg.path("locus_summary"))
    return summary
