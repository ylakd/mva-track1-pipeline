"""Stage 11 (supplementary): allele-fraction (B-allele frequency) and depth screen for chromosome imbalance.

Heterozygous-site allele fractions are a second, partly independent signal besides depth: in cells with an
extra or missing chromosome they drift away from 0.5. Both are summarised per chromosome and per 5 Mb
window as robust z-scores against the genome-wide distribution, and the smallest UNIFORM mosaic fraction
the depth signal could detect is estimated.

LIMITATION: in mosaic variegated aneuploidy the affected chromosomes differ from cell to cell, so a
single chromosome is imbalanced in few cells. A null result bounds what a bulk test could see; it does
not exclude the condition.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from . import config as cfg
from .utils import write_json

AUTOSOMES = [str(i) for i in range(1, 23)]


def robust_z(x: pd.Series) -> pd.Series:
    med = x.median()
    mad = (x - med).abs().median() * 1.4826
    return (x - med) / mad if mad and mad > 0 else x * 0.0


def trimmed_mean(s: pd.Series, lo: float = 0.05, hi: float = 0.95) -> float:
    """Mean between the lo and hi quantiles; unlike the median of integer depths it is not quantised."""
    a, b = s.quantile([lo, hi])
    return float(s[(s >= a) & (s <= b)].mean())


def het_snv_table(qc: pd.DataFrame) -> pd.DataFrame:
    d = qc[(qc["gt"] == 1) & (qc["ref"].str.len() == 1) & (qc["alt"].str.len() == 1)
           & (~qc["alt"].str.contains(",")) & (qc["dp"] >= cfg.BAF_MIN_DP)].copy()
    tot = d["ad_ref"] + d["ad_alt"]
    d = d[(d["ad_ref"] >= 0) & (tot > 0)]
    d["alt_frac"] = d["ad_alt"] / (d["ad_ref"] + d["ad_alt"])
    d = d[d["dp"] <= 2 * d["dp"].median()]
    d["dev"] = (d["alt_frac"] - 0.5).abs()
    return d


def stage11_baf(qc_path=None, bin_size: int = cfg.BIN_SIZE) -> dict:
    qc = pd.read_parquet(qc_path or cfg.path("qc"))
    het = het_snv_table(qc)
    chroms = AUTOSOMES + ["X"]

    g = het[het["chrom"].isin(chroms)].groupby("chrom")
    chrom_tab = pd.DataFrame({
        "n_het_snv": g.size(),
        "dp_trimmed_mean": g["dp"].apply(trimmed_mean),
        "median_alt_frac": g["alt_frac"].median(),
        "mean_abs_dev_from_0.5": g["dev"].mean(),
        "frac_outside_0.4_0.6": g["alt_frac"].apply(lambda s: float(((s < 0.4) | (s > 0.6)).mean())),
    }).reindex(chroms)
    auto = chrom_tab.loc[AUTOSOMES]
    chrom_tab["depth_ratio"] = chrom_tab["dp_trimmed_mean"] / auto["dp_trimmed_mean"].median()
    chrom_tab["z_depth"] = robust_z(chrom_tab.loc[AUTOSOMES, "depth_ratio"]).reindex(chroms)
    chrom_tab["z_baf_dev"] = robust_z(chrom_tab.loc[AUTOSOMES, "mean_abs_dev_from_0.5"]).reindex(chroms)
    chrom_tab["flag"] = (chrom_tab["z_depth"].abs() > 3.5) | (chrom_tab["z_baf_dev"].abs() > 3.5)
    chrom_tab.to_csv(cfg.path("baf_chrom"))

    h = het[het["chrom"].isin(AUTOSOMES)].copy()
    h["bin"] = (h["pos"] // bin_size).astype(int)
    gb = h.groupby(["chrom", "bin"])
    bins = pd.DataFrame({"n_het_snv": gb.size(), "dp_trimmed_mean": gb["dp"].apply(trimmed_mean),
                         "mean_abs_dev": gb["dev"].mean()}).reset_index()
    bins = bins[bins["n_het_snv"] >= cfg.BAF_MIN_SITES_PER_BIN].copy()
    bins["depth_ratio"] = bins["dp_trimmed_mean"] / bins["dp_trimmed_mean"].median()
    bins["z_depth"] = robust_z(bins["depth_ratio"])
    bins["z_baf_dev"] = robust_z(bins["mean_abs_dev"])
    bins["start"] = bins["bin"] * bin_size
    bins["flag"] = (bins["z_depth"].abs() > 3.5) | (bins["z_baf_dev"].abs() > 3.5)
    bins.to_csv(cfg.path("baf_bins"), index=False)

    # uniform mosaic fraction f shifts depth ratio by ~ f/2 (trisomy +f/2, monosomy -f/2)
    ratios = auto["dp_trimmed_mean"] / auto["dp_trimmed_mean"].median()
    spread = max(float(ratios.std()), float((ratios - 1).abs().median() * 1.4826))
    summary = {
        "het_snv_used": int(len(het)), "bin_size_bp": bin_size, "bins_tested": int(len(bins)),
        "chromosomes_flagged": chrom_tab.index[chrom_tab["flag"].fillna(False)].tolist(),
        "bins_flagged": int(bins["flag"].sum()),
        "flagged_bins_detail": bins[bins["flag"]][["chrom", "start", "n_het_snv", "depth_ratio", "z_depth",
                                                   "z_baf_dev"]].round(3).to_dict("records")[:20],
        "autosome_depth_ratio_spread": round(spread, 4),
        "approx_min_detectable_uniform_mosaic_fraction_depth": round(min(1.0, 2 * 3.5 * spread), 4),
        "x_depth_ratio": round(float(chrom_tab.loc["X", "depth_ratio"]), 3),
        "x_frac_het_among_all_calls": None,
        "caveat": "Bulk sample; MVA aneuploidies vary cell to cell, so each chromosome is imbalanced in few "
                  "cells. Only variants passing QC contribute. The detection limit is approximate and "
                  "assumes a uniform mosaic fraction.",
    }
    xs = qc[qc["chrom"] == "X"]
    if len(xs):
        summary["x_frac_het_among_all_calls"] = round(float((xs["gt"] == 1).mean()), 3)
    write_json(summary, cfg.path("baf_summary"))
    return summary
