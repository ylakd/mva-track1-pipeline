"""Stages 13-14: tests of the method itself.

Stage 13  Multiverse (sensitivity) analysis. The damage weights of the main ranking were fixed after the
          candidate was known. Here the ranking is recomputed over a grid of reasonable analyst choices
          (damage weights, rarity threshold, treatment of variants without a population frequency,
          term-weighting scheme) and we report how often, and by what margin, the top gene stays first.

Stage 14  In-silico recovery benchmark. A compound-heterozygous pair is planted in a randomly chosen gene
          of the HPO resource, a phenotype is simulated from that gene's own annotations (plus noise terms),
          the proband's real background of damaging variants is kept, and we ask at which rank the method
          returns the planted gene. Compared with two baselines (damage only, phenotype only) and with a
          no-signal control. The phenotype is generated from the same resource used for ranking, so the
          recovery rate is an upper bound on real-world performance.

Neither stage reads sequence data; both use only the annotated variant table and the HPO table, and write
only aggregate files that contain no variant coordinates.
"""
from __future__ import annotations

import itertools
import math
import time

import numpy as np
import pandas as pd

from . import config as cfg
from .stage06_phenotype import (damage_scores, damage_base, gene_term_matrix, information_content,
                                load_hpo)
from .utils import read_json, write_json

# ---------------------------------------------------------------- stage 13

GRID = {
    "w_other": [0.25, 0.5, 0.75],
    "w_hom_mult": [1.0, 1.5, 2.0, 3.0],
    "w_single": [0.0, 0.1, 0.25, 0.5],
    "rare_af": [1e-4, 1e-3, 1e-2],
    "no_value": ["keep", "drop"],
    "term_weight": ["ic", "uniform", "sqrt_ic", "ic_squared"],
}


def damage_with(t: pd.DataFrame, w_other: float, w_hom_mult: float, w_single: float,
                rare_af: float, no_value: str) -> pd.Series:
    """Damage score per gene under one parameter setting (same logic as stage 6, parametrised)."""
    af = t["gnomad_af"]
    rare = (af < rare_af) | (af.isna() if no_value == "keep" else False)
    base = t[t["alto_mod"] & rare & ~t["clinvar_benigno"] & t["gen"].notna()].copy()
    if base.empty:
        return pd.Series(dtype=float)
    grave = base["efecto"].fillna("").apply(lambda e: bool(set(e.split("&")) & cfg.TRUNCATING))
    base["w"] = np.where(grave | base["clinvar_patog"], cfg.W_STRONG, w_other)
    out = {}
    for gen, g in base.groupby("gen"):
        hom = g.loc[g["gt"] == 3, "w"]
        het = g.loc[g["gt"] == 1, "w"].sort_values(ascending=False)
        s_hom = w_hom_mult * hom.max() if len(hom) else 0.0
        s_ch = het.iloc[0] + het.iloc[1] if len(het) >= 2 else 0.0
        s_one = w_single * het.iloc[0] if len(het) else 0.0
        out[gen] = max(s_hom, s_ch, s_one)
    return pd.Series(out, dtype=float)


def reweight(ic: dict, scheme: str) -> dict:
    f = {"ic": lambda v: v, "uniform": lambda v: 1.0, "sqrt_ic": math.sqrt,
         "ic_squared": lambda v: v * v}[scheme]
    return {h: f(v) for h, v in ic.items()}


def stage13_robustness(annotated_path=None, hpo_path=None, grid: dict | None = None) -> dict:
    grid = grid or GRID
    t = pd.read_parquet(annotated_path or cfg.path("annotated"))
    hpo = load_hpo(hpo_path)
    terms = cfg.load_hpo_terms()
    ic, _, _ = information_content(hpo, list(terms))
    M = gene_term_matrix(hpo, ic)
    top = read_json(cfg.workdir() / "phenotype_summary.json")["top_gene"]

    dano_cache: dict = {}
    rows = []
    keys = list(grid)
    for combo in itertools.product(*[grid[k] for k in keys]):
        p = dict(zip(keys, combo))
        dkey = (p["w_other"], p["w_hom_mult"], p["w_single"], p["rare_af"], p["no_value"])
        if dkey not in dano_cache:
            dano_cache[dkey] = damage_with(t, *dkey)
        dano = dano_cache[dkey]
        w = reweight(ic, p["term_weight"])
        ps = (M[list(w)] * pd.Series(w)).sum(axis=1)
        common = ps.index.intersection(dano.index)
        sc = (ps[common] * dano[common]).sort_values(ascending=False)
        sc = sc[sc > 0]
        if top in sc.index:
            rank = int(list(sc.index).index(top) + 1)
            others = sc.drop(top)
            margin = float(sc[top] / others.iloc[0]) if len(others) and others.iloc[0] > 0 else float("inf")
        else:
            rank, margin = None, float("nan")
        rows.append({**p, "rank_of_top_gene": rank, "margin_vs_runner_up": margin,
                     "top_gene_scored": rank is not None})
    df = pd.DataFrame(rows)
    df.to_csv(cfg.workdir() / "robustness_grid.csv", index=False)

    ranks = df["rank_of_top_gene"].dropna()
    margins = df.loc[df["rank_of_top_gene"] == 1, "margin_vs_runner_up"].replace(np.inf, np.nan).dropna()
    by_factor = {}
    for k in keys:
        by_factor[k] = {str(v): float((df.loc[df[k] == v, "rank_of_top_gene"] == 1).mean())
                        for v in grid[k]}
    summary = {
        "top_gene_of_default_analysis": top,
        "n_settings": int(len(df)),
        "fraction_rank1": float((df["rank_of_top_gene"] == 1).mean()),
        "fraction_top3": float((df["rank_of_top_gene"] <= 3).mean()),
        "fraction_not_scored": float(df["rank_of_top_gene"].isna().mean()),
        "worst_rank": None if ranks.empty else int(ranks.max()),
        "median_margin_when_rank1": None if margins.empty else float(margins.median()),
        "min_margin_when_rank1": None if margins.empty else float(margins.min()),
        "fraction_rank1_by_factor_level": by_factor,
        "note": "A multiverse analysis over analyst choices; it addresses the non-blind choice of weights "
                "but not the dependence on the HPO resource itself.",
    }
    write_json(summary, cfg.workdir() / "robustness_summary.json")
    return summary


# ---------------------------------------------------------------- stage 14

def wilson(k: int, n: int, z: float = 1.96) -> tuple[float, float]:
    if n == 0:
        return (float("nan"), float("nan"))
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return (max(0.0, c - h), min(1.0, c + h))


def _rank(score: dict, gene: str) -> float:
    """Pessimistic rank: 1 + number of OTHER genes with a score >= the planted gene's."""
    s = score.get(gene, 0.0)
    if s <= 0:
        return float("inf")
    return 1.0 + sum(1 for g, v in score.items() if g != gene and v >= s)


def stage14_benchmark(annotated_path=None, hpo_path=None, n_targets: int = cfg.BENCH_N,
                      seed: int = cfg.BENCH_SEED, n_terms: int = 8,
                      signal_levels=(0, 2, 3, 4), planted_damage: float = cfg.W_STRONG + cfg.W_OTHER) -> dict:
    t0 = time.time()
    rng = np.random.default_rng(seed)
    t = pd.read_parquet(annotated_path or cfg.path("annotated"))
    hpo = load_hpo(hpo_path)
    N = hpo["gene_symbol"].nunique()
    n_genes = hpo.groupby("hpo_id")["gene_symbol"].nunique()
    ic_all = (-np.log2(n_genes / N)).to_dict()
    root = "HP:0000118"
    informative = [h for h, v in ic_all.items() if h != root and v >= 1.0]   # exclude near-root terms

    # background of damaging variants: the proband's real table, without the real candidate gene
    top = None
    f = cfg.workdir() / "phenotype_summary.json"
    if f.exists():
        top = read_json(f).get("top_gene")
    dano = damage_scores(damage_base(t))
    if top in dano.index:
        dano = dano.drop(top)
    background = dano.to_dict()
    bg_genes = set(background)

    term_genes = hpo.groupby("hpo_id")["gene_symbol"].apply(set).to_dict()
    gene_terms = hpo.groupby("gene_symbol")["hpo_id"].apply(
        lambda s: [h for h in set(s) if h in ic_all and h != root and ic_all[h] >= 1.0]).to_dict()
    pool = [g for g, ts in gene_terms.items() if len(ts) >= max(signal_levels) and g != top]
    targets = list(rng.choice(pool, size=min(n_targets, len(pool)), replace=False))

    rows = []
    for g in targets:
        own = gene_terms[g]
        for m in signal_levels:
            signal = list(rng.choice(own, size=min(m, len(own)), replace=False)) if m else []
            noise = list(rng.choice(informative, size=n_terms - len(signal), replace=False))
            patient = signal + noise
            w = {h: ic_all[h] for h in patient}
            universe = bg_genes | {g}
            ph = {}
            for h in patient:
                for gg in term_genes.get(h, ()):
                    if gg in universe:
                        ph[gg] = ph.get(gg, 0.0) + w[h]
            dmg = dict(background); dmg[g] = planted_damage
            combined = {gg: ph.get(gg, 0.0) * dmg[gg] for gg in universe}
            rows.append({
                "planted_gene": g, "signal_terms": m,
                "rank_combined": _rank(combined, g),
                "rank_phenotype_only": _rank({gg: ph.get(gg, 0.0) for gg in universe}, g),
                "rank_damage_only": _rank({gg: dmg[gg] for gg in universe}, g),
            })
    df = pd.DataFrame(rows)
    df.to_csv(cfg.workdir() / "benchmark_replicates.csv", index=False)

    out = {"n_targets": len(targets), "n_terms_per_patient": n_terms, "seed": seed,
           "background_genes_with_damage": len(bg_genes), "planted_damage": planted_damage,
           "by_signal_terms": {}}
    for m, g in df.groupby("signal_terms"):
        res = {"n": int(len(g))}
        for col, name in (("rank_combined", "combined"), ("rank_phenotype_only", "phenotype_only"),
                          ("rank_damage_only", "damage_only")):
            r = g[col].replace(np.inf, np.nan)
            res[name] = {
                "top1": float((g[col] <= 1).mean()), "top1_ci95": wilson(int((g[col] <= 1).sum()), len(g)),
                "top3": float((g[col] <= 3).mean()), "top10": float((g[col] <= 10).mean()),
                "median_rank": None if r.dropna().empty else float(r.median()),
                "not_returned": float(r.isna().mean()),
            }
        out["by_signal_terms"][str(int(m))] = res
    out["seconds"] = round(time.time() - t0, 1)
    out["seconds_per_simulated_patient"] = round((time.time() - t0) / max(1, len(df)), 4)
    out["caveats"] = [
        "The phenotype is simulated from the same HPO resource used for ranking: the recovery rate is an "
        "upper bound for real patients, whose annotations are incomplete.",
        "signal_terms = 0 is a no-signal control: it estimates the rate at which the planted gene would be "
        "first by chance.",
        "Ranks are pessimistic (ties count against the planted gene).",
        "Only compound-heterozygous truncating + missense genotypes are simulated.",
    ]
    write_json(out, cfg.workdir() / "benchmark_summary.json")
    return out
