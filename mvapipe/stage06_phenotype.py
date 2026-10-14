"""Stage 6: specificity-weighted phenotype ranking, leave-one-symptom-out, permutation null,
and rule-based candidate selection."""
from __future__ import annotations

import numpy as np
import pandas as pd

from . import config as cfg
from .utils import write_json


def load_hpo(path=None) -> pd.DataFrame:
    hpo = pd.read_csv(path or cfg.path("hpo"), sep="\t")
    return hpo.drop_duplicates(["hpo_id", "gene_symbol"])


def information_content(hpo: pd.DataFrame, terms) -> tuple[dict, pd.Series, int]:
    """IC(t) = -log2(n_genes(t) / N)."""
    N = hpo["gene_symbol"].nunique()
    n_genes = hpo.groupby("hpo_id")["gene_symbol"].nunique()
    ic = {t: float(-np.log2(n_genes[t] / N)) for t in terms if t in n_genes.index}
    return ic, n_genes, N


def gene_term_matrix(hpo: pd.DataFrame, ic: dict) -> pd.DataFrame:
    return (hpo[hpo["hpo_id"].isin(ic)].assign(v=1)
            .pivot_table(index="gene_symbol", columns="hpo_id", values="v", fill_value=0)
            .reindex(columns=list(ic), fill_value=0))


def damage_base(t: pd.DataFrame) -> pd.DataFrame:
    base = t[t["alto_mod"] & t["rara_o_sin_dato"] & ~t["clinvar_benigno"] & t["gen"].notna()].copy()
    grave = base["efecto"].fillna("").apply(lambda e: bool(set(e.split("&")) & cfg.TRUNCATING))
    base["w"] = np.where(grave | base["clinvar_patog"], cfg.W_STRONG, cfg.W_OTHER)
    return base


def damage_scores(base: pd.DataFrame) -> pd.Series:
    out = {}
    for gen, g in base.groupby("gen"):
        hom = g.loc[g["gt"] == 3, "w"]
        het = g.loc[g["gt"] == 1, "w"].sort_values(ascending=False)
        s_hom = cfg.W_HOM_MULT * hom.max() if len(hom) else 0.0
        s_ch = het.iloc[0] + het.iloc[1] if len(het) >= 2 else 0.0
        s_one = cfg.W_SINGLE_ALLELE * het.iloc[0] if len(het) else 0.0
        out[gen] = max(s_hom, s_ch, s_one)
    return pd.Series(out, dtype=float)


def make_ranker(M: pd.DataFrame, dano: pd.Series):
    def pheno_score(weights: dict) -> pd.Series:
        return (M[list(weights)] * pd.Series(weights)).sum(axis=1)

    def ranking(weights: dict) -> pd.Series:
        ps = pheno_score(weights)
        common = ps.index.intersection(dano.index)
        sc = (ps[common] * dano[common]).sort_values(ascending=False)
        return sc[sc > 0]

    return pheno_score, ranking


def permutation_null(hpo, ic: dict, n_genes: pd.Series, N: int, dano: pd.Series,
                     observed_gene_fit: float, observed_best: float,
                     gene: str, n_perm: int = cfg.PERM_N, seed: int = cfg.PERM_SEED) -> dict:
    """Replace every proband term by a random term of similar specificity and recompute.

    p_gene_fit: probability that `gene` alone reaches its observed phenotype fit.
    p_any_gene_combined: probability that ANY gene reaches the observed best combined score
    (corrects for picking the top gene after the fact)."""
    rng = np.random.default_rng(seed)
    ic_all = -np.log2(n_genes / N)
    term_genes = (hpo[hpo["gene_symbol"].isin(dano.index)]
                  .groupby("hpo_id")["gene_symbol"].apply(list).to_dict())
    twins = {h: ic_all[(ic_all - v).abs() < cfg.IC_TOL].index.difference(list(ic)).tolist()
             for h, v in ic.items()}
    gene_null, best_null = [], []
    for _ in range(n_perm):
        acc: dict[str, float] = {}
        for h in ic:
            if not twins[h]:
                continue
            term = rng.choice(twins[h])
            for g in term_genes.get(term, []):
                acc[g] = acc.get(g, 0.0) + float(ic_all[term])
        best_null.append(max((v * dano[g] for g, v in acc.items()), default=0.0))
        gene_null.append(acc.get(gene, 0.0))
    return {
        "n_permutations": n_perm, "seed": seed, "gene": gene,
        "observed_gene_phenotype_fit": observed_gene_fit,
        "observed_best_combined_score": observed_best,
        "p_gene_fit": float(np.mean(np.array(gene_null) >= observed_gene_fit)),
        "p_any_gene_combined": float(np.mean(np.array(best_null) >= observed_best)),
        "twin_terms_per_symptom": {h: len(v) for h, v in twins.items()},
    }


def select_candidates(t: pd.DataFrame, ranked: pd.Series) -> dict:
    """Primary: top-ranked gene with >= 2 heterozygous candidate alleles (two best form the pair).
    Secondary: truncating ClinVar-pathogenic variants outside that gene, gnomAD missing or rare."""
    base = damage_base(t)
    het = base[base["gt"] == 1]
    primary = None
    for gene in ranked.index:
        g = het[het["gen"] == gene].sort_values(["w", "pos"], ascending=[False, True])
        if len(g) >= 2:
            top = g.head(2).sort_values("pos")
            primary = {"gene": gene,
                       "variants": [{"chrom": str(r.chrom), "pos": int(r.pos), "ref": r.ref,
                                     "alt": r.alt, "efecto": r.efecto, "hgvs_p": r.hgvs_p}
                                    for r in top.itertuples()]}
            break
    secondary = []
    if primary:
        sec = t[t["clinvar_patog"] & ~t["clinvar_benigno"] & (t["gen"] != primary["gene"])
                & t["efecto"].fillna("").apply(lambda e: bool(set(e.split("&")) & cfg.TRUNCATING))
                & (t["gnomad_af"].isna() | (t["gnomad_af"] < cfg.RARE_AF))]
        for r in sec.sort_values("pos").itertuples():
            secondary.append({"gene": r.gen, "chrom": str(r.chrom), "pos": int(r.pos),
                              "ref": r.ref, "alt": r.alt, "efecto": r.efecto, "hgvs_p": r.hgvs_p,
                              "gt": int(r.gt)})
    return {"primary": primary, "secondary": secondary}


def stage06_phenotype(annotated_path=None, hpo_path=None, n_perm: int = cfg.PERM_N) -> dict:
    t = pd.read_parquet(annotated_path or cfg.path("annotated"))
    hpo = load_hpo(hpo_path)
    terms = cfg.load_hpo_terms()
    ic, n_genes, N = information_content(hpo, list(terms))
    M = gene_term_matrix(hpo, ic)
    base = damage_base(t)
    dano = damage_scores(base)
    pheno_score, ranking = make_ranker(M, dano)

    ranked = ranking(ic)
    ranked.rename("combined_score").rename_axis("gene").to_csv(cfg.path("ranking"))

    top_gene = ranked.index[0] if len(ranked) else None
    root = n_genes.get("HP:0000118")

    loo_rows = []
    for h in ic:
        rr = ranking({k: v for k, v in ic.items() if k != h})
        pos = int(list(rr.index).index(top_gene) + 1) if top_gene in rr.index else None
        loo_rows.append({"removed": terms[h], "hpo_id": h,
                         "rank_of_top_gene": pos, "n_scored_genes": len(rr),
                         "top_gene_without_this_symptom": rr.index[0] if len(rr) else None,
                         "score_of_top_gene": float(rr.get(top_gene, np.nan))})
    pd.DataFrame(loo_rows).to_csv(cfg.path("loo"), index=False)

    perm = {}
    if top_gene is not None:
        perm = permutation_null(hpo, ic, n_genes, N, dano,
                                observed_gene_fit=float(pheno_score(ic)[top_gene]),
                                observed_best=float(ranked.iloc[0]),
                                gene=top_gene, n_perm=n_perm)
    write_json(perm, cfg.path("perm"))

    selected = select_candidates(t, ranked)
    write_json(selected, cfg.path("selected"))

    summary = {
        "hpo_genes_N": int(N), "genes_under_root_term": None if root is None else int(root),
        "ic_weights": {terms[h]: round(v, 2) for h, v in ic.items()},
        "genes_scored": int(len(ranked)), "top_gene": top_gene,
        "top_score": float(ranked.iloc[0]) if len(ranked) else None,
        "runner_up": ranked.index[1] if len(ranked) > 1 else None,
        "runner_up_score": float(ranked.iloc[1]) if len(ranked) > 1 else None,
        "top_gene_rank_in_all_leave_one_out": all(r["rank_of_top_gene"] == 1 for r in loo_rows),
        "p_gene_fit": perm.get("p_gene_fit"), "p_any_gene_combined": perm.get("p_any_gene_combined"),
        "primary": selected["primary"], "secondary_n": len(selected["secondary"]),
    }
    write_json(summary, cfg.workdir() / "phenotype_summary.json")
    return summary
