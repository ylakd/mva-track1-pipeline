"""Stages 4-5: batched annotation (MyVariant.info) and the pre-specified, phenotype-agnostic scan."""
from __future__ import annotations

import json
import time

import numpy as np
import pandas as pd
import requests

from . import config as cfg
from .utils import GENE_RE, as_list, flat_nums, gtf_iter, hgvs_genomic, maxnum, write_json


def add_hgvs(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    df["hgvs"] = [hgvs_genomic(c, p, r, a) for c, p, r, a in
                  zip(df["chrom"], df["pos"], df["ref"], df["alt"])]
    return df


def stage04_annotate(coding_path=None, force: bool = False, session=None) -> dict:
    """POST variants to MyVariant.info in batches; the raw JSON is cached on disk."""
    coding_path = coding_path or cfg.path("coding")
    raw_path = cfg.path("myvariant_raw")
    df = add_hgvs(pd.read_parquet(coding_path))
    df.to_parquet(cfg.workdir() / "variants_coding_hgvs.parquet")
    n_multi = int(df["hgvs"].isna().sum())

    if raw_path.exists() and not force:
        with open(raw_path) as fh:
            raw = json.load(fh)
    else:
        http = session or requests
        ids = df["hgvs"].dropna().unique().tolist()
        raw = []
        for i in range(0, len(ids), cfg.MYVARIANT_BATCH):
            batch = ids[i:i + cfg.MYVARIANT_BATCH]
            for attempt in range(3):
                r = http.post("https://myvariant.info/v1/variant",
                              data={"ids": ",".join(batch), "assembly": "hg38",
                                    "fields": cfg.MYVARIANT_FIELDS}, timeout=120)
                if r.status_code == 200:
                    raw.extend(r.json())
                    break
                time.sleep(3 * (attempt + 1))
            else:
                raise RuntimeError(f"MyVariant batch starting at {i} failed after 3 attempts")
        with open(raw_path, "w") as fh:
            json.dump(raw, fh)
    stats = {"coding_variants": int(len(df)), "multiallelic_not_annotated": n_multi,
             "responses": len(raw), "notfound": sum(1 for x in raw if x.get("notfound"))}
    write_json(stats, cfg.workdir() / "annotation_stats.json")
    return stats


def parse_record(x: dict) -> dict:
    """One MyVariant record -> flat dict (worst effect over transcripts; max AF; ClinVar flags)."""
    out = {"hgvs": x.get("query")}
    best = None
    for a in as_list((x.get("snpeff") or {}).get("ann")):
        eff = str(a.get("effect", ""))
        parts = set(eff.split("&"))
        rank = 2 if parts & cfg.IMPACT_HIGH else 1 if parts & cfg.IMPACT_MODERATE else 0
        if best is None or rank > best[0]:
            best = (rank, eff, a.get("genename") or a.get("gene_name"), a.get("hgvs_p"))
    out["impacto"], out["efecto"], out["gen_api"], out["hgvs_p"] = best if best else (np.nan,) * 4

    afs: list[float] = []
    for k in ("gnomad_genome", "gnomad_exome"):
        afs += flat_nums(((x.get(k) or {}).get("af") or {}).get("af"))
    out["gnomad_af"] = max(afs) if afs else np.nan

    cs = [str(r.get("clinical_significance", ""))
          for r in as_list((x.get("clinvar") or {}).get("rcv"))]
    out["clinvar_patog"] = any(("athogenic" in c) and ("Conflicting" not in c) for c in cs)
    out["clinvar_benigno"] = any("enign" in c for c in cs)
    dbn = x.get("dbnsfp") or {}
    cadd = dbn.get("cadd")
    out["cadd"] = maxnum(cadd.get("phred")) if isinstance(cadd, dict) else float("nan")
    out["revel"] = maxnum(dbn.get("revel"))
    return out


def assign_genes_local(t: pd.DataFrame, gtf_path) -> pd.Series:
    """Assign a gene symbol from CDS intervals in the GTF (no web service needed)."""
    rows = []
    for chrom, _f, s, e, _st, attrs in gtf_iter(gtf_path, {"CDS"}):
        m = GENE_RE.search(attrs)
        if m:
            rows.append((chrom, s - cfg.CDS_PAD, e + cfg.CDS_PAD, m.group(1)))
    cds = pd.DataFrame(rows, columns=["chrom", "s", "e", "gene"]).drop_duplicates()

    gene = np.full(len(t), None, dtype=object)
    assigned = np.zeros(len(t), dtype=bool)
    chrom_arr = t["chrom"].to_numpy()
    pos_arr = t["pos"].to_numpy()
    for ch, sub in cds.groupby("chrom"):
        sub = sub.sort_values("s")
        s, e, g = sub["s"].to_numpy(), sub["e"].to_numpy(), sub["gene"].to_numpy()
        ix = np.where(chrom_arr == ch)[0]
        if len(ix) == 0:
            continue
        pos = pos_arr[ix]
        k = np.searchsorted(s, pos, side="right") - 1
        for back in range(4):
            kk = k - back
            ok = (kk >= 0) & (pos <= e[np.clip(kk, 0, None)]) & ~assigned[ix]
            for j in np.where(ok)[0]:
                gene[ix[j]] = g[kk[j]]
                assigned[ix[j]] = True
    return pd.Series(gene, index=t.index)


def stage05_scan(hgvs_path=None, raw_path=None, gtf_path=None) -> dict:
    hgvs_path = hgvs_path or cfg.workdir() / "variants_coding_hgvs.parquet"
    raw_path = raw_path or cfg.path("myvariant_raw")
    gtf_path = gtf_path or cfg.path("gtf")

    exo = pd.read_parquet(hgvs_path)
    with open(raw_path) as fh:
        raw = json.load(fh)
    ann = pd.DataFrame([parse_record(x) for x in raw if not x.get("notfound")]) \
        .drop_duplicates("hgvs")
    t = exo.merge(ann, on="hgvs", how="left")

    t["gen_local"] = assign_genes_local(t, gtf_path)
    t["gen"] = t["gen_local"].fillna(t["gen_api"])

    t["rara"] = t["gnomad_af"].lt(cfg.RARE_AF)
    t["sin_dato_gnomad"] = t["gnomad_af"].isna()
    t["rara_o_sin_dato"] = t["rara"] | t["sin_dato_gnomad"]
    t["alto_mod"] = t["impacto"].isin([1, 2])
    for c in ("clinvar_patog", "clinvar_benigno", "rara", "alto_mod", "rara_o_sin_dato"):
        t[c] = t[c].astype("boolean").fillna(False).astype(bool)
    t.to_parquet(cfg.path("annotated"))

    cols = ["gen", "chrom", "pos", "ref", "alt", "gt", "efecto", "hgvs_p", "gnomad_af",
            "clinvar_patog", "cadd", "revel", "dp", "gq", "ad_ref", "ad_alt"]
    base = t[t["alto_mod"] & t["rara_o_sin_dato"] & ~t["clinvar_benigno"]]

    A = base[base["gt"] == 3].sort_values("gen")
    het = base[base["gt"] == 1]
    n = het.groupby("gen").size()
    B = het[het["gen"].isin(n[n >= 2].index)].sort_values(["gen", "pos"])
    efecto_grave = t["efecto"].fillna("").apply(lambda e: bool(set(e.split("&")) & cfg.TRUNCATING))
    genes_trunc = sorted(set(B[efecto_grave.reindex(B.index)]["gen"]))
    cp = t[t["clinvar_patog"] & ~t["clinvar_benigno"]].sort_values("gnomad_af")

    A[cols].to_csv(cfg.path("scan_A"), index=False)
    B[cols].to_csv(cfg.path("scan_B"), index=False)
    cp[cols].to_csv(cfg.path("scan_clinvar"), index=False)

    summary = {
        "annotated_with_effect": int(t["efecto"].notna().sum()),
        "with_gnomad_value": int((~t["sin_dato_gnomad"]).sum()),
        "without_gnomad_value": int(t["sin_dato_gnomad"].sum()),
        "gene_assigned_api": int(t["gen_api"].notna().sum()),
        "gene_assigned_local": int(t["gen_local"].notna().sum()),
        "scan_A_rows": int(len(A)),
        "scan_B_genes": int(B["gen"].nunique()),
        "scan_B_rows": int(len(B)),
        "scan_B_genes_with_truncating": genes_trunc,
        "clinvar_pathogenic_coding": int(len(cp)),
    }
    write_json(summary, cfg.path("scan_summary"))
    return summary
