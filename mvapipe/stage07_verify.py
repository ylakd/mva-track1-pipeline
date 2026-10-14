"""Stage 7: technical and functional verification of the selected candidates."""
from __future__ import annotations

import re

import pandas as pd
import requests

from . import config as cfg
from .utils import as_list, fmt, hgvs_genomic, maxnum, read_json, write_json

MV_FIELDS = "clinvar,gnomad_genome.af,gnomad_exome.af,dbnsfp,snpeff.ann"


def vcf_call(vcf_path, chrom: str, pos: int, ref: str, alt: str) -> dict:
    """Read one call straight from the VCF (contig names without the 'chr' prefix)."""
    from cyvcf2 import VCF
    chrom = str(chrom).replace("chr", "")
    for v in VCF(str(vcf_path))(f"{chrom}:{pos}-{pos}"):
        if v.POS == pos and v.REF == ref and alt in v.ALT:
            ad = fmt(v, "AD")
            ad_ref, ad_alt = (int(ad[0][0]), int(ad[0][1])) if ad is not None else (None, None)
            tot = (ad_ref + ad_alt) if ad_ref is not None else None
            g = v.genotypes[0]
            return {"found": True, "chrom": chrom, "pos": pos, "ref": v.REF, "alt": list(v.ALT),
                    "filter": v.FILTER or "PASS", "genotype": g[:2], "phased": bool(g[2]),
                    "ad_ref": ad_ref, "ad_alt": ad_alt,
                    "alt_fraction": round(ad_alt / tot, 3) if tot else None,
                    "dp": int(v.gt_depths[0]), "gq": int(v.gt_quals[0]),
                    "has_PGT": fmt(v, "PGT") is not None, "has_PID": fmt(v, "PID") is not None}
    return {"found": False, "chrom": chrom, "pos": pos}


def intervening_hets(qc_df: pd.DataFrame, chrom: str, p1: int, p2: int) -> dict:
    lo, hi = sorted((p1, p2))
    w = qc_df[(qc_df["chrom"] == str(chrom)) & qc_df["pos"].between(lo, hi) & (qc_df["gt"] == 1)] \
        .sort_values("pos")
    gaps = w["pos"].diff().dropna().astype(int).tolist()
    return {"distance_bp": hi - lo, "heterozygous_calls_in_window": int(len(w)),
            "positions": w["pos"].astype(int).tolist(),
            "not_in_dbsnp": w.loc[~w["in_dbsnp"], "pos"].astype(int).tolist(),
            "gaps_bp": gaps, "max_gap_bp": max(gaps) if gaps else None,
            "short_read_phasing_feasible": bool(gaps) and max(gaps) < 300}


def myvariant_detail(chrom, pos, ref, alt, session=None) -> dict:
    http = session or requests
    h = hgvs_genomic(chrom, pos, ref, alt)
    try:
        r = http.get(f"https://myvariant.info/v1/variant/{h}",
                     params={"assembly": "hg38", "fields": MV_FIELDS}, timeout=60)
        if r.status_code != 200:
            return {"hgvs": h, "http": r.status_code}
        j = r.json()
    except Exception as e:
        return {"hgvs": h, "error": str(e)}
    cv = []
    for rcv in as_list((j.get("clinvar") or {}).get("rcv")):
        cond = rcv.get("conditions")
        names = [cond.get("name")] if isinstance(cond, dict) else [c.get("name") for c in as_list(cond)]
        cv.append({"significance": rcv.get("clinical_significance"), "conditions": names,
                   "review_status": rcv.get("review_status"),
                   "submitters": rcv.get("number_submitters")})
    db = j.get("dbnsfp") or {}
    cadd = db.get("cadd") if isinstance(db.get("cadd"), dict) else {}
    ann = as_list((j.get("snpeff") or {}).get("ann"))
    return {
        "hgvs": h, "clinvar": cv,
        "gnomad_genome_af": maxnum(((j.get("gnomad_genome") or {}).get("af") or {}).get("af")),
        "gnomad_exome_af": maxnum(((j.get("gnomad_exome") or {}).get("af") or {}).get("af")),
        "predictors": {
            "cadd_phred": maxnum(cadd.get("phred")), "revel": maxnum(db.get("revel")),
            "alphamissense": maxnum(db.get("alphamissense")),
            "polyphen2_hdiv": maxnum((db.get("polyphen2") or {}).get("hdiv")),
            "sift": maxnum(db.get("sift")), "metarnn": maxnum(db.get("metarnn")),
            "phylop_100way": maxnum((db.get("phylop") or {}).get("100way_vertebrate")),
            "phastcons_100way": maxnum((db.get("phastcons") or {}).get("100way_vertebrate"))},
        "uniprot_acc": [u.get("acc") for u in as_list(db.get("uniprot")) if isinstance(u, dict)],
        "snpeff_effects": [a.get("effect") for a in ann][:6],
    }


def uniprot_features(accession: str, residue: int, window: int = 15, session=None) -> dict:
    http = session or requests
    try:
        r = http.get(f"https://rest.uniprot.org/uniprotkb/{accession}.json", timeout=60)
        if r.status_code != 200:
            return {"accession": accession, "http": r.status_code}
        j = r.json()
    except Exception as e:
        return {"accession": accession, "error": str(e)}
    feats = []
    for ft in j.get("features", []):
        s, e = ft["location"]["start"].get("value"), ft["location"]["end"].get("value")
        if s is None or e is None:
            continue
        if ft["type"] in ("Domain", "Region", "Motif", "Active site", "Binding site",
                          "Natural variant", "Mutagenesis", "Modified residue") \
                and (s <= residue <= e or abs(s - residue) <= window or abs(e - residue) <= window):
            feats.append({"type": ft["type"], "start": s, "end": e,
                          "description": (ft.get("description") or "")[:140],
                          "contains_residue": s <= residue <= e})
    return {"accession": accession, "length": j["sequence"]["length"], "residue": residue,
            "features_near_residue": feats}


def residue_from_hgvs_p(h: str | None) -> int | None:
    m = re.search(r"p\.[A-Za-z]{3}(\d+)", h or "")
    return int(m.group(1)) if m else None


def stage07_verify(session=None) -> dict:
    sel = read_json(cfg.path("selected"))
    qc_df = pd.read_parquet(cfg.path("qc"))
    out: dict = {"primary": None, "secondary": []}
    prim = sel.get("primary")
    if prim:
        calls = [vcf_call(cfg.path("vcf"), v["chrom"], v["pos"], v["ref"], v["alt"])
                 for v in prim["variants"]]
        detail = [myvariant_detail(v["chrom"], v["pos"], v["ref"], v["alt"], session)
                  for v in prim["variants"]]
        v1, v2 = prim["variants"][0], prim["variants"][1]
        phase = intervening_hets(qc_df, v1["chrom"], v1["pos"], v2["pos"])
        phase["any_phase_tags"] = any(c.get("phased") or c.get("has_PGT") or c.get("has_PID")
                                      for c in calls if c.get("found"))
        uni = []
        for v, d in zip(prim["variants"], detail):
            res = residue_from_hgvs_p(v.get("hgvs_p"))
            acc = (d.get("uniprot_acc") or [None])[0]
            if res and acc:
                uni.append(uniprot_features(acc.split("-")[0], res, session=session))
        out["primary"] = {"gene": prim["gene"], "variants": prim["variants"], "vcf_calls": calls,
                          "annotation": detail, "pair_geometry_and_phase": phase, "uniprot": uni}
    for s in sel.get("secondary", []):
        out["secondary"].append({
            "variant": s,
            "vcf_call": vcf_call(cfg.path("vcf"), s["chrom"], s["pos"], s["ref"], s["alt"]),
            "annotation": myvariant_detail(s["chrom"], s["pos"], s["ref"], s["alt"], session)})
    write_json(out, cfg.path("verification"))
    return out
