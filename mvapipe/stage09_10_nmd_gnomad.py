"""Stages 9-10 (supplementary): nonsense-mediated decay position check and gnomAD v4 lookup."""
from __future__ import annotations

import requests

from . import config as cfg
from .utils import gtf_iter, parse_attrs, read_json, write_json

NMD_RULE_NT = 50


def pick_transcript(gtf_path, gene: str) -> dict | None:
    """MANE Select; else Ensembl canonical; else the transcript with most exon bases."""
    tx: dict[str, dict] = {}
    for chrom, feat, s, e, strand, attrs in gtf_iter(gtf_path, {"transcript", "exon"}):
        if f'gene_name "{gene}"' not in attrs:
            continue
        a = parse_attrs(attrs)
        tid = a.get("transcript_id")
        d = tx.setdefault(tid, {"id": tid, "chrom": chrom, "strand": strand, "exons": [], "tags": []})
        if feat == "transcript":
            tags = a.get("tag", [])
            d["tags"] = tags if isinstance(tags, list) else [tags]
        else:
            d["exons"].append((s, e))
    if not tx:
        return None
    for want in ("MANE_Select", "Ensembl_canonical"):
        hit = [d for d in tx.values() if want in d["tags"]]
        if hit:
            hit[0]["selected_by"] = want
            return hit[0]
    best = max(tx.values(), key=lambda d: sum(e - s + 1 for s, e in d["exons"]))
    best["selected_by"] = "longest_exon_total"
    return best


def nmd_position(tx: dict, pos: int) -> dict:
    """Exon index and distance from `pos` to the last exon-exon junction, in spliced transcript space."""
    ex = sorted(tx["exons"], key=lambda x: x[0], reverse=(tx["strand"] == "-"))
    n = len(ex)
    offset, running, hit_idx = None, 0, None
    for i, (s, e) in enumerate(ex):
        if s <= pos <= e:
            offset = running + ((pos - s) if tx["strand"] == "+" else (e - pos))
            hit_idx = i + 1
        running += e - s + 1
    if offset is None:
        return {"in_exon": False}
    last_len = ex[-1][1] - ex[-1][0] + 1
    last_junction = running - last_len
    dist = last_junction - offset
    in_last = hit_idx == n
    return {"in_exon": True, "exon_number": hit_idx, "total_exons": n,
            "transcript_offset_nt": offset, "mrna_length_nt": running,
            "distance_to_last_junction_nt": dist, "in_last_exon": in_last,
            "predicted_nmd": bool((not in_last) and dist > NMD_RULE_NT),
            "interpretation": ("last exon: expected to escape NMD" if in_last else
                               f"{dist} nt upstream of the last junction: "
                               + ("NMD expected (>50 nt rule)" if dist > NMD_RULE_NT
                                  else "within 50 nt of the last junction: may escape NMD"))}


def stage09_nmd(gtf_path=None) -> dict:
    sel = read_json(cfg.path("selected"))
    gtf_path = gtf_path or cfg.path("gtf")
    results, targets = [], []
    if sel.get("primary"):
        targets += [dict(v, gene=sel["primary"]["gene"]) for v in sel["primary"]["variants"]]
    targets += sel.get("secondary", [])
    for v in targets:
        if not set(str(v.get("efecto", "")).split("&")) & cfg.TRUNCATING:
            continue
        tx = pick_transcript(gtf_path, v["gene"])
        if tx is None:
            results.append({"variant": v, "error": "gene not found"})
            continue
        r = nmd_position(tx, int(v["pos"]))
        r.update({"variant": v, "transcript": tx["id"], "transcript_selected_by": tx["selected_by"],
                  "caveat": "approximate: uses the variant position as the premature stop position "
                            "(error of at most 2 nt); frameshift stops lie further downstream"})
        results.append(r)
    write_json({"nmd_rule_nt": NMD_RULE_NT, "results": results}, cfg.path("nmd"))
    return {"n_checked": len(results)}


GNOMAD_URL = "https://gnomad.broadinstitute.org/api"
GNOMAD_QUERY = """
query($id: String!, $ds: DatasetId!) {
  variant(variantId: $id, dataset: $ds) {
    variant_id
    exome { ac an homozygote_count }
    genome { ac an homozygote_count }
  }
}"""


def gnomad_lookup(chrom: str, pos: int, ref: str, alt: str, dataset: str = "gnomad_r4", session=None) -> dict:
    """Query the gnomAD GraphQL API. Never raises: returns a status instead."""
    http = session or requests
    vid = f"{str(chrom).replace('chr', '')}-{pos}-{ref}-{alt}"
    try:
        r = http.post(GNOMAD_URL, json={"query": GNOMAD_QUERY, "variables": {"id": vid, "ds": dataset}},
                      timeout=60)
    except Exception as e:
        return {"variant_id": vid, "status": "unavailable", "detail": f"{type(e).__name__}: {e}"}
    if r.status_code != 200:
        return {"variant_id": vid, "status": "unavailable", "detail": f"HTTP {r.status_code}"}
    j = r.json()
    v = (j.get("data") or {}).get("variant")
    if v is None:
        msgs = [e.get("message", "") for e in j.get("errors", [])]
        if any("not found" in m.lower() for m in msgs):
            return {"variant_id": vid, "status": "absent_from_gnomad_v4", "detail": msgs}
        return {"variant_id": vid, "status": "unavailable", "detail": msgs or "empty response"}
    out = {"variant_id": vid, "status": "present"}
    for k in ("exome", "genome"):
        d = v.get(k)
        out[k] = None if not d else {"ac": d["ac"], "an": d["an"],
                                     "af": (d["ac"] / d["an"]) if d["an"] else None,
                                     "homozygotes": d.get("homozygote_count")}
    return out


def stage10_gnomad(session=None) -> dict:
    sel = read_json(cfg.path("selected"))
    variants = []
    if sel.get("primary"):
        variants += sel["primary"]["variants"]
    variants += sel.get("secondary", [])
    res = [gnomad_lookup(v["chrom"], v["pos"], v["ref"], v["alt"], session=session) for v in variants]
    write_json({"dataset": "gnomad_r4", "results": res,
                "note": "'unavailable' means the query could not be completed, NOT that the variant is absent. "
                        "Check unavailable entries manually in the gnomAD browser."}, cfg.path("gnomad"))
    return {"n": len(res), "statuses": [r["status"] for r in res]}
