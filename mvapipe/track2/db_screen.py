"""Database screen for candidate drugs (run where Internet is available, e.g. Kaggle with Internet ON).

Queries three public resources for the gene sets of the mechanism (spindle checkpoint, interferon / cGAS-STING axis,
senescence, NAD+/sirtuin axis) and writes a table of drug-gene pairs plus a pathway-enrichment table:
  * DGIdb (GraphQL)            drug-gene interactions, approval flag
  * Open Targets (GraphQL)     clinical-stage drugs per target (field names have changed between releases; both tried)
  * Enrichr (REST)             pathway enrichment of the gene sets
Every request failure is recorded, never silently skipped. No patient data are used or sent: only gene symbols.
Usage:  python -m mvapipe.track2.db_screen --out OUTDIR
"""
from __future__ import annotations

import argparse
import json
import time
import urllib.parse
import urllib.request
from pathlib import Path

GENE_SETS = {
    "spindle_checkpoint": ["BUB1B", "BUB1", "BUB3", "MAD2L1", "MAD1L1", "CDC20", "TTK", "AURKB", "PPP2R5A", "PLK1", "KNL1"],
    "innate_immune_axis": ["CGAS", "STING1", "TBK1", "IRF3", "IFNAR1", "IFNAR2", "JAK1", "JAK2", "TYK2", "STAT1", "STAT2"],
    "senescence_sasp": ["CDKN2A", "CDKN1A", "TP53", "BCL2", "BCL2L1", "SRC", "IL6", "MMP3"],
    "nad_sirtuin": ["SIRT1", "SIRT2", "NAMPT", "NMNAT1", "NAPRT", "PARP1"],
}
UA = {"User-Agent": "mva-track2-screen/1.0", "Content-Type": "application/json"}


def _post(url: str, payload: dict, timeout: int = 60):
    req = urllib.request.Request(url, data=json.dumps(payload).encode(), headers=UA)
    with urllib.request.urlopen(req, timeout=timeout) as r:
        return json.load(r)


def dgidb(genes: list[str]) -> list[dict]:
    q = """query($g:[String!]){genes(names:$g){nodes{name interactions{drug{name approved}
           interactionScore interactionTypes{type directionality}}}}}"""
    d = _post("https://dgidb.org/api/graphql", {"query": q, "variables": {"g": genes}})
    rows = []
    for g in d["data"]["genes"]["nodes"]:
        for it in g["interactions"]:
            rows.append({"source": "DGIdb", "gene": g["name"], "drug": it["drug"]["name"], "approved": it["drug"].get("approved"),
                         "type": ";".join(t["type"] or "" for t in it.get("interactionTypes", [])), "score": it.get("interactionScore")})
    return rows


def ensembl_id(symbol: str) -> str | None:
    q = 'query($q:String!){search(queryString:$q,entityNames:["target"]){hits{id name}}}'
    d = _post("https://api.platform.opentargets.org/api/v4/graphql", {"query": q, "variables": {"q": symbol}})
    for h in d["data"]["search"]["hits"]:
        if h["name"].upper() == symbol.upper():
            return h["id"]
    return None


def opentargets(symbol: str) -> list[dict]:
    eid = ensembl_id(symbol)
    if not eid:
        return []
    url = "https://api.platform.opentargets.org/api/v4/graphql"
    attempts = [
        ("drugAndClinicalCandidates", """query($id:String!){target(ensemblId:$id){drugAndClinicalCandidates{rows{maxClinicalStage drug{name}}}}}"""),
        ("knownDrugs", """query($id:String!){target(ensemblId:$id){knownDrugs{rows{phase status drug{name} disease{name}}}}}"""),
    ]
    last = None
    for name, q in attempts:
        try:
            d = _post(url, {"query": q, "variables": {"id": eid}})
            if d.get("errors"):
                last = d["errors"][0]["message"]; continue
            rows = d["data"]["target"][name]["rows"]
            return [{"source": "OpenTargets", "gene": symbol, "drug": r["drug"]["name"],
                     "stage": r.get("maxClinicalStage") or r.get("phase"), "disease": (r.get("disease") or {}).get("name")} for r in rows]
        except Exception as e:                       # record and try the next schema
            last = repr(e)
    raise RuntimeError(f"Open Targets: no schema variant worked ({last})")


def enrichr(genes: list[str], libs=("Reactome_2022", "GO_Biological_Process_2023", "KEGG_2021_Human")) -> list[dict]:
    boundary = "----mva"
    body = (f'--{boundary}\r\nContent-Disposition: form-data; name="list"\r\n\r\n' + "\n".join(genes) +
            f'\r\n--{boundary}\r\nContent-Disposition: form-data; name="description"\r\n\r\nmva\r\n--{boundary}--\r\n').encode()
    req = urllib.request.Request("https://maayanlab.cloud/Enrichr/addList", data=body,
                                 headers={"Content-Type": f"multipart/form-data; boundary={boundary}", "User-Agent": UA["User-Agent"]})
    uid = json.load(urllib.request.urlopen(req, timeout=60))["userListId"]
    rows = []
    for lib in libs:
        url = f"https://maayanlab.cloud/Enrichr/enrich?userListId={uid}&backgroundType={urllib.parse.quote(lib)}"
        d = json.load(urllib.request.urlopen(urllib.request.Request(url, headers={"User-Agent": UA["User-Agent"]}), timeout=60))
        for r in d[lib][:8]:      # rank, term, p, z, combined, genes, adj_p, ...
            rows.append({"library": lib, "term": r[1], "p": r[2], "adj_p": r[6], "genes": ";".join(r[5])})
        time.sleep(0.5)
    return rows


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", default="work_track2")
    a = ap.parse_args(argv)
    out = Path(a.out); out.mkdir(parents=True, exist_ok=True)
    import pandas as pd
    errors, drugs, enr = [], [], []
    allg = sorted({g for v in GENE_SETS.values() for g in v})
    try:
        drugs += dgidb(allg)
    except Exception as e:
        errors.append(f"DGIdb: {e!r}")
    for g in allg:
        try:
            drugs += opentargets(g)
        except Exception as e:
            errors.append(f"OpenTargets {g}: {e!r}")
    for name, genes in GENE_SETS.items():
        try:
            enr += [{"set": name, **r} for r in enrichr(genes)]
        except Exception as e:
            errors.append(f"Enrichr {name}: {e!r}")
    pd.DataFrame(drugs).to_csv(out / "db_drug_gene.csv", index=False)
    pd.DataFrame(enr).to_csv(out / "db_enrichment.csv", index=False)
    (out / "db_errors.json").write_text(json.dumps(errors, indent=1))
    print("drug-gene rows:", len(drugs), "| enrichment rows:", len(enr), "| errors:", len(errors))
    for e in errors[:10]:
        print(" ", e)
    if drugs:
        df = pd.DataFrame(drugs)
        print(df.groupby("drug").agg(genes=("gene", lambda s: ",".join(sorted(set(s)))), n=("gene", "nunique"))
              .sort_values("n", ascending=False).head(40).to_string())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
