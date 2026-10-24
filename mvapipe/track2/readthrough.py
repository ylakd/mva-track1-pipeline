"""Variant-specific check of eligibility for translational readthrough (nonsense variants only).

What is computed: which stop codon the variant creates, the base at +4 (immediately 3' of the codon), the
qualitative readthrough propensity from the published hierarchies (UGA > UAG > UAA; a C at +4 favours
readthrough, a G disfavours it), and whether nonsense-mediated decay is predicted to remove the transcript
first. These hierarchies are population-level tendencies, not predictions for a given allele.
"""
from __future__ import annotations

import itertools
import json
import urllib.request

BASES = "TCAG"
AA = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
CODON = {a + b + c: AA[i] for i, (a, b, c) in enumerate(itertools.product(BASES, repeat=3))}
THREE = {"Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C", "Gln": "Q", "Glu": "E", "Gly": "G", "His": "H",
         "Ile": "I", "Leu": "L", "Lys": "K", "Met": "M", "Phe": "F", "Pro": "P", "Ser": "S", "Thr": "T",
         "Trp": "W", "Tyr": "Y", "Val": "V"}
COMP = str.maketrans("ACGT", "TGCA")


def candidate_stop_codons(ref_aa3: str, ref: str, alt: str, strand: str) -> list[dict]:
    """All codons for amino acid `ref_aa3` in which substituting `ref`>`alt` (genomic, given strand) makes a
    stop. Returns codon, index of the changed base in the codon (0-2, 5'->3' on the transcript) and the stop."""
    aa = THREE[ref_aa3]
    r, a = (ref, alt) if strand == "+" else (ref.translate(COMP), alt.translate(COMP))
    out = []
    for cod, x in CODON.items():
        if x != aa:
            continue
        for i in range(3):
            if cod[i] == r:
                new = cod[:i] + a + cod[i + 1:]
                if CODON[new] == "*":
                    out.append({"codon": cod, "index": i, "stop": new.replace("T", "U")})
    return out


def plus4_coordinate(pos: int, index: int, strand: str) -> int:
    """Genomic (1-based) coordinate of the base at +4, i.e. 3' of the stop codon on the transcript."""
    return pos - index + 3 if strand == "+" else pos + index - 3


def minus1_coordinate(pos: int, index: int, strand: str) -> int:
    """Genomic coordinate of the base immediately 5' of the stop codon on the transcript (position -1)."""
    return pos - index - 1 if strand == "+" else pos + index + 1


def propensity(stop: str, plus4: str | None, minus1: str | None = None) -> dict:
    """Heuristic, additive score from the context features reported by Floquet et al. 2012 (PLoS Genet):
    basal readthrough UGA > UAG > UAA; a C at +4 raises readthrough; a U at -1 is the main determinant of the
    gentamicin response. It orders contexts qualitatively and is not a prediction of efficacy."""
    s = {"UGA": 2, "UAG": 1, "UAA": 0}[stop]
    p = {"C": 1, "U": 0, "T": 0, "A": 0, "G": -1}.get(plus4 or "", 0)
    m = 1 if (minus1 or "") in ("T", "U") else 0
    total = s + p + m
    level = "higher" if total >= 3 else "intermediate" if total >= 1 else "lower"
    return {"stop_score": s, "plus4_score": p, "minus1_score": m, "level": level}


def fetch_base(chrom: str, pos: int, genome: str = "hg38", timeout: int = 20) -> str:
    """One reference base from the UCSC sequence API (needs internet). Returns upper-case A/C/G/T."""
    c = chrom if str(chrom).startswith("chr") else f"chr{chrom}"
    url = f"https://api.genome.ucsc.edu/getData/sequence?genome={genome};chrom={c};start={pos - 1};end={pos}"
    with urllib.request.urlopen(url, timeout=timeout) as r:
        return json.load(r)["dna"].upper()


def assess(ref_aa3: str, ref: str, alt: str, strand: str, nmd_predicted: bool,
           plus4_base: str | None = None, minus1_base: str | None = None) -> dict:
    cands = candidate_stop_codons(ref_aa3, ref, alt, strand)
    res = {"candidates": cands, "unique": len(cands) == 1, "nmd_predicted": bool(nmd_predicted)}
    if len(cands) != 1:
        res["verdict"] = "stop codon identity ambiguous from the amino-acid change alone"
        return res
    c = cands[0]
    plus4 = plus4_base.upper() if plus4_base else None
    if plus4 and strand == "-":
        plus4 = plus4.translate(COMP)
    minus1 = minus1_base.upper() if minus1_base else None
    if minus1 and strand == "-":
        minus1 = minus1.translate(COMP)
    prop = propensity(c["stop"], plus4, minus1)
    res.update({"stop": c["stop"], "plus4": plus4, "minus1": minus1, "propensity": prop})
    if nmd_predicted:
        res["verdict"] = ("limited: the transcript is predicted to be degraded by NMD before translation, so "
                          "readthrough alone cannot restore expression (no approved NMD inhibitor)")
    else:
        res["verdict"] = f"candidate for readthrough ({prop['level']} propensity)"
    return res
