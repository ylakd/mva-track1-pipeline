"""FMEA-style triage of intervention hypotheses, with hard stops and a rank-stability analysis.

Scales (1-10, as in classical FMEA; lower RPN is better):
  S  severity   : harm to the patient if the intervention is given and the hypothesis is wrong.
  O  occurrence : likelihood that the hypothesis is wrong or ineffective. Derived, not typed in:
                  O = round(10 - 1.2*E - 2*c), where E is the evidence grade (1-5) and c the coverage (0-1).
  D  detection  : difficulty of noticing failure or harm in time (1 = a validated biomarker exists).
  RPN = S * O * D.

Hard stops (applied before ranking): evidence grade below 2; no approved route (availability); severity >= 9.
The S and D values are expert judgements stored with a rationale in the knowledge file; the rank-stability
analysis shows how far they can move before the order changes.
"""
from __future__ import annotations

import numpy as np

E_MIN = 2
S_STOP = 9
APPROVED = ("approved_paediatric", "approved_adult", "approved_off_label")


def clamp(x: float, lo: int = 1, hi: int = 10) -> int:
    return int(max(lo, min(hi, round(x))))


def occurrence(evidence: int, coverage: float) -> int:
    return clamp(10 - 1.2 * evidence - 2.0 * coverage)


def hard_stops(item: dict) -> list[str]:
    """All hard-stop reasons that apply (an item can fail on several grounds)."""
    why = []
    if item["E"] < E_MIN:
        why.append(f"insufficient evidence (grade {item['E']} < {E_MIN})")
    if item.get("approval", "none") not in APPROVED + ("supplement",):
        why.append("no approved or available agent")
    if item["S"] >= S_STOP:
        why.append(f"severity {item['S']} >= {S_STOP}")
    if item.get("opposite_mechanism"):
        why.append("acts in the opposite direction to the mechanism (" + item["opposite_mechanism"] + ")")
    return why


def hard_stop(item: dict) -> str | None:
    w = hard_stops(item)
    return "; ".join(w) if w else None


def score(item: dict, cov: float) -> dict:
    o = occurrence(item["E"], cov)
    return {"id": item["id"], "label": item["label"], "tier": item["tier"], "coverage": round(cov, 2),
            "E": item["E"], "S": item["S"], "O": o, "D": item["D"], "RPN": item["S"] * o * item["D"],
            "stop": hard_stop(item)}


def rank(rows: list[dict]) -> tuple[list[dict], list[dict]]:
    kept = sorted([r for r in rows if r["stop"] is None], key=lambda r: (r["RPN"], r["S"]))
    out = [r for r in rows if r["stop"] is not None]
    return kept, out


def rank_stability(rows: list[dict], n: int = 10000, jitter: int = 1, seed: int = 2026) -> dict:
    """Perturb S, O and D of every candidate by an integer in [-jitter, +jitter] and recompute the order."""
    rng = np.random.default_rng(seed)
    ids = [r["id"] for r in rows]
    base = np.array([[r["S"], r["O"], r["D"]] for r in rows], dtype=float)
    ranks = np.zeros((n, len(rows)), dtype=int)
    same = 0
    base_order = list(np.argsort([r["RPN"] for r in rows], kind="stable"))
    for k in range(n):
        noise = rng.integers(-jitter, jitter + 1, size=base.shape)
        v = np.clip(base + noise, 1, 10)
        rpn = v[:, 0] * v[:, 1] * v[:, 2]
        order = list(np.argsort(rpn, kind="stable"))
        for pos, idx in enumerate(order):
            ranks[k, idx] = pos + 1
        same += order == base_order
    return {
        "n": n, "jitter": jitter, "p_same_order": same / n,
        "per_candidate": {i: {"mean_rank": float(ranks[:, j].mean()),
                              "p_rank1": float((ranks[:, j] == 1).mean()),
                              "p_same_rank_as_base": float((ranks[:, j] == (base_order.index(j) + 1)).mean())}
                          for j, i in enumerate(ids)},
    }
