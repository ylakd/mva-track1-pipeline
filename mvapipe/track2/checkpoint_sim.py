"""Qualitative, parameter-free-of-data simulation of the spindle assembly checkpoint (SAC) with BUBR1 dosage.

THIS IS A TOY MODEL, NOT A PREDICTOR. Its parameters are not fitted to data. What it is used for is a *direction-of-effect*
test: for broad, random parameter sets, does an intervention that acts on mitosis raise or lower chromosome
mis-segregation in a checkpoint-weakened cell? A conclusion is reported only if it holds across the parameter ranges.

Model (per cell, expected values). N chromosomes each start with an attachment error with probability e0; an error is
corrected at rate kc. Checkpoint signal while E(t) = N e0 exp(-kc t) errors remain:  M(t) = g(b) * (1 - m) * E(t),
with g(b) = b**h (b = BUBR1 level relative to normal, h = dosage sensitivity) and m = fractional MPS1 inhibition.
Anaphase starts when M(t) < theta (checkpoint satisfied) but not before a fixed minimum time T0, and the cell slips out of
mitosis at T_max whatever the state. Mis-segregation = expected number of uncorrected errors at anaphase onset.
Not modelled: DNA damage, cell death, spatial effects, the specific biochemistry of MCC assembly, tissues.
"""
from __future__ import annotations

import numpy as np

N_CHR = 46


def draw_params(n: int, seed: int = 2026) -> dict:
    rng = np.random.default_rng(seed)
    lu = lambda lo, hi: np.exp(rng.uniform(np.log(lo), np.log(hi), n))
    return {"e0": lu(0.03, 0.3), "kc": lu(0.02, 0.3), "theta": lu(0.1, 3.0), "T0": lu(5, 30),
            "Tmax": lu(100, 600), "h": rng.uniform(1.0, 3.0, n)}


def missegregation(p: dict, b, mps1_inhib=0.0, err_mult=1.0, corr_mult=1.0):
    """Expected uncorrected errors at anaphase, and mitotic duration, for each parameter set."""
    e0 = np.clip(p["e0"] * err_mult, 0, 1)
    kc = p["kc"] * corr_mult
    s = (np.asarray(b, dtype=float) ** p["h"]) * (1.0 - mps1_inhib)
    n_err0 = N_CHR * e0
    need = np.where(n_err0 * s > p["theta"], np.log(np.maximum(n_err0 * s, 1e-12) / p["theta"]) / kc, 0.0)
    T = np.minimum(np.maximum(p["T0"], need), p["Tmax"])
    return n_err0 * np.exp(-kc * T), T


# interventions on mitosis: (label, kwargs applied to the DISEASE state)
INTERVENTIONS = {
    "BUBR1 level +50% (I1-like)": {"b_mult": 1.5},
    "BUBR1 level +100% (I1-like, upper bound)": {"b_mult": 2.0},
    "MPS1 inhibitor, 30% inhibition": {"mps1_inhib": 0.3},
    "MPS1 inhibitor, 70% inhibition": {"mps1_inhib": 0.7},
    "Microtubule agent, low dose (X2-like)": {"err_mult": 1.5, "corr_mult": 0.7},
    "Microtubule agent, moderate dose (X2-like)": {"err_mult": 2.5, "corr_mult": 0.5},
}


def run(n: int = 5000, b_disease: float = 0.35, b_healthy: float = 1.0, seed: int = 2026) -> dict:
    p = draw_params(n, seed)
    base_d, T_d = missegregation(p, b_disease)
    base_h, T_h = missegregation(p, b_healthy)
    out = {"n_parameter_sets": n, "b_disease": b_disease, "model_note": "toy model, not fitted to data",
           "baseline": {"misseg_disease_median": float(np.median(base_d)), "misseg_healthy_median": float(np.median(base_h)),
                        "p_disease_worse_than_healthy": float(np.mean(base_d > base_h + 1e-9)),
                        "p_disease_not_better": float(np.mean(base_d >= base_h - 1e-9))},
           "interventions": {}}
    for name, k in INTERVENTIONS.items():
        bb = min(b_disease * k.get("b_mult", 1.0), 1.0)
        m_d, T = missegregation(p, bb, k.get("mps1_inhib", 0.0), k.get("err_mult", 1.0), k.get("corr_mult", 1.0))
        m_h, _ = missegregation(p, b_healthy, k.get("mps1_inhib", 0.0), k.get("err_mult", 1.0), k.get("corr_mult", 1.0))
        eps = 1e-9
        d_dis = m_d - base_d
        d_hea = m_h - base_h
        tied = np.abs(d_dis) <= eps
        out["interventions"][name] = {
            "p_worse": float(np.mean(d_dis > eps)), "p_better": float(np.mean(d_dis < -eps)), "p_no_change": float(np.mean(tied)),
            "median_change_disease": float(np.median(d_dis)),
            "median_relative_change_disease": float(np.median(d_dis / np.maximum(base_d, 1e-9))),
            "median_change_healthy": float(np.median(d_hea)),
            "p_disease_harmed_more_than_healthy": float(np.mean(d_dis > d_hea + eps)),
        }
    out["_raw"] = {"base": base_d, "params": p}
    return out


def dose_response(p: dict, bs=None) -> dict:
    bs = bs if bs is not None else np.linspace(0.1, 1.0, 10)
    return {float(b): float(np.median(missegregation(p, b)[0])) for b in bs}
