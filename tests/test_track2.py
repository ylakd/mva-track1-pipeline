"""Unit tests for the Track 2 framework (no network, no patient data)."""
import json
import tempfile
from pathlib import Path

from mvapipe.track2 import faulttree as ft, fmea, readthrough as rt, run as t2


def small_tree():
    return {"top": "T", "nodes": {
        "T": {"type": "AND", "label": "t", "children": ["A", "B"]},
        "A": {"type": "OR", "label": "a", "children": ["a1", "a2"]},
        "B": {"type": "OR", "label": "b", "children": ["b1", "b2"]},
        "a1": {"type": "BASIC", "label": ""}, "a2": {"type": "BASIC", "label": ""},
        "b1": {"type": "BASIC", "label": ""}, "b2": {"type": "BASIC", "label": ""}}}


def test_cut_sets_and_or():
    cuts = ft.cut_sets(small_tree())
    assert len(cuts) == 4 and all(len(c) == 2 for c in cuts)
    assert frozenset({"a1", "b1"}) in cuts


def test_minimisation_removes_supersets():
    tree = {"top": "T", "nodes": {"T": {"type": "OR", "label": "", "children": ["x", "G"]},
            "G": {"type": "AND", "label": "", "children": ["x", "y"]},
            "x": {"type": "BASIC", "label": ""}, "y": {"type": "BASIC", "label": ""}}}
    assert ft.cut_sets(tree) == [frozenset({"x"})]


def test_coverage_and_hitting_sets():
    cuts = ft.cut_sets(small_tree())
    assert ft.coverage(cuts, {"a1"}) == 0.5
    assert ft.coverage(cuts, {"a1", "a2"}) == 1.0
    hs = ft.fully_mitigating_sets(cuts)
    assert frozenset({"a1", "a2"}) in hs and frozenset({"b1", "b2"}) in hs


def test_validate_flags_errors():
    bad = {"top": "Z", "nodes": {"T": {"type": "AND", "children": ["q"]}}}
    assert ft.validate(bad)


def test_occurrence_monotone():
    assert fmea.occurrence(5, 1.0) < fmea.occurrence(2, 0.0)
    assert 1 <= fmea.occurrence(5, 1.0) <= 10


def test_hard_stops_list_all_reasons():
    item = {"E": 1, "S": 10, "approval": "none"}
    assert len(fmea.hard_stops(item)) == 3
    ok = {"E": 3, "S": 5, "approval": "approved_adult"}
    assert fmea.hard_stop(ok) is None
    assert fmea.hard_stop({**ok, "opposite_mechanism": "x"})


def test_rank_orders_by_rpn_and_separates_stopped():
    rows = [{"id": "a", "RPN": 100, "S": 5, "stop": None}, {"id": "b", "RPN": 50, "S": 5, "stop": None},
            {"id": "c", "RPN": 1, "S": 9, "stop": "x"}]
    kept, out = fmea.rank(rows)
    assert [r["id"] for r in kept] == ["b", "a"] and out[0]["id"] == "c"


def test_rank_stability_deterministic_and_bounded():
    rows = [{"id": "a", "S": 4, "O": 5, "D": 3, "RPN": 60}, {"id": "b", "S": 8, "O": 6, "D": 5, "RPN": 240}]
    s1, s2 = fmea.rank_stability(rows, n=500), fmea.rank_stability(rows, n=500)
    assert s1 == s2 and 0 <= s1["p_same_order"] <= 1
    assert s1["per_candidate"]["a"]["p_rank1"] > 0.9


def test_stop_codon_candidates():
    # Leu TTA -> TGA on + strand (T>G at codon position 2)
    c = rt.candidate_stop_codons("Leu", "T", "G", "+")
    assert any(x["codon"] == "TTA" and x["stop"] == "UGA" for x in c)
    # reverse strand: complement
    c2 = rt.candidate_stop_codons("Leu", "A", "C", "-")
    assert any(x["stop"] == "UGA" for x in c2)


def test_plus4_coordinate():
    assert rt.plus4_coordinate(100, 0, "+") == 103
    assert rt.plus4_coordinate(100, 2, "+") == 101
    assert rt.plus4_coordinate(100, 0, "-") == 97


def test_assess_nmd_limits_readthrough():
    r = rt.assess("Leu", "T", "G", "+", True, plus4_base="C")
    assert r["verdict"].startswith("limited")
    r2 = rt.assess("Leu", "T", "G", "+", False, plus4_base="C")
    assert "candidate" in r2["verdict"]
    assert rt.propensity("UGA", "C")["level"] == "higher"
    assert rt.propensity("UAA", "G")["level"] == "lower"


def test_end_to_end_on_knowledge_base():
    kb = t2.load_kb()
    res = t2.analyse(kb)
    ids = [r["id"] for r in res["ranking"]]
    assert ids and all(i.startswith("I") for i in ids)
    assert {x["id"] for x in res["excluded"]} >= {"X1", "X8"}
    x8 = next(x for x in res["excluded"] if x["id"] == "X8")
    assert "opposite" in x8["stop"]
    assert res["fully_mitigating_sets_root"]
    with tempfile.TemporaryDirectory() as d:
        assert json.dumps(res)


def test_kb_has_no_coordinates():
    txt = Path(t2.__file__).resolve().parents[2].joinpath("knowledge", "bub1b_mva1_kb.json").read_text()
    import re
    assert not re.search(r"chr\d+[:_]\d+", txt) and not re.search(r"\bc\.\d+", txt)


def test_minus1_context_raises_propensity():
    assert rt.propensity("UGA", "C", "T")["level"] == "higher"
    assert rt.propensity("UAA", "A", "G")["minus1_score"] == 0
    assert rt.minus1_coordinate(100, 0, "+") == 99 and rt.minus1_coordinate(100, 0, "-") == 101


def test_engine_runs_on_a_second_knowledge_file():
    """Scalability: the same engine, a different genotype class, no code change."""
    kb = t2.load_kb(Path(t2.__file__).resolve().parents[2] / "knowledge" / "TEMPLATE_new_genotype.json")
    res = t2.analyse(kb)
    assert [r["id"] for r in res["ranking"]] == ["I1", "I2"] or [r["id"] for r in res["ranking"]] == ["I2", "I1"]
    assert res["excluded"][0]["id"] == "X1"


def test_every_cited_reference_exists():
    kb = t2.load_kb()
    for it in kb["interventions"]:
        assert set(it.get("refs", [])) <= set(kb["references"]), it["id"]


def test_scenarios_relevance_and_modifier():
    from mvapipe.track2 import scenarios as scn
    kb = t2.load_kb()
    cuts = {"root": ft.cut_sets(kb["root_tree"]), "consequence": ft.cut_sets(kb["consequence_tree"])}
    i1 = next(i for i in kb["interventions"] if i["id"] == "I1")
    i2 = next(i for i in kb["interventions"] if i["id"] == "I2")
    base = scn.score_under(i1, cuts, 1.0)
    low = scn.score_under(i1, cuts, 0.2)
    assert base["RPN"] <= low["RPN"] and base["RPN"] == 84          # root tier gets worse when r falls
    assert scn.score_under(i2, cuts, 0.2)["O"] == scn.score_under(i2, cuts, 1.0)["O"]   # consequence tier unchanged
    mod = kb["uncertainty"]["modifier"]
    assert scn.score_under(i2, cuts, 1.0, mod)["stop"]                # severity 9 -> hard stop
    assert scn.score_under(i1, cuts, 1.0, None)["S"] == 4
    tab = scn.scenario_table(kb, True)
    assert tab["r_grid"] and "I2" not in tab["minimax_order"]


def test_checkpoint_sim_directions():
    from mvapipe.track2 import checkpoint_sim as cs
    r = cs.run(2000)
    raw = r["_raw"]
    # dose-response: less BUBR1 never gives fewer errors (median curve non-increasing in b)
    dr = cs.dose_response(raw["params"])
    vals = [dr[k] for k in sorted(dr)]
    assert all(a >= b - 1e-9 for a, b in zip(vals, vals[1:]))
    assert r["baseline"]["p_disease_not_better"] == 1.0
    iv = r["interventions"]
    assert all(v["p_better"] == 0 for k, v in iv.items() if k.startswith(("MPS1", "Microtubule")))
    assert all(v["p_worse"] == 0 for k, v in iv.items() if k.startswith("BUBR1"))
    assert cs.run(500, seed=1)["interventions"] == cs.run(500, seed=1)["interventions"]


def test_db_triage_rules():
    import pandas as pd
    from mvapipe.track2 import db_triage as dt
    def row(src, gene, drug, approved=None, type=None, stage=None):
        return {"source": src, "gene": gene, "drug": drug, "approved": approved, "type": type, "stage": stage}
    df = pd.DataFrame([
        row("DGIdb", "JAK1", "baricitinib", True, "inhibitor"),
        row("DGIdb", "STING1", "agonistdrug", True, "agonist"),
        row("DGIdb", "AURKB", "mitoticdrug", True, "inhibitor"),
        row("DGIdb", "SRC", "dasatinib", True, "inhibitor"), row("DGIdb", "AURKB", "dasatinib", True, None),
        row("DGIdb", "NAMPT", "nampti", True, "inhibitor"),
        row("DGIdb", "SIRT1", "niacinlike", True, "activator"),
        row("DGIdb", "TP53", "textmined", True, None),
        row("DGIdb", "JAK2", "unapproved", False, "inhibitor"),
        row("DGIdb", "IFNAR1", "biologic", True, None), row("OpenTargets", "IFNAR1", "biologic", None, None, "PHASE_4")])
    r = dt.triage(df, t2.load_kb()).set_index("drug")
    assert r.loc["BARICITINIB", "status"] == "PASS" and r.loc["BARICITINIB", "in_knowledge_file"] == "I2"
    assert "opposite" in r.loc["AGONISTDRUG", "stop"]
    assert "mitotic" in r.loc["MITOTICDRUG", "stop"]
    assert r.loc["DASATINIB", "status"] == "PASS"                       # untyped AURKB row ignored ...
    assert "mitotic" in dt.triage(df, t2.load_kb(), True).set_index("drug").loc["DASATINIB", "stop"]   # ... unless strict
    assert "opposite" in r.loc["NAMPTI", "stop"] and r.loc["NIACINLIKE", "status"] == "PASS"
    assert "primary target" in r.loc["TEXTMINED", "stop"] and "not approved" in r.loc["UNAPPROVED", "stop"]
    assert r.loc["BIOLOGIC", "status"].startswith("REVIEW")
