"""Minimal fault-tree engine (AND / OR gates over basic events) with minimal cut sets and coverage.

A cut set is a set of basic events whose joint occurrence produces the top event. An intervention that
clears basic event(s) removes every cut set that contains them. 'Coverage' of an intervention is the
fraction of minimal cut sets it breaks; an intervention that breaks all of them fully mitigates the tree.
"""
from __future__ import annotations

from itertools import product


def _minimise(sets: list[frozenset]) -> list[frozenset]:
    sets = sorted(set(sets), key=len)
    out: list[frozenset] = []
    for s in sets:
        if not any(o <= s for o in out):
            out.append(s)
    return out


def cut_sets(tree: dict, node: str | None = None) -> list[frozenset]:
    """tree = {"top": id, "nodes": {id: {"type": "AND"|"OR"|"BASIC", "children": [ids], "label": str}}}"""
    node = node or tree["top"]
    spec = tree["nodes"][node]
    kind = spec["type"]
    if kind == "BASIC":
        return [frozenset([node])]
    child = [cut_sets(tree, c) for c in spec["children"]]
    if kind == "OR":
        return _minimise([s for c in child for s in c])
    if kind == "AND":
        combos = [frozenset().union(*pick) for pick in product(*child)]
        return _minimise(combos)
    raise ValueError(f"unknown gate type {kind!r}")


def coverage(cuts: list[frozenset], cleared: set[str]) -> float:
    """Fraction of minimal cut sets broken when the basic events in `cleared` are removed."""
    if not cuts:
        return 0.0
    return sum(1 for c in cuts if c & cleared) / len(cuts)


def fully_mitigating_sets(cuts: list[frozenset]) -> list[frozenset]:
    """Smallest sets of basic events that hit every cut set (minimal hitting sets), by brute force.
    Trees here are tiny, so exhaustive search is fine."""
    events = sorted(set().union(*cuts)) if cuts else []
    best: list[frozenset] = []
    for r in range(1, len(events) + 1):
        for comb in __import__("itertools").combinations(events, r):
            cs = frozenset(comb)
            if all(c & cs for c in cuts) and not any(b <= cs for b in best):
                best.append(cs)
        if best:
            break
    return best


def validate(tree: dict) -> list[str]:
    problems = []
    nodes = tree.get("nodes", {})
    if tree.get("top") not in nodes:
        problems.append("top node missing")
    for k, v in nodes.items():
        if v["type"] in ("AND", "OR"):
            for c in v.get("children", []):
                if c not in nodes:
                    problems.append(f"{k}: unknown child {c}")
        elif v["type"] != "BASIC":
            problems.append(f"{k}: bad type {v['type']}")
    return problems
