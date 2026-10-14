"""Shared helpers."""
from __future__ import annotations

import gzip
import hashlib
import json
import os
import re
from typing import Iterator

import numpy as np

ATTR_RE = re.compile(r'(\w+) "([^"]*)"')
GENE_RE = re.compile(r'gene_name "([^"]+)"')


def as_list(x):
    if x is None:
        return []
    return x if isinstance(x, list) else [x]


def flat_nums(x) -> list[float]:
    """Flatten any mixture of numbers / lists / dicts with a 'score' key into floats."""
    out: list[float] = []
    for v in as_list(x):
        if isinstance(v, dict):
            out += flat_nums(v.get("score"))
        elif isinstance(v, list):
            out += flat_nums(v)
        elif isinstance(v, (int, float)) and not isinstance(v, bool):
            out.append(float(v))
    return out


def maxnum(x) -> float:
    n = flat_nums(x)
    return max(n) if n else float("nan")


def hgvs_genomic(chrom: str, pos: int, ref: str, alt: str) -> str | None:
    """Genomic HGVS string (hg38) for MyVariant.info. Multi-allelic records return None."""
    chrom = str(chrom).replace("chr", "")
    if "," in alt:
        return None
    if len(ref) == len(alt) == 1:
        return f"chr{chrom}:g.{pos}{ref}>{alt}"
    if len(ref) > len(alt) and ref.startswith(alt):
        s, e = pos + len(alt), pos + len(ref) - 1
        return f"chr{chrom}:g.{s}del" if s == e else f"chr{chrom}:g.{s}_{e}del"
    if len(alt) > len(ref) and alt.startswith(ref):
        q = pos + len(ref) - 1
        return f"chr{chrom}:g.{q}_{q + 1}ins{alt[len(ref):]}"
    if len(ref) == 1:
        return f"chr{chrom}:g.{pos}delins{alt}"
    return f"chr{chrom}:g.{pos}_{pos + len(ref) - 1}delins{alt}"


def merge_intervals(intervals: list[tuple[int, int]]) -> tuple[np.ndarray, np.ndarray]:
    if not intervals:
        return np.array([], dtype=np.int64), np.array([], dtype=np.int64)
    intervals = sorted(intervals)
    merged = [list(intervals[0])]
    for s, e in intervals[1:]:
        if s <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], e)
        else:
            merged.append([s, e])
    a = np.array(merged, dtype=np.int64)
    return a[:, 0], a[:, 1]


def in_intervals(starts: np.ndarray, ends: np.ndarray, pos: np.ndarray) -> np.ndarray:
    if len(starts) == 0:
        return np.zeros(len(pos), dtype=bool)
    i = np.searchsorted(starts, pos, side="right") - 1
    return (i >= 0) & (pos <= ends[np.clip(i, 0, None)])


def gtf_iter(gtf_path, features: set[str] | None = None) -> Iterator[tuple]:
    opener = gzip.open if str(gtf_path).endswith(".gz") else open
    with opener(gtf_path, "rt") as f:
        for line in f:
            if line[0] == "#":
                continue
            c = line.rstrip("\n").split("\t")
            if len(c) < 9:
                continue
            if features is not None and c[2] not in features:
                continue
            yield c[0], c[2], int(c[3]), int(c[4]), c[6], c[8]


def parse_attrs(attr_string: str) -> dict:
    d: dict = {}
    for k, v in ATTR_RE.findall(attr_string):
        if k in d:
            if not isinstance(d[k], list):
                d[k] = [d[k]]
            d[k].append(v)
        else:
            d[k] = v
    return d


def sha256_file(path, chunk: int = 1 << 20) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        while True:
            b = f.read(chunk)
            if not b:
                break
            h.update(b)
    return h.hexdigest()


def write_json(obj, path) -> None:
    with open(path, "w") as f:
        json.dump(obj, f, indent=2, default=str)


def read_json(path):
    with open(path) as f:
        return json.load(f)


def hf_token() -> str | None:
    tok = os.environ.get("HF_TOKEN")
    if tok:
        return tok
    try:
        from kaggle_secrets import UserSecretsClient  # type: ignore
        return UserSecretsClient().get_secret("HF_TOKEN")
    except Exception:
        return None


def fmt(v, tag: str):
    """Read a FORMAT field from a cyvcf2 Variant; None when the tag is not declared in the header
    (cyvcf2 raises KeyError for undeclared tags, e.g. PGT/PID in an unphased callset)."""
    try:
        return v.format(tag)
    except KeyError:
        return None
