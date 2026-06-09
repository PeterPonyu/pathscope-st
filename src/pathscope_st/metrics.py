from __future__ import annotations
from math import sqrt


def pearson(xs: list[float], ys: list[float]) -> float:
    if len(xs) != len(ys) or len(xs) < 2:
        raise ValueError("pearson requires paired vectors of length >= 2")
    mx = sum(xs) / len(xs)
    my = sum(ys) / len(ys)
    num = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    denx = sqrt(sum((x - mx) ** 2 for x in xs))
    deny = sqrt(sum((y - my) ** 2 for y in ys))
    return num / (denx * deny) if denx and deny else 0.0


def interval_coverage(truth: list[float], lower: list[float] | None, upper: list[float] | None) -> float | None:
    if lower is None or upper is None:
        return None
    if len(truth) != len(lower) or len(truth) != len(upper):
        raise ValueError("interval vectors must align")
    return sum(lo <= y <= hi for y, lo, hi in zip(truth, lower, upper)) / len(truth)
