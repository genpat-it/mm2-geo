#!/usr/bin/env python3
"""Stratify small-variant F1 by indel length from rtg vcfeval output dirs.
For each eval dir we read tp-baseline.vcf.gz (TP, truth side), fn.vcf.gz (FN, truth side),
fp.vcf.gz (FP, call side) and bin every record by indel length L=abs(len(REF)-len(ALT)).
SNV = both alleles length 1. Bins: 1, 2-5, 6-15, 16-50 bp."""
import gzip, sys, os

BINS = [("SNV", None), ("1", (1, 1)), ("2-5", (2, 5)), ("6-15", (6, 15)), ("16-50", (16, 50))]

def classify(ref, alt):
    # take first ALT allele (rtg outputs are decomposed/biallelic in practice)
    alt = alt.split(",")[0]
    if len(ref) == 1 and len(alt) == 1:
        return "SNV", 0
    return "INDEL", abs(len(ref) - len(alt))

def bin_of(L):
    for name, rng in BINS:
        if name == "SNV":
            continue
        lo, hi = rng
        if lo <= L <= hi:
            return name
    return None  # >50 or 0-length complex

def count(path):
    """return dict: {'SNV': n, '1': n, '2-5': n, ...}"""
    d = {name: 0 for name, _ in BINS}
    if not os.path.exists(path):
        return d
    with gzip.open(path, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            c = line.rstrip("\n").split("\t")
            if len(c) < 5:
                continue
            ref, alt = c[3], c[4]
            kind, L = classify(ref, alt)
            if kind == "SNV":
                d["SNV"] += 1
            else:
                b = bin_of(L)
                if b:
                    d[b] += 1
    return d

def f1(tp_call, fp, tp_base, fn):
    # GA4GH convention (matches the other scripts in this suite):
    #   precision uses call-side TP (tp.vcf.gz); recall uses truth-side TP (tp-baseline.vcf.gz).
    p = tp_call / (tp_call + fp) if (tp_call + fp) else 0.0
    r = tp_base / (tp_base + fn) if (tp_base + fn) else 0.0
    return (2 * p * r / (p + r)) if (p + r) else 0.0, p, r

def eval_dir(ev):
    tp_base = count(os.path.join(ev, "tp-baseline.vcf.gz"))  # truth side -> recall
    tp_call = count(os.path.join(ev, "tp.vcf.gz"))           # call side  -> precision
    fp = count(os.path.join(ev, "fp.vcf.gz"))
    fn = count(os.path.join(ev, "fn.vcf.gz"))
    out = {}
    for name, _ in BINS:
        F, P, R = f1(tp_call[name], fp[name], tp_base[name], fn[name])
        # report truth-side TP with call-side FP and truth-side FN (TRUTH.TP / QUERY.FP / TRUTH.FN)
        out[name] = dict(tp=tp_base[name], fp=fp[name], fn=fn[name], F1=F, P=P, R=R)
    return out

def main():
    root = sys.argv[1] if len(sys.argv) > 1 else "strat"
    techs = ["hifi", "ont"]
    aligners = ["mm2", "geo"]
    print(f"{'tech':6} {'bin':7} | {'mm2 F1':>8} {'geo F1':>8} {'dF1':>8} | "
          f"{'mm2 TP/FP/FN':>16} {'geo TP/FP/FN':>16}")
    print("-" * 96)
    rows = []
    for tech in techs:
        res = {}
        for a in aligners:
            ev = os.path.join(root, f"{tech}_{a}", "ev")
            res[a] = eval_dir(ev)
        for name, _ in BINS:
            m, g = res["mm2"][name], res["geo"][name]
            d = g["F1"] - m["F1"]
            print(f"{tech:6} {name:7} | {m['F1']:8.4f} {g['F1']:8.4f} {d:+8.4f} | "
                  f"{m['tp']:5d}/{m['fp']:4d}/{m['fn']:4d}    "
                  f"{g['tp']:5d}/{g['fp']:4d}/{g['fn']:4d}")
            rows.append((tech, name, m, g, d))
        print("-" * 96)
    # LaTeX-friendly dump
    print("\n=== LaTeX rows (tech & bin & mm2 F1 & geo F1 & dF1) ===")
    for tech, name, m, g, d in rows:
        print(f"{tech} & {name} & {m['F1']:.4f} & {g['F1']:.4f} & {d:+.4f} \\\\")

if __name__ == "__main__":
    main()
