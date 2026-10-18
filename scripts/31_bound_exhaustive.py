#!/usr/bin/env python3
# Exhaustive check of the score certificate (Supplementary, "Proof of the score certificate"): for small gaps,
# a DP over ALL paths (dual-affine gaps, runs split arbitrarily, every aligned pair scored +a, i.e. the most
# favourable case) finds the best path that leaves the band |i-j|<=b, and compares it with U(b).
# Scorings: HiFi preset, ONT preset, single affine; gaps with the read longer AND shorter than the reference.
import sys
NEG = -10**9
def run(a, O, label):
    g = lambda l: min(o + e*l for o, e in O)
    def U(q, t, b):
        d = abs(q - t); L = (b+1) + (b+1-d)
        return a*((q + t - L)//2) - g(b+1) - g(b+1-d)
    def best_out(q, t, b):
        S = {}; states = ['M'] + [f'I{k}' for k in range(len(O))] + [f'D{k}' for k in range(len(O))]
        get = lambda i, j, st, f: S.get((i, j, st, f), NEG)
        for i in range(q+1):
            for j in range(t+1):
                out = 1 if abs(i-j) > b else 0
                for f in (0, 1):
                    prevfs = [pf for pf in (0, 1) if (pf or out) == f]
                    for st in states:
                        v = 0 if (i == 0 and j == 0 and st == 'M' and f == out) else NEG
                        for pf in prevfs:
                            if st == 'M' and i > 0 and j > 0:
                                v = max(v, max(get(i-1, j-1, s2, pf) for s2 in states) + a)
                            for k, (o, e) in enumerate(O):
                                if st == f'I{k}' and i > 0:
                                    v = max(v, max(get(i-1, j, s2, pf) for s2 in states if s2 != st) - o - e, get(i-1, j, st, pf) - e)
                                if st == f'D{k}' and j > 0:
                                    v = max(v, max(get(i, j-1, s2, pf) for s2 in states if s2 != st) - o - e, get(i, j-1, st, pf) - e)
                        S[(i, j, st, f)] = v
        return max(get(q, t, st, 1) for st in states)
    n = bad = tight = 0
    for q in range(20, 46, 5):
        for t in list(range(q-10, q+1, 3)) + list(range(q+1, q+11, 3)):
            d = abs(q - t)
            for b in range(max(d, 1), min(q, t)//2):
                bo, u = best_out(q, t, b), U(q, t, b); n += 1
                if bo > u: bad += 1; print("VIOLATION", label, q, t, b, bo, u)
                if bo == u: tight += 1
    print(f"{label}\tcases={n}\tviolations={bad}\ttight={tight}", flush=True)
run(1, [(6, 2), (26, 1)], "hifi_preset")
run(2, [(4, 2), (24, 1)], "ont_preset")
run(1, [(5, 4)], "single_affine")
