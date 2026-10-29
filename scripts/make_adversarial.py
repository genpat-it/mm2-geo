#!/usr/bin/env python3
"""Plant dense indels into chr22 to stress mm2-geo's DP band. Emits a truth VCF
(GT 1/1) whose REF/ALT are built against the reference, so `bcftools consensus`
turns it into the adversarial sample genome and the same VCF is the truth."""
import random, sys, os

REF_FA = "chr22_named.fa"      # single contig named "chr22"
OUT_VCF = os.environ.get("ADV_OUT", "adv_chr22/truth_adv.vcf")
CONTIG = "chr22"
STEP = int(os.environ.get("ADV_STEP", "800"))   # spacing (> max L, no overlap)
LMIN = int(os.environ.get("ADV_LMIN", "10"))
LMAX = int(os.environ.get("ADV_LMAX", "60"))
SEED = int(os.environ.get("ADV_SEED", "22"))

random.seed(SEED)

# read the single-contig fasta
seq = []
with open(REF_FA) as fh:
    for line in fh:
        if line.startswith(">"):
            continue
        seq.append(line.strip())
seq = "".join(seq).upper()
N = len(seq)

# find the first and last non-N base to stay inside the assembled region
lo = 0
while lo < N and seq[lo] == "N":
    lo += 1
hi = N - 1
while hi > 0 and seq[hi] == "N":
    hi -= 1
lo += 1000; hi -= 1000  # margin

recs = []
p = lo
bases = "ACGT"
while p < hi:
    L = random.randint(LMIN, LMAX)
    # require a clean (N-free) window
    win = seq[p:p + L + 2]
    if "N" in win or len(win) < L + 2:
        p += STEP
        continue
    if random.random() < 0.5:
        # deletion of L bases after anchor p: REF = anchor + L deleted, ALT = anchor
        ref = seq[p:p + L + 1]
        alt = seq[p]
    else:
        # insertion of L random bases after anchor p
        ref = seq[p]
        alt = seq[p] + "".join(random.choice(bases) for _ in range(L))
    recs.append((p + 1, ref, alt))  # VCF is 1-based
    p += STEP

with open(OUT_VCF, "w") as out:
    out.write("##fileformat=VCFv4.2\n")
    out.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
    out.write(f"##contig=<ID={CONTIG},length={N}>\n")
    out.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n")
    for pos, ref, alt in recs:
        out.write(f"{CONTIG}\t{pos}\t.\t{ref}\t{alt}\t.\t.\t.\tGT\t1/1\n")

print(f"planted {len(recs)} indels into {CONTIG} (non-N region {lo}-{hi}), L in [{LMIN},{LMAX}], step {STEP}")
