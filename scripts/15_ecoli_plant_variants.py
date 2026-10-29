#!/usr/bin/env python3
# Plant known SNPs + small indels into the E. coli assembly (strain's own genome),
# producing a mutated reference R' and a truth VCF in R' coordinates.
# Real SRR9900640 reads (which match the ORIGINAL assembly) mapped to R' must then
# call exactly these planted variants (reversed): a real-read benchmark with known truth.
import random, sys

random.seed(20260716)                      # fixed seed (reproducible)
ASM = "ecoli_real2/flye/assembly.fasta"
OUT_REF = "ecoli_vc/ref_mut.fa"
OUT_VCF = "ecoli_vc/truth.vcf"
SNP_EVERY = 700                            # ~1 SNP per 700 bp
INDEL_EVERY = 3500                         # ~1 indel per 3500 bp
MARGIN = 60                                # keep away from contig ends / each other
BASES = "ACGT"

def read_fasta(p):
    seqs, name = {}, None
    for line in open(p):
        line = line.rstrip()
        if line.startswith(">"):
            name = line[1:].split()[0]; seqs[name] = []
        else:
            seqs[name].append(line.upper())
    return {k: "".join(v) for k, v in seqs.items()}

asm = read_fasta(ASM)
import os
os.makedirs("ecoli_vc", exist_ok=True)
ref_out = open(OUT_REF, "w")
vcf = open(OUT_VCF, "w")
vcf.write("##fileformat=VCFv4.2\n")
for c, s in asm.items():
    vcf.write(f"##contig=<ID={c},length={len(s)}>\n")
vcf.write('##FILTER=<ID=PASS,Description="All filters passed">\n')
vcf.write('##INFO=<ID=.,Number=0,Type=Flag,Description=".">\n')
vcf.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
vcf.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n")

records = []   # (contig, mutpos_1based, REF, ALT)
for c, s in asm.items():
    L = len(s)
    if L < 5000:                           # skip tiny contigs
        ref_out.write(f">{c}\n{s}\n"); continue
    # choose sites in ORIGINAL coords, spaced, non-overlapping
    sites = {}   # origpos -> ("snp") or ("del",k) or ("ins",seq)
    i = MARGIN
    while i < L - MARGIN:
        r = random.random()
        step = SNP_EVERY
        if r < SNP_EVERY / INDEL_EVERY:    # occasionally an indel instead
            if random.random() < 0.5:
                k = random.randint(1, 5); sites[i] = ("del", k); step = INDEL_EVERY; i += k
            else:
                ins = "".join(random.choice(BASES) for _ in range(random.randint(1, 5)))
                sites[i] = ("ins", ins); step = INDEL_EVERY
        else:
            b = s[i]
            if b in BASES:
                nb = random.choice([x for x in BASES if x != b]); sites[i] = ("snp", nb)
        i += step + random.randint(0, 200)
    # build mutated sequence, emit truth in MUTATED coords
    out = []; mut_len = 0
    def emit(ch):
        global mut_len
        out.append(ch); mut_len += len(ch)
    j = 0
    while j < L:
        if j in sites:
            kind = sites[j]
            if kind[0] == "snp":
                emit(kind[1])                                  # mutated base in R'
                records.append((c, mut_len, kind[1], s[j]))    # REF=R', ALT=reads(orig)
                j += 1
            elif kind[0] == "del":                             # delete k orig bases from R'
                k = kind[1]                                     # reads have them -> INSERTION in reads
                anchor = out[-1] if out else s[j-1]            # last emitted base (anchor, already in R')
                ins_bases = s[j:j+k]
                # VCF insertion at anchor (mut coord = current mut_len): REF=anchor, ALT=anchor+ins
                records.append((c, mut_len, anchor, anchor + ins_bases))
                j += k                                          # skip these in R'
            else:                                              # insert bases into R' not in reads
                ins = kind[1]                                   # reads lack them -> DELETION in reads
                anchor = s[j-1] if j > 0 else s[0]
                # emit anchor already emitted; now add inserted bases to R'
                start_anchor_mutpos = mut_len                   # anchor is at mut_len (last emitted)
                emit(ins)
                # VCF deletion: REF=anchor+ins (in R'), ALT=anchor
                records.append((c, start_anchor_mutpos, anchor + ins, anchor))
                emit(s[j]); j += 1
        else:
            emit(s[j]); j += 1
    seq = "".join(out)
    ref_out.write(f">{c}\n")
    for x in range(0, len(seq), 80):
        ref_out.write(seq[x:x+80] + "\n")
ref_out.close()

# sort records by contig order then pos, write VCF
order = {c: idx for idx, c in enumerate(asm)}
records.sort(key=lambda r: (order[r[0]], r[1]))
nsnp = nind = 0
for c, pos, ref, alt in records:
    if len(ref) == 1 and len(alt) == 1: nsnp += 1
    else: nind += 1
    vcf.write(f"{c}\t{pos}\t.\t{ref}\t{alt}\t50\tPASS\t.\tGT\t1/1\n")
vcf.close()
print(f"planted SNPs={nsnp} indels={nind} total={len(records)}")
