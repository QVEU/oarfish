'''
Offline test for polyprotein_peptides.py: builds a synthetic genome with NCBI-style GFF3 holding
(1) a plus-strand polyprotein, (2) a -1 ribosomal frameshift polyprotein (pp1ab-like join with a
1-nt overlap, and a mature protein spanning the frameshift) sharing its start with a pp1a-like CDS,
and (3) a minus-strand polyprotein whose mature regions have no Parent attribute.
Checks every mature protein comes out with the right aa range and sequence, from both the
--genome-fasta (translation) and --protein-fasta routes.
Run: python tests/test_polyprotein_peptides.py   (or pytest)
'''
import os, random, subprocess, sys, tempfile
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, ".."))
import polyprotein_peptides as pp

TABLE = pp.codonTable(1)
SENSE = [c for c, a in TABLE.items() if a != "*"]

def randCodons(rng, n, first=None):
    out = []
    for i in range(n):
        pool = [c for c in SENSE if first is None or i > 0 or c[0] == first]
        out.append(rng.choice(pool))
    return out

def aa(codons):
    return "".join(TABLE[c] for c in codons)

def gffLine(seqid, ftype, s, e, strand, attrs, phase="."):
    return "\t".join([seqid, "RefSeq", ftype, str(s), str(e), ".", strand, phase, attrs])

def build(rng):
    seqid = "NC_TEST.1"
    lines, truth, prots = ["##gff-version 3"], [], {}
    genome = "".join(rng.choice("ACGT") for _ in range(50))          # 5' UTR

    # (1) plus-strand polyprotein: 3 mature proteins, VP-like product names with spaces/%2C
    start = len(genome) + 1
    sizes = [30, 45, 25]
    cod = randCodons(rng, sum(sizes))
    P = aa(cod)
    genome += "".join(cod) + "TAA"
    end = len(genome)
    prots["NP_041277.1"] = P
    lines.append(gffLine(seqid, "CDS", start, end, "+", "ID=cds-NP_041277.1;Name=NP_041277.1;product=polyprotein;protein_id=NP_041277.1", "0"))
    a = 1
    for k, n in enumerate(sizes):
        s = start + 3 * (a - 1)
        e = s + 3 * n - 1
        name = ["VP4", "VP2%2C capsid", "2A protease"][k]
        lines.append(gffLine(seqid, "mature_protein_region_of_CDS", s, e, "+",
                             f"ID=id-NP_041277.1:{a}..{a+n-1};Parent=cds-NP_041277.1;product={name};protein_id=YP_00000{k}.1"))
        truth.append(("NP_041277.1", name.replace("%2C", ","), a, a + n - 1, P[a - 1:a + n - 1]))
        a += n
    genome += "".join(rng.choice("ACGT") for _ in range(20))

    # (2) frameshift: pp1a = A + stop ; pp1ab = A (ending at slippery nt) then -1 frame B
    s1 = len(genome) + 1
    codA = randCodons(rng, 40)
    codB = randCodons(rng, 35, first=codA[-1][-1])                   # first codon of B rereads last nt of A
    seg1 = "".join(codA)
    seg2 = "".join(codB)
    genome += seg1 + seg2[1:] + "TAG"
    e1 = s1 + len(seg1) - 1
    s2 = e1                                                             # 1-nt overlap
    e2 = s2 + len(seg2) + 3 - 1
    PAB = aa(codA) + aa(codB)
    prots["YP_009724389.1"] = PAB
    lines.append(gffLine(seqid, "CDS", s1, e1, "+", "ID=cds-YP_009724389.1;Name=YP_009724389.1;product=ORF1ab polyprotein;protein_id=YP_009724389.1", "0"))
    lines.append(gffLine(seqid, "CDS", s2, e2, "+", "ID=cds-YP_009724389.1;Name=YP_009724389.1;product=ORF1ab polyprotein;protein_id=YP_009724389.1", "0"))
    # pp1a: same start, reads A then continues in frame 0 (we just annotate the first 30 aa as nsp1)
    PA = aa(codA)
    prots["YP_009725295.1"] = PA
    lines.append(gffLine(seqid, "CDS", s1, e1, "+", "ID=cds-YP_009725295.1;Name=YP_009725295.1;product=ORF1a polyprotein;protein_id=YP_009725295.1", "0"))
    lines.append(gffLine(seqid, "mature_protein_region_of_CDS", s1, s1 + 89, "+",
                         "ID=id-YP_009725295.1:1..30;Parent=cds-YP_009725295.1;product=nsp1"))
    truth.append(("YP_009725295.1", "nsp1", 1, 30, PA[:30]))
    lines.append(gffLine(seqid, "mature_protein_region_of_CDS", s1, s1 + 89, "+",
                         "ID=id-YP_009724389.1:1..30;Parent=cds-YP_009724389.1;product=nsp1"))
    truth.append(("YP_009724389.1", "nsp1", 1, 30, PAB[:30]))
    # nsp12-like: aa 31..60, spanning the frameshift (aa 40 is last of A, 41 first of B)
    ms = s1 + 90
    lines.append(gffLine(seqid, "mature_protein_region_of_CDS", ms, e1, "+",
                         "ID=id-YP_009724389.1:31..60;Parent=cds-YP_009724389.1;product=RNA-dependent RNA polymerase"))
    lines.append(gffLine(seqid, "mature_protein_region_of_CDS", s2, s2 + 3 * 20 - 1, "+",
                         "ID=id-YP_009724389.1:31..60;Parent=cds-YP_009724389.1;product=RNA-dependent RNA polymerase"))
    truth.append(("YP_009724389.1", "RNA-dependent RNA polymerase", 31, 60, PAB[30:60]))
    genome += "".join(rng.choice("ACGT") for _ in range(20))

    # (3) minus-strand polyprotein, mature regions without Parent
    codM = randCodons(rng, 50)
    PM = aa(codM)
    prots["NP_000002.1"] = PM
    coding = "".join(codM) + "TGA"
    ms = len(genome) + 1
    genome += pp.revcomp(coding)
    me = len(genome)
    lines.append(gffLine(seqid, "CDS", ms, me, "-", "ID=cds-NP_000002.1;protein_id=NP_000002.1", "0"))
    for name, a0, a1 in [("capsid", 1, 20), ("prM", 21, 50)]:
        hi = me - 3 * (a0 - 1)
        lo = me - 3 * a1 + 1
        lines.append(gffLine(seqid, "mat_peptide", lo, hi, "-", f"ID=mp-{name};product={name}"))
        truth.append(("NP_000002.1", name, a0, a1, PM[a0 - 1:a1]))
    genome += "".join(rng.choice("ACGT") for _ in range(30))
    return seqid, genome, "\n".join(lines) + "\n", prots, truth

def readTsv(path):
    with open(path) as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        return [dict(zip(hdr, l.rstrip("\n").split("\t"))) for l in fh]

def run(tmp, extra, accession=None):
    out = os.path.join(tmp, "out", "res")
    cmd = [sys.executable, os.path.join(HERE, "..", "polyprotein_peptides.py")] + ([accession] if accession else []) \
        + ["--gff3", os.path.join(tmp, "g.gff3"), "-o", out] + extra
    subprocess.run(cmd, check=True, capture_output=True, text=True)
    return readTsv(out + ".tsv")

def check(rows, truth):
    got = sorted((r["polyprotein_id"], r["product"], int(r["aa_start"]), int(r["aa_end"]), r["sequence"]) for r in rows)
    assert got == sorted(truth), (got, sorted(truth))

def test_all():
    rng = random.Random(1)
    seqid, genome, gff, prots, truth = build(rng)
    with tempfile.TemporaryDirectory() as tmp:
        with open(os.path.join(tmp, "g.gff3"), "w") as fh:
            fh.write(gff)
        with open(os.path.join(tmp, "g.fasta"), "w") as fh:
            fh.write(f">{seqid} synthetic\n{genome}\n")
        with open(os.path.join(tmp, "p.faa"), "w") as fh:
            fh.write("".join(f">{k} polyprotein\n{v}\n" for k, v in prots.items()))
        # translation route must reproduce the polyproteins exactly
        cdsById, _ = pp.parseGff3(gff)
        g = {seqid: genome}
        for cid, c in cdsById.items():
            assert pp.translateCds(c, g) == prots[c["attrs"]["protein_id"]], cid
        check(run(tmp, ["--genome-fasta", os.path.join(tmp, "g.fasta")]), truth)
        check(run(tmp, ["--protein-fasta", os.path.join(tmp, "p.faa")]), truth)
        # a protein accession restricts output to that polyprotein (version-insensitive)
        check(run(tmp, ["--protein-fasta", os.path.join(tmp, "p.faa")], accession="YP_009724389"),
              [t for t in truth if t[0] == "YP_009724389.1"])
    print(f"OK: {len(truth)} mature proteins across plus-strand, frameshifted and minus-strand polyproteins")

if __name__ == "__main__":
    test_all()
