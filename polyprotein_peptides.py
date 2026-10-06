#!/usr/bin/env python
'''
polyprotein_peptides.py
v1.0

oarfish

    ><IIIII\\
 //IIIIIIIII/
 \IIIIº>

Patrick T. Dolan
Unit Chief, Quantitative Virology and Evolution Unit

Splits a viral polyprotein into its mature proteins using the GFF3 annotation.

Give it an NCBI accession for either the genome (e.g. NC_001612.1) or the polyprotein itself
(e.g. NP_041277.1). It downloads the GFF3 for the genome from NCBI, finds every mature protein
region (mature_protein_region_of_CDS / mat_peptide) and its parent polyprotein CDS, converts the
genomic coordinates to positions in the polyprotein, and writes the protein sequence of each region.

Coordinates are mapped through the spliced CDS, so ribosomal frameshifts (e.g. coronavirus pp1ab,
where the CDS and nsp12 are join()s that overlap at the slippery site) come out right.

Polyprotein sequences are the NCBI protein records for the CDS protein_id (so readthrough and
transl_except are already applied). Offline, give --gff3 plus --protein-fasta, or --gff3 plus
--genome-fasta to translate the CDS instead (standard code, or the CDS's transl_table if set).

Output:
    <prefix>.faa   one FASTA record per mature protein:
                   >polyprotein_id|product|aaStart-aaEnd [mature_id=...] [region=seqid:start..end(strand)]
    <prefix>.tsv   the same as a table, with the sequence

USAGE:
    python polyprotein_peptides.py NC_001612.1 -o out/EV-A71
    python polyprotein_peptides.py NP_041277.1 -o out/EV-A71_pp     (only that polyprotein)
    python polyprotein_peptides.py --gff3 genome.gff3 --genome-fasta genome.fasta -o out/local
Set NCBI_API_KEY (and optionally NCBI_EMAIL) for higher NCBI rate limits.
'''

##### Imports #####
import argparse
import json
import os
import re
import sys
import time
import urllib.parse
import urllib.request

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
MATURE_TYPES = {"mature_protein_region_of_CDS", "mature_peptide", "mat_peptide", "mature_protein_region"}
COMP = str.maketrans("ACGTUNacgtun", "TGCAANtgcaan")

def revcomp(s):
    return s.translate(COMP)[::-1]

##### NCBI #####
def eutils(tool, **params):
    '''Call an E-utility and return the response text. Retries a few times on errors/429s.'''
    key, email = os.environ.get("NCBI_API_KEY"), os.environ.get("NCBI_EMAIL")
    if key:
        params["api_key"] = key
    if email:
        params["email"] = email
    params["tool"] = "stickleback_polyprotein_peptides"
    data = urllib.parse.urlencode(params).encode()
    for attempt in range(4):
        try:
            with urllib.request.urlopen(EUTILS + tool, data=data, timeout=60) as r:
                return r.read().decode()
        except Exception as e:
            if attempt == 3:
                raise SystemExit(f"NCBI {tool} failed for {params.get('id')}: {e}")
            time.sleep(2 ** (attempt + 1))

def isProteinAccession(acc):
    '''RefSeq protein prefixes (NP_, YP_, XP_, AP_, WP_) or INSDC protein IDs (3 letters + 5 digits,
    e.g. QHD43415). Everything else is treated as nucleotide.'''
    return bool(re.match(r"^([NYXAW]P_\d+|[A-Z]{3}\d{5,7})(\.\d+)?$", acc))

def nuccoreForProtein(acc):
    '''Nucleotide record that codes for a protein accession.'''
    js = json.loads(eutils("elink.fcgi", dbfrom="protein", db="nuccore", id=acc,
                           linkname="protein_nuccore", retmode="json"))
    for ls in js.get("linksets", []):
        for db in ls.get("linksetdbs", []):
            if db.get("links"):
                return db["links"][0]
    raise SystemExit(f"No nucleotide record linked to protein {acc}")

def fetchGff3(nucId):
    return eutils("efetch.fcgi", db="nuccore", id=nucId, rettype="gff3", retmode="text")

def fetchFasta(db, ids):
    return parseFasta(eutils("efetch.fcgi", db=db, id=",".join(ids), rettype="fasta", retmode="text").splitlines())

##### Parsing #####
def parseFasta(lines):
    '''{first word of header: sequence}. Also keyed without the version (NP_041277.1 -> NP_041277).'''
    seqs, name, buf = {}, None, []
    for l in lines:
        l = l.strip()
        if l.startswith(">"):
            if name:
                seqs[name] = "".join(buf)
            name, buf = l[1:].split()[0], []
        elif l:
            buf.append(l)
    if name:
        seqs[name] = "".join(buf)
    for k in list(seqs):
        # NCBI headers can look like ref|NP_041277.1| ; keep the accession part too
        for part in k.split("|"):
            if part:
                seqs.setdefault(part, seqs[k])
                seqs.setdefault(part.split(".")[0], seqs[k])
    return seqs

def parseAttrs(s):
    attrs = {}
    for kv in s.strip().strip(";").split(";"):
        if "=" in kv:
            k, v = kv.split("=", 1)
            attrs[k.strip()] = urllib.parse.unquote(v.strip())
    return attrs

def parseGff3(text):
    '''Returns (cdsById, matureFeatures). Lines sharing an ID (join()s) are merged into one feature
    with an ordered list of segments. Stops at ##FASTA.'''
    feats, order = {}, []
    for n, line in enumerate(text.splitlines()):
        if line.startswith("##FASTA"):
            break
        if not line.strip() or line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 9:
            continue
        seqid, ftype, start, end, strand, phase, attrs = f[0], f[2], int(f[3]), int(f[4]), f[6], f[7], parseAttrs(f[8])
        if ftype != "CDS" and ftype not in MATURE_TYPES:
            continue
        fid = (ftype, attrs.get("ID") or f"line{n}")
        if fid not in feats:
            feats[fid] = {"type": ftype, "seqid": seqid, "strand": strand, "attrs": attrs,
                          "segs": [], "phase": int(phase) if phase.isdigit() else 0}
            order.append(fid)
        feats[fid]["segs"].append((start, end))
    for feat in feats.values():
        # transcript order: ascending on +, descending on - (GFF3 lines may come in either order)
        feat["segs"].sort(key=lambda se: se[0], reverse=(feat["strand"] == "-"))
    cds = {feats[k]["attrs"].get("ID", k[1]): feats[k] for k in order if k[0] == "CDS"}
    mature = [feats[k] for k in order if k[0] != "CDS"]
    return cds, mature

##### Coordinate mapping #####
def positions(feat):
    '''Genomic positions (1-based) of a feature in transcript order.'''
    out = []
    for s, e in feat["segs"]:
        out.extend(range(s, e + 1) if feat["strand"] != "-" else range(e, s - 1, -1))
    return out

def proteinRange(mat, cds):
    '''1-based inclusive amino-acid range of a mature region within its parent CDS's protein.'''
    cdsPos = positions(cds)[cds["phase"]:]
    matPos = positions(mat)
    try:
        i = cdsPos.index(matPos[0])
        j = cdsPos.index(matPos[-1], i)
    except ValueError:
        return None
    return i // 3 + 1, j // 3 + 1

def findParent(mat, cdsById):
    for p in mat["attrs"].get("Parent", "").split(","):
        if p in cdsById:
            return p
    # no usable Parent: smallest CDS on the same sequence/strand containing the region
    lo, hi = min(min(s) for s in mat["segs"]), max(max(s) for s in mat["segs"])
    hits = [(max(e for _, e in c["segs"]) - min(s for s, _ in c["segs"]), cid) for cid, c in cdsById.items()
            if c["seqid"] == mat["seqid"] and c["strand"] == mat["strand"]
            and min(s for s, _ in c["segs"]) <= lo and hi <= max(e for _, e in c["segs"])
            and proteinRange(mat, c)]
    return min(hits)[1] if hits else None

##### Translation (offline fallback) #####
def codonTable(n):
    bases = "TCAG"
    aas = {1:  "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
           2:  "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSS**VVVVAAAADDEEGGGG",
           4:  "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
           5:  "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSSSVVVVAAAADDEEGGGG",
           11: "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"}
    if n not in aas:
        print(f"WARNING: transl_table {n} not built in; using the standard code", file=sys.stderr)
    aa = aas.get(n, aas[1])
    return {a + b + c: aa[16 * i + 4 * j + k] for i, a in enumerate(bases) for j, b in enumerate(bases) for k, c in enumerate(bases)}

def translateCds(cds, genome):
    seq = genome[cds["seqid"]].upper().replace("U", "T")
    nt = "".join(seq[s - 1:e] if cds["strand"] != "-" else revcomp(seq[s - 1:e]) for s, e in cds["segs"])[cds["phase"]:]
    table = codonTable(int(cds["attrs"].get("transl_table", 1)))
    prot = "".join(table.get(nt[i:i + 3], "X") for i in range(0, len(nt) - 2, 3))
    return prot[:-1] if prot.endswith("*") else prot

##### Main #####
def regionString(feat):
    return feat["seqid"] + ":" + ",".join(f"{s}..{e}" for s, e in feat["segs"]) + f"({feat['strand']})"

def main():
    ap = argparse.ArgumentParser(description="Protein sequences of the mature proteins in a viral polyprotein, from its GFF3.")
    ap.add_argument("accession", nargs="?", help="NCBI nucleotide (genome) or protein (polyprotein) accession")
    ap.add_argument("-o", "--out", help="output prefix (writes <prefix>.faa and <prefix>.tsv); default: FASTA to stdout")
    ap.add_argument("--gff3", help="local GFF3 instead of downloading it")
    ap.add_argument("--protein-fasta", help="local polyprotein FASTA (headers = CDS protein_id)")
    ap.add_argument("--genome-fasta", help="local genome FASTA; the polyprotein is translated from the CDS")
    ap.add_argument("--keep-gff3", action="store_true", help="also save the downloaded GFF3 as <prefix>.gff3")
    args = ap.parse_args()
    if not args.accession and not args.gff3:
        ap.error("give an accession or --gff3")

    # 1. GFF3
    wantProtein = args.accession if args.accession and isProteinAccession(args.accession) else None
    if args.gff3:
        with open(args.gff3) as fh:
            gff = fh.read()
    else:
        nuc = nuccoreForProtein(args.accession) if wantProtein else args.accession
        gff = fetchGff3(nuc)
    cdsById, mature = parseGff3(gff)
    if not mature:
        raise SystemExit("No mature protein regions (mature_protein_region_of_CDS / mat_peptide) in the GFF3.")

    # 2. assign each region to its polyprotein and get aa coordinates
    rows = []
    for m in mature:
        pid = findParent(m, cdsById)
        if pid is None:
            print(f"WARNING: no parent CDS for {m['attrs'].get('product', m['attrs'].get('ID'))}; skipped", file=sys.stderr)
            continue
        cds = cdsById[pid]
        polyId = cds["attrs"].get("protein_id") or cds["attrs"].get("Name") or pid
        if wantProtein and polyId.split(".")[0] != wantProtein.split(".")[0]:
            continue
        rng = proteinRange(m, cds)
        if rng is None:
            print(f"WARNING: {m['attrs'].get('product')} does not map onto {polyId}; skipped", file=sys.stderr)
            continue
        rows.append({"poly": polyId, "cds": cds, "mat": m, "aaStart": rng[0], "aaEnd": rng[1]})
    if not rows:
        raise SystemExit(f"No mature protein regions found{' for ' + wantProtein if wantProtein else ''}.")

    # 3. polyprotein sequences
    polyIds = list(dict.fromkeys(r["poly"] for r in rows))
    if args.protein_fasta:
        with open(args.protein_fasta) as fh:
            prots = parseFasta(fh)
    elif args.genome_fasta:
        with open(args.genome_fasta) as fh:
            genome = parseFasta(fh)
        prots = {r["poly"]: translateCds(r["cds"], genome) for r in rows}
    else:
        prots = fetchFasta("protein", polyIds)

    # 4. write
    if args.out:
        os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    faa = open(args.out + ".faa", "w") if args.out else sys.stdout
    tsv = open(args.out + ".tsv", "w") if args.out else None
    if tsv:
        tsv.write("polyprotein_id\tproduct\tmature_protein_id\taa_start\taa_end\tlength\tgenomic_region\tsequence\n")
    for r in rows:
        prot = prots.get(r["poly"]) or prots.get(r["poly"].split(".")[0])
        if prot is None:
            print(f"WARNING: no sequence for polyprotein {r['poly']}; skipped", file=sys.stderr)
            continue
        if r["aaEnd"] > len(prot):
            print(f"WARNING: {r['poly']} is {len(prot)} aa but region ends at {r['aaEnd']}", file=sys.stderr)
        a = r["mat"]["attrs"]
        product = a.get("product") or a.get("Name") or a.get("ID", "")
        matId = a.get("protein_id", "")
        pep = prot[r["aaStart"] - 1:r["aaEnd"]]
        header = f">{r['poly']}|{product.replace(' ', '_')}|{r['aaStart']}-{r['aaEnd']}"
        header += (f" mature_id={matId}" if matId else "") + f" region={regionString(r['mat'])}"
        faa.write(header + "\n" + "\n".join(pep[i:i + 60] for i in range(0, len(pep), 60)) + "\n")
        if tsv:
            tsv.write(f"{r['poly']}\t{product}\t{matId}\t{r['aaStart']}\t{r['aaEnd']}\t{len(pep)}\t{regionString(r['mat'])}\t{pep}\n")
    if args.out:
        faa.close()
        tsv.close()
        if args.keep_gff3 and not args.gff3:
            with open(args.out + ".gff3", "w") as fh:
                fh.write(gff)
        print(f"{len(rows)} mature proteins from {len(polyIds)} polyprotein(s) -> {args.out}.faa, {args.out}.tsv", file=sys.stderr)

if __name__ == "__main__":
    main()
