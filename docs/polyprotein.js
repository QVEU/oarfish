/*
 * polyprotein.js — browser port of polyprotein_peptides.py
 * Splits a viral polyprotein into its mature proteins using the GFF3 annotation.
 * Pure functions (parsing, coordinate mapping, translation, output) plus NCBI E-utilities fetchers.
 * Works in the browser (window.Polyprotein) and in Node (module.exports) for testing.
 */
(function (root) {
  "use strict";

  const EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/";
  const MATURE_TYPES = new Set(["mature_protein_region_of_CDS", "mature_peptide", "mat_peptide", "mature_protein_region"]);
  const COMP = { A: "T", C: "G", G: "C", T: "A", U: "A", N: "N", a: "t", c: "g", g: "c", t: "a", u: "a", n: "n" };

  function revcomp(s) {
    let out = "";
    for (let i = s.length - 1; i >= 0; i--) out += COMP[s[i]] || s[i];
    return out;
  }

  // ---------- NCBI ----------
  function isProteinAccession(acc) {
    // RefSeq protein prefixes (NP_, YP_, XP_, AP_, WP_) or INSDC protein IDs (3 letters + 5-7 digits)
    return /^([NYXAW]P_\d+|[A-Z]{3}\d{5,7})(\.\d+)?$/.test(acc);
  }

  async function eutils(tool, params, apiKey) {
    const p = new URLSearchParams(Object.assign({ tool: "oarfish_polyprotein_peptides" }, params));
    if (apiKey) p.set("api_key", apiKey);
    let lastErr;
    for (let attempt = 0; attempt < 4; attempt++) {
      try {
        const r = await fetch(EUTILS + tool + "?" + p.toString());
        if (r.ok) return await r.text();
        lastErr = new Error(`NCBI ${tool} returned HTTP ${r.status}`);
        if (r.status !== 429 && r.status < 500) break;
      } catch (e) {
        lastErr = e;
      }
      await new Promise(res => setTimeout(res, 1000 * 2 ** attempt));
    }
    throw new Error(`NCBI ${tool} failed for ${params.id}: ${lastErr && lastErr.message}`);
  }

  async function nuccoreForProtein(acc, apiKey) {
    const js = JSON.parse(await eutils("elink.fcgi", { dbfrom: "protein", db: "nuccore", id: acc, linkname: "protein_nuccore", retmode: "json" }, apiKey));
    for (const ls of js.linksets || []) for (const db of ls.linksetdbs || []) if (db.links && db.links.length) return db.links[0];
    throw new Error(`No nucleotide record linked to protein ${acc}`);
  }

  const fetchGff3 = (nucId, apiKey) => eutils("efetch.fcgi", { db: "nuccore", id: nucId, rettype: "gff3", retmode: "text" }, apiKey);
  const fetchProteinFasta = async (ids, apiKey) =>
    parseFasta(await eutils("efetch.fcgi", { db: "protein", id: ids.join(","), rettype: "fasta", retmode: "text" }, apiKey));

  // ---------- Parsing ----------
  function parseFasta(text) {
    // {first word of header: sequence}; also keyed by each |-separated part and without version
    const seqs = {};
    let name = null, buf = [];
    for (let l of text.split(/\r?\n/)) {
      l = l.trim();
      if (l.startsWith(">")) {
        if (name) seqs[name] = buf.join("");
        name = l.slice(1).split(/\s+/)[0];
        buf = [];
      } else if (l) buf.push(l);
    }
    if (name) seqs[name] = buf.join("");
    for (const k of Object.keys(seqs)) {
      for (const part of k.split("|")) {
        if (!part) continue;
        if (!(part in seqs)) seqs[part] = seqs[k];
        const nov = part.split(".")[0];
        if (!(nov in seqs)) seqs[nov] = seqs[k];
      }
    }
    return seqs;
  }

  function parseAttrs(s) {
    const attrs = {};
    for (const kv of s.trim().replace(/^;+|;+$/g, "").split(";")) {
      const i = kv.indexOf("=");
      if (i < 0) continue;
      let v = kv.slice(i + 1).trim();
      try { v = decodeURIComponent(v); } catch (e) { /* leave as is */ }
      attrs[kv.slice(0, i).trim()] = v;
    }
    return attrs;
  }

  function parseGff3(text) {
    // Returns {cdsById, mature}. Lines sharing an ID (join()s) are merged into one feature with ordered segments.
    const feats = new Map();
    const lines = text.split(/\r?\n/);
    for (let n = 0; n < lines.length; n++) {
      const line = lines[n];
      if (line.startsWith("##FASTA")) break;
      if (!line.trim() || line.startsWith("#")) continue;
      const f = line.split("\t");
      if (f.length < 9) continue;
      const ftype = f[2];
      if (ftype !== "CDS" && !MATURE_TYPES.has(ftype)) continue;
      const attrs = parseAttrs(f[8]);
      const key = ftype + "\u0000" + (attrs.ID || "line" + n);
      if (!feats.has(key)) {
        feats.set(key, { type: ftype, seqid: f[0], strand: f[6], attrs, segs: [], phase: /^\d+$/.test(f[7]) ? parseInt(f[7], 10) : 0 });
      }
      feats.get(key).segs.push([parseInt(f[3], 10), parseInt(f[4], 10)]);
    }
    const cdsById = {}, mature = [];
    for (const [key, feat] of feats) {
      feat.segs.sort((a, b) => (feat.strand === "-" ? b[0] - a[0] : a[0] - b[0]));
      if (feat.type === "CDS") cdsById[feat.attrs.ID || key.split("\u0000")[1]] = feat;
      else mature.push(feat);
    }
    return { cdsById, mature };
  }

  // ---------- Coordinate mapping ----------
  function positions(feat) {
    const out = [];
    for (const [s, e] of feat.segs) {
      if (feat.strand !== "-") for (let p = s; p <= e; p++) out.push(p);
      else for (let p = e; p >= s; p--) out.push(p);
    }
    return out;
  }

  function proteinRange(mat, cds) {
    const cdsPos = positions(cds).slice(cds.phase);
    const matPos = positions(mat);
    const i = cdsPos.indexOf(matPos[0]);
    if (i < 0) return null;
    const j = cdsPos.indexOf(matPos[matPos.length - 1], i);
    if (j < 0) return null;
    return [Math.floor(i / 3) + 1, Math.floor(j / 3) + 1];
  }

  const span = f => [Math.min(...f.segs.map(s => s[0])), Math.max(...f.segs.map(s => s[1]))];

  function findParent(mat, cdsById) {
    for (const p of (mat.attrs.Parent || "").split(",")) if (p in cdsById) return p;
    // no usable Parent: smallest CDS on the same sequence/strand containing the region
    const [lo, hi] = span(mat);
    let best = null;
    for (const [cid, c] of Object.entries(cdsById)) {
      const [cs, ce] = span(c);
      if (c.seqid !== mat.seqid || c.strand !== mat.strand || cs > lo || hi > ce || !proteinRange(mat, c)) continue;
      if (!best || ce - cs < best[0] || (ce - cs === best[0] && cid < best[1])) best = [ce - cs, cid];
    }
    return best ? best[1] : null;
  }

  // ---------- Translation (offline fallback) ----------
  const TABLES = {
    1: "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    2: "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSS**VVVVAAAADDEEGGGG",
    4: "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    5: "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSSSVVVVAAAADDEEGGGG",
    11: "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
  };

  function codonTable(n, warn) {
    if (!(n in TABLES) && warn) warn(`transl_table ${n} not built in; using the standard code`);
    const aa = TABLES[n] || TABLES[1], b = "TCAG", t = {};
    for (let i = 0; i < 4; i++) for (let j = 0; j < 4; j++) for (let k = 0; k < 4; k++) t[b[i] + b[j] + b[k]] = aa[16 * i + 4 * j + k];
    return t;
  }

  function translateCds(cds, genome, warn) {
    const seq = (genome[cds.seqid] || genome[cds.seqid.split(".")[0]] || "").toUpperCase().replace(/U/g, "T");
    if (!seq) throw new Error(`Sequence ${cds.seqid} not found in the genome FASTA`);
    const nt = cds.segs.map(([s, e]) => (cds.strand !== "-" ? seq.slice(s - 1, e) : revcomp(seq.slice(s - 1, e)))).join("").slice(cds.phase);
    const table = codonTable(parseInt(cds.attrs.transl_table || "1", 10), warn);
    let prot = "";
    for (let i = 0; i + 2 < nt.length; i += 3) prot += table[nt.slice(i, i + 3)] || "X";
    return prot.endsWith("*") ? prot.slice(0, -1) : prot;
  }

  // ---------- Assembly ----------
  const regionString = f => `${f.seqid}:${f.segs.map(([s, e]) => `${s}..${e}`).join(",")}(${f.strand})`;

  function matureRegions(gffText, wantProtein, warn) {
    // -> rows [{poly, cds, mat, aaStart, aaEnd}]
    warn = warn || (() => {});
    const { cdsById, mature } = parseGff3(gffText);
    if (!mature.length) throw new Error("No mature protein regions (mature_protein_region_of_CDS / mat_peptide) in the GFF3.");
    const rows = [];
    for (const m of mature) {
      const label = m.attrs.product || m.attrs.ID;
      const pid = findParent(m, cdsById);
      if (pid === null) { warn(`no parent CDS for ${label}; skipped`); continue; }
      const cds = cdsById[pid];
      const poly = cds.attrs.protein_id || cds.attrs.Name || pid;
      if (wantProtein && poly.split(".")[0] !== wantProtein.split(".")[0]) continue;
      const rng = proteinRange(m, cds);
      if (!rng) { warn(`${label} does not map onto ${poly}; skipped`); continue; }
      rows.push({ poly, cds, mat: m, aaStart: rng[0], aaEnd: rng[1] });
    }
    if (!rows.length) throw new Error(`No mature protein regions found${wantProtein ? " for " + wantProtein : ""}.`);
    return rows;
  }

  function attachSequences(rows, prots, warn) {
    // -> records [{poly, product, matureId, aaStart, aaEnd, length, region, sequence, polyLength}]
    warn = warn || (() => {});
    const out = [];
    for (const r of rows) {
      const prot = prots[r.poly] || prots[r.poly.split(".")[0]];
      if (prot === undefined) { warn(`no sequence for polyprotein ${r.poly}; skipped`); continue; }
      if (r.aaEnd > prot.length) warn(`${r.poly} is ${prot.length} aa but region ends at ${r.aaEnd}`);
      const a = r.mat.attrs;
      const seq = prot.slice(r.aaStart - 1, r.aaEnd);
      out.push({
        poly: r.poly, product: a.product || a.Name || a.ID || "", matureId: a.protein_id || "",
        aaStart: r.aaStart, aaEnd: r.aaEnd, length: seq.length, region: regionString(r.mat), sequence: seq, polyLength: prot.length,
      });
    }
    return out;
  }

  function toFasta(records) {
    return records.map(r => {
      let h = `>${r.poly}|${r.product.replace(/ /g, "_")}|${r.aaStart}-${r.aaEnd}`;
      if (r.matureId) h += ` mature_id=${r.matureId}`;
      h += ` region=${r.region}`;
      const lines = [];
      for (let i = 0; i < r.sequence.length; i += 60) lines.push(r.sequence.slice(i, i + 60));
      return h + "\n" + lines.join("\n") + "\n";
    }).join("");
  }

  function toTsv(records) {
    const head = "polyprotein_id\tproduct\tmature_protein_id\taa_start\taa_end\tlength\tgenomic_region\tsequence\n";
    return head + records.map(r => [r.poly, r.product, r.matureId, r.aaStart, r.aaEnd, r.length, r.region, r.sequence].join("\t") + "\n").join("");
  }

  // ---------- End to end ----------
  async function run(opts) {
    // opts: {accession, gff3Text, proteinFastaText, genomeFastaText, apiKey, warn, progress}
    const warn = opts.warn || (() => {}), progress = opts.progress || (() => {});
    const acc = (opts.accession || "").trim();
    const wantProtein = acc && isProteinAccession(acc) ? acc : null;
    let gff = opts.gff3Text;
    if (!gff) {
      if (!acc) throw new Error("Give an accession or a GFF3 file.");
      let nuc = acc;
      if (wantProtein) { progress(`Finding the genome that codes for ${acc}…`); nuc = await nuccoreForProtein(acc, opts.apiKey); }
      progress(`Downloading GFF3 for ${wantProtein ? "nuccore " + nuc : nuc}…`);
      gff = await fetchGff3(nuc, opts.apiKey);
      if (/^\s*(Error|<\?xml|<ERROR)/i.test(gff)) throw new Error(`NCBI did not return a GFF3 for ${acc}: ${gff.slice(0, 200)}`);
    }
    const rows = matureRegions(gff, wantProtein, warn);
    const polyIds = [...new Set(rows.map(r => r.poly))];
    let prots;
    if (opts.proteinFastaText) prots = parseFasta(opts.proteinFastaText);
    else if (opts.genomeFastaText) {
      const genome = parseFasta(opts.genomeFastaText);
      prots = {};
      for (const r of rows) if (!(r.poly in prots)) prots[r.poly] = translateCds(r.cds, genome, warn);
    } else {
      progress(`Downloading ${polyIds.length} polyprotein sequence(s)…`);
      prots = await fetchProteinFasta(polyIds, opts.apiKey);
    }
    return { records: attachSequences(rows, prots, warn), gff3: gff };
  }

  const api = {
    revcomp, isProteinAccession, parseFasta, parseAttrs, parseGff3, positions, proteinRange, findParent,
    codonTable, translateCds, regionString, matureRegions, attachSequences, toFasta, toTsv, run,
    nuccoreForProtein, fetchGff3, fetchProteinFasta,
  };
  if (typeof module !== "undefined" && module.exports) module.exports = api;
  else root.Polyprotein = api;
})(typeof window !== "undefined" ? window : globalThis);
