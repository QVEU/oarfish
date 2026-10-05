/* UI for docs/index.html. Core logic lives in polyprotein.js (window.Polyprotein). */
(function () {
  "use strict";
  const P = window.Polyprotein;
  const $ = id => document.getElementById(id);
  const COLORS = ["#2a7f8a", "#c0582b", "#5a6fc4", "#b8860b", "#8b4f9e", "#3d8b3d", "#c4436b", "#6b7f1f", "#2f6db0", "#a0522d", "#4f8f7a", "#9b5f2a"];

  let state = { records: [], gff3: "", label: "", filter: null };

  // ---------- helpers ----------
  const el = (tag, attrs, ...kids) => {
    const n = document.createElement(tag);
    for (const [k, v] of Object.entries(attrs || {})) {
      if (k === "class") n.className = v;
      else if (k.startsWith("on")) n.addEventListener(k.slice(2), v);
      else n.setAttribute(k, v);
    }
    for (const k of kids) if (k !== null && k !== undefined) n.append(k);
    return n;
  };
  const readFile = input => (input.files[0] ? input.files[0].text() : Promise.resolve(null));
  const safeName = s => (s || "polyprotein").replace(/[^\w.-]+/g, "_");

  function setStatus(msg, isError) {
    const s = $("status");
    s.hidden = !msg;
    s.textContent = msg || "";
    s.classList.toggle("error", !!isError);
  }
  function setWarnings(list) {
    const w = $("warnings");
    w.replaceChildren(...list.map(t => el("li", {}, "⚠ " + t)));
    w.hidden = !list.length;
  }
  function download(name, text, type) {
    const a = el("a", { href: URL.createObjectURL(new Blob([text], { type })), download: name });
    document.body.append(a);
    a.click();
    setTimeout(() => { URL.revokeObjectURL(a.href); a.remove(); }, 0);
  }
  async function copy(text, btn) {
    try {
      await navigator.clipboard.writeText(text);
      const t = btn.textContent;
      btn.textContent = "Copied";
      setTimeout(() => (btn.textContent = t), 1200);
    } catch (e) {
      setStatus("Could not copy to the clipboard in this browser.", true);
    }
  }
  function busy(on) {
    document.querySelectorAll("button.primary").forEach(b => (b.disabled = on));
  }

  // ---------- tabs ----------
  function selectTab(which) {
    for (const t of ["ncbi", "local"]) {
      $("tab-" + t).setAttribute("aria-selected", String(t === which));
      $("panel-" + t).hidden = t !== which;
    }
  }
  $("tab-ncbi").addEventListener("click", () => selectTab("ncbi"));
  $("tab-local").addEventListener("click", () => selectTab("local"));

  // ---------- API key ----------
  try { $("apikey").value = localStorage.getItem("ncbiApiKey") || ""; } catch (e) { /* storage unavailable */ }
  $("apikey").addEventListener("change", () => {
    try { localStorage.setItem("ncbiApiKey", $("apikey").value.trim()); } catch (e) { /* ignore */ }
  });

  // ---------- run ----------
  async function go(opts, label) {
    const warnings = [];
    busy(true);
    setWarnings([]);
    $("results").hidden = true;
    setStatus("Working…");
    try {
      const res = await P.run(Object.assign({ apiKey: $("apikey").value.trim(), warn: w => warnings.push(w), progress: m => setStatus(m) }, opts));
      state = { records: res.records, gff3: res.gff3, label, filter: null };
      setStatus("");
      render();
    } catch (e) {
      const cors = e instanceof TypeError || /Failed to fetch|NetworkError/i.test(e.message);
      setStatus(cors ? `Could not reach NCBI (${e.message}). Check your connection, or use the Local files tab.` : e.message, true);
    } finally {
      setWarnings(warnings);
      busy(false);
    }
  }

  $("panel-ncbi").addEventListener("submit", ev => {
    ev.preventDefault();
    const acc = $("acc").value.trim();
    if (!acc) return;
    const u = new URL(location.href);
    u.searchParams.set("acc", acc);
    history.replaceState(null, "", u);
    go({ accession: acc }, acc);
  });
  document.querySelectorAll(".examples button").forEach(b =>
    b.addEventListener("click", () => {
      $("acc").value = b.dataset.acc;
      $("panel-ncbi").requestSubmit();
    }));

  $("panel-local").addEventListener("submit", async ev => {
    ev.preventDefault();
    const [gff, prot, genome] = await Promise.all([readFile($("gffFile")), readFile($("protFile")), readFile($("genomeFile"))]);
    if (!gff) { setStatus("Choose a GFF3 file.", true); return; }
    const acc = $("localAcc").value.trim();
    const label = ($("gffFile").files[0].name || "local").replace(/\.(gff3?|txt)$/i, "");
    go({ accession: acc, gff3Text: gff, proteinFastaText: prot, genomeFastaText: prot ? null : genome }, acc || label);
  });

  // ---------- render ----------
  function colorFor(polyRecords) {
    const m = new Map();
    polyRecords.forEach((r, i) => m.set(r, COLORS[i % COLORS.length]));
    return m;
  }

  function highlight(key, on) {
    document.querySelectorAll(`[data-key="${key}"]`).forEach(n => n.classList.toggle("hl", on));
  }

  function render() {
    const recs = state.records;
    const polys = [...new Set(recs.map(r => r.poly))];
    const shown = state.filter ? recs.filter(r => r.poly === state.filter) : recs;
    const colors = new Map();
    for (const p of polys) for (const [r, c] of colorFor(recs.filter(r => r.poly === p))) colors.set(r, c);

    $("resTitle").textContent = `Mature proteins · ${state.label}`;
    $("resSummary").textContent = `${recs.length} mature protein${recs.length === 1 ? "" : "s"} from ${polys.length} polyprotein${polys.length === 1 ? "" : "s"}`;
    $("dlGff").hidden = !state.gff3;

    // filter chips
    const chips = $("polyFilter");
    chips.replaceChildren();
    if (polys.length > 1) {
      const mk = (val, text) => el("button", {
        type: "button", "aria-pressed": String(state.filter === val),
        onclick: () => { state.filter = val; render(); },
      }, text);
      chips.append(mk(null, "All"), ...polys.map(p => mk(p, p)));
    }

    // maps
    const maps = $("maps");
    maps.replaceChildren();
    for (const p of state.filter ? [state.filter] : polys) {
      const pr = recs.filter(r => r.poly === p);
      const L = pr[0].polyLength || Math.max(...pr.map(r => r.aaEnd));
      const track = el("div", { class: "track" });
      pr.forEach((r, i) => {
        const key = recs.indexOf(r);
        const seg = el("div", {
          class: "seg", "data-key": key,
          title: `${r.product}\n${r.aaStart}–${r.aaEnd} (${r.length} aa)`,
          style: `left:${((r.aaStart - 1) / L) * 100}%;width:${(r.length / L) * 100}%;background:${colors.get(r)}`,
          onmouseenter: () => highlight(key, true),
          onmouseleave: () => highlight(key, false),
          onclick: () => document.querySelector(`tr[data-key="${key}"]`).scrollIntoView({ behavior: "smooth", block: "center" }),
        }, r.product);
        track.append(seg);
      });
      maps.append(el("div", { class: "map" },
        el("div", { class: "map-head" },
          el("span", {}, el("strong", {}, p), ` · ${pr.length} region${pr.length === 1 ? "" : "s"}`),
          el("a", { href: `https://www.ncbi.nlm.nih.gov/protein/${encodeURIComponent(p)}`, target: "_blank", rel: "noopener" }, "NCBI ↗")),
        track,
        el("div", { class: "scale" }, el("span", {}, "1"), el("span", {}, `${L} aa`))));
    }

    // table
    const tbody = $("table").querySelector("tbody");
    tbody.replaceChildren(...shown.map(r => {
      const key = recs.indexOf(r);
      const seq = el("span", { class: "seq", title: "Click to expand" }, r.sequence);
      seq.addEventListener("click", () => seq.classList.toggle("open"));
      const fasta = P.toFasta([r]);
      return el("tr", {
        "data-key": key,
        onmouseenter: () => highlight(key, true),
        onmouseleave: () => highlight(key, false),
      },
        el("td", {}, el("span", { class: "swatch", style: `background:${colors.get(r)}` }), r.product),
        el("td", { class: "mono" }, r.poly),
        el("td", { class: "num" }, `${r.aaStart}–${r.aaEnd}`),
        el("td", { class: "num" }, String(r.length)),
        el("td", { class: "mono" }, r.matureId
          ? el("a", { href: `https://www.ncbi.nlm.nih.gov/protein/${encodeURIComponent(r.matureId)}`, target: "_blank", rel: "noopener" }, r.matureId)
          : "—"),
        el("td", { class: "mono" }, r.region),
        el("td", {}, el("div", { class: "seqcell" }, seq,
          el("button", { type: "button", title: "Copy as FASTA", onclick: ev => copy(fasta, ev.currentTarget) }, "Copy"))));
    }));
    $("results").hidden = false;
  }

  const current = () => (state.filter ? state.records.filter(r => r.poly === state.filter) : state.records);
  const base = () => safeName(state.filter || state.label);
  $("dlFasta").addEventListener("click", () => download(base() + "_mature.faa", P.toFasta(current()), "text/plain"));
  $("dlTsv").addEventListener("click", () => download(base() + "_mature.tsv", P.toTsv(current()), "text/tab-separated-values"));
  $("dlGff").addEventListener("click", () => download(safeName(state.label) + ".gff3", state.gff3, "text/plain"));
  $("copyFasta").addEventListener("click", ev => copy(P.toFasta(current()), ev.currentTarget));

  // ---------- deep link ?acc= ----------
  const q = new URLSearchParams(location.search).get("acc");
  if (q) {
    $("acc").value = q;
    go({ accession: q }, q);
  }
})();
