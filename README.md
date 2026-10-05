# oarfish

<p align="center"><img src="docs/oarfish.jpg" alt="Hand-coloured plate of an oarfish, Gymnetrus gladius" width="720"></p>
<p align="center"><sub><i>Le Gymnètre épée</i> (<i>Gymnetrus gladius</i>), Acarie Baron del., Sebin sc. Cuvier, <i>Le Règne Animal</i>, Poissons pl. 69.</sub></p>

Viral genome annotation utilities from the Quantitative Virology and Evolution Unit.

## Polyprotein → mature proteins: `polyprotein_peptides.py`
Give an NCBI accession for a viral genome (e.g. `NC_001612.1`) or a polyprotein (e.g. `NP_041277.1`). The script downloads the genome's GFF3 from NCBI, finds each mature protein region (`mature_protein_region_of_CDS` / `mat_peptide`), maps it onto its parent polyprotein CDS, and writes the protein sequence of each region. Only the standard library is needed.

```
python polyprotein_peptides.py NC_001612.1 -o out/genome           # every polyprotein in the genome
python polyprotein_peptides.py NP_041277.1 -o out/pp               # just that polyprotein
python polyprotein_peptides.py NC_045512.2                        # FASTA to stdout
```

Outputs:
- `out/X.faa`: one record per mature protein, `>polyprotein_id|product|aaStart-aaEnd mature_id=... region=seqid:start..end(strand)`.
- `out/X.tsv`: `polyprotein_id, product, mature_protein_id, aa_start, aa_end, length, genomic_region, sequence`.

Genomic coordinates are mapped through the spliced CDS, so ribosomal frameshifts (e.g. coronavirus pp1ab/nsp12) map correctly. A mature protein shared by pp1a and pp1ab is listed under each. Polyprotein sequences come from NCBI's protein records. Offline, use `--gff3 file.gff3` with either `--protein-fasta polyproteins.faa` (headers = CDS `protein_id`) or `--genome-fasta genome.fasta`, which translates the CDS instead. `--keep-gff3` saves the downloaded GFF3. Set `NCBI_API_KEY` for higher rate limits.

Test: `python tests/test_polyprotein_peptides.py`

## Web interface (GitHub Pages)
`docs/` is a static site that runs the same steps in the browser. Enter an accession, or upload a GFF3 plus a protein or genome FASTA. It draws a map of each polyprotein and shows a table of mature proteins, with FASTA/TSV/GFF3 downloads. `?acc=NC_045512.2` in the URL runs that accession directly.

To publish: in the repo's **Settings → Pages**, set **Source: Deploy from a branch**, branch `main`, folder `/docs`. To try it locally, open `docs/index.html` in a browser.

The parsing and mapping logic in `docs/polyprotein.js` mirrors `polyprotein_peptides.py`. The test checks that both give identical output (needs `node`).
