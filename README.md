# LEGUS-SIGNALS-NGC628

Analysis code and figure-generating scripts for:

> **Teh et al. (2023)**, "Constraining the LyC escape fraction from LEGUS star clusters
> with SIGNALS H II region observations: a pilot study of NGC 628", *MNRAS*, **524**, 1191.
>
> [ADS](https://ui.adsabs.harvard.edu/abs/2023MNRAS.524.1191T/abstract) ·
> [arXiv:2306.05457](https://arxiv.org/abs/2306.05457) ·
> [doi:10.1093/mnras/stad1780](https://doi.org/10.1093/mnras/stad1780)

<p><h2>Repo Structure</h2></p>
<ul>
    <li><code>./lib:</code> Input catalogues and model libraries. The small catalogues are tracked in Git; the large SLUG2 model library is kept locally (see <a href="DATA.md"><code>DATA.md</code></a>).
    <ul class="square">
          <li><code>/LEGUS:</code> LEGUS NGC 628 star-cluster catalogue. </li>
          <li><code>/SIGNALS:</code> SIGNALS Hα H II-region catalogue. </li>
          <li><code>/SLUG2:</code> SLUG2 model libraries (<code>cluster_slug</code> Bayesian PDF tables, <code>cluster_lib</code> cluster library) — git-ignored. </li>
        </ul></li>
    <li><code>./src:</code> Source code for the paper.
        <ul class="square">
          <li><code>paths.py:</code> Central path constants, all resolved relative to the repo root.</li>
          <li><code>/tools:</code> Catalogue builders (<code>create_*_table</code>) and analysis helpers (geometry, stats, plotting).</li>
          <li><code>/reg_files:</code> Region and selection definitions (FOV, classes, associations, H II-region masks).</li>
          <li><code>/dat:</code> Post-processed arrays (<code>.npy</code>) the figures read; not versioned (see <a href="#data">Data</a>).</li>
          <li><code>/plots:</code> Scripts that produce the paper figures (<code>fg*</code>), with <code>alternative_params/</code>, <code>dark/</code>, and <code>extras/</code> variants.</li>
          <li><code>/figs:</code> Generated figures; output from the <code>/plots</code> scripts.</li>
        </ul></li>
    <li><code>make_figures.py:</code> Master runner at the repo root — regenerates <strong>every</strong> paper figure by running each <code>src/plots/fg*.py</code> in turn. Start here (see <a href="#setup--run">Setup &amp; run</a>).</li>
    <li><code>fesc_handoff/:</code> Standalone, <code>numpy</code>-only bundle that extracts the matched H II-region × star-cluster catalogue and the escape-fraction inputs to plain CSVs — for sharing with collaborators who don't need the full pipeline (see <a href="#extracting-h-ii-region--cluster-data-escape-fraction">Extracting H II region × cluster data</a>).</li>
</ul>

## Setup & run

Install dependencies (Python 3.10+):

```
python -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
```

(Optional: a system LaTeX install is needed for `usetex` figure rendering.)

### Reproduce all figures (start here)

The quickest way to regenerate every paper figure is the master runner at the repo root:

```
python make_figures.py            # run all figure scripts -> src/figs/
python make_figures.py fg5 fg7    # only scripts whose name matches
python make_figures.py --list     # list what would run, then exit
```

It runs each `src/plots/fg*.py` in its own subprocess, continues past any failure, and
prints a summary at the end. This needs the post-processed arrays under `src/dat/` to be
in place (see [Data](#data) below); `fg10_clusterstatistics.py` in particular regenerates
heavy Monte-Carlo arrays and is slow.

### Run a single figure

Each `src/plots/fg*.py` script reproduces one figure from the paper and writes it to
`src/figs/`. Run them **from the repository root**, e.g.:

```
python src/plots/fg5_LHa_vs_QH0.py
```

Each script adds the repo root to `sys.path` at startup, so no `PYTHONPATH` or `-m` is
needed; running inside Spyder (with the folder as the project root) also works.

The scripts run as-is and regenerate the paper figures **once the post-processed data
arrays are in place under `src/dat/`** (see [Data](#data) below).

## Data

The figure scripts read **post-processed data arrays** (`.npy`) from `src/dat/`. The bulk
of these are the Monte-Carlo Q(H⁰) distributions (`*LHa_QH0_noCons*.npy`, ~850 MB) and the
clusterslug PDF tables (~180 MB). They are too large to track in Git and are **not
included** in this repository.

Only part of this can be rebuilt from what *is* bundled in Git:

- `combined_catalogue.npy` is rebuilt by `src/tools/create_combined_table` from the bundled
  LEGUS and SIGNALS catalogues — no extra data required.
- The clusterslug tables and the large `*LHa_QH0_noCons*.npy` Monte-Carlo arrays require the
  **raw SLUG2 model libraries** (`lib/SLUG2/`), which are also not tracked in Git. Note the
  heavy arrays are produced by the *plot scripts themselves* (e.g. `fg5_LHa_vs_QH0.py`, via
  the `reRun` flag near the top), not by `src/tools/`.

So a fresh clone reproduces the figures only once `src/dat/` is populated. Two ways to get
there:

- **Easiest — request the data.** The full `src/dat/` bundle (and the SLUG2 libraries) can be
  shared on request — contact **Jia Wei Teh** at **jiaweiteh.astro@gmail.com**. Drop the
  arrays into `src/dat/` and every `fg*` script reproduces its figure as-is.
- **Full regeneration.** Obtain the SLUG2 libraries and place them under `lib/SLUG2/`, build
  the catalogue and clusterslug tables in `src/tools/`, then run the `fg*` scripts with their
  `reRun` flag set to `True` to regenerate the Monte-Carlo arrays.

See [`DATA.md`](DATA.md) for a full inventory of what lives where. Paths are configured
centrally in `src/paths.py`.

## Extracting H II region × cluster data (escape fraction)

`fesc_handoff/` is a standalone bundle for handing the matched catalogue and the
escape-fraction inputs to a collaborator who does **not** need the full pipeline, the
SLUG2 libraries, or any of the multi-GB data. It depends only on `numpy` and does not
import the `src` package.

It reads two post-processed arrays from `src/dat/`:

- `combined_catalogue.npy` (~1.4 MB) — every H II region matched to the star cluster(s)
  inside its radius, with both catalogues' columns.
- `GenKroupc1234LHa_QH0_noCons.npy` (~222 MB) — per-region Monte-Carlo Q(H⁰) and f_esc
  PDFs for the final sample (Geneva tracks, Kroupa IMF, classes 1–4).

Run it from the repo root:

```
cd fesc_handoff
python extract_fesc_data.py --dat ../src/dat --out .
```

This writes two CSVs:

- **`matched_pairs.csv`** — one row per (H II region, matched cluster): the ID linkage
  plus region RA/DEC/L_Hα/size/class and cluster RA/DEC/age/mass/E(B−V)/class
  (1031 pairs across 493 regions / 897 clusters).
- **`fesc_per_region.csv`** — one row per H II region in the final f_esc sample (139
  regions): L_Hα, Q(H⁰) percentiles, and f_esc percentiles.

The escape fraction is defined per region (Teh et al. 2023):

```
Q_Hα  = L_Hα × 7.31e11             # ionising-photon rate implied by observed Hα (Kennicutt)
f_esc = (Q(H⁰) − Q_Hα) / Q(H⁰)     # = 1 − Q_Hα / Q(H⁰)
```

where Q(H⁰) is the **summed** intrinsic ionising-photon budget of the region's matched
clusters (10⁵ Monte-Carlo realisations of each cluster's cluster_slug Q(H⁰) posterior,
added across clusters). The `fesc_median_recomputed` column applies this formula to the
stored Q(H⁰) samples and reproduces the pipeline's `fesc_median` bit-for-bit, so the
L_Hα→Q_Hα conversion or the f_esc definition can be changed and recomputed directly.

**Quality cut.** The 139-region table is the *pre-cleaning* sample; poorly constrained
regions carry large-negative f_esc. The paper (`fg5`/`fg7`) keeps only regions whose 1σ
log-Q(H⁰) width is < 0.5 dex — column **`passes_QH0_cut`** (46 pass). Filter on
`passes_QH0_cut == 1` before averaging or histogramming f_esc.

If a collaborator only wants the already-computed values, the two CSVs are
self-sufficient (no `.npy` needed). See [`fesc_handoff/README.md`](fesc_handoff/README.md)
for the full column dictionary.
