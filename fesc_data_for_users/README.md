# NGC 628 escape-fraction data — for users

A small, standalone way to get the matched HII-region/star-cluster catalogue and the
escape-fraction inputs from the NGC 628 LEGUS × SIGNALS pilot study (Teh et al. 2023,
MNRAS 524, 1191). You do **not** need the analysis pipeline, the SLUG2 model libraries,
or any of the multi-GB data — only Python with `numpy` and the two input files below.

## What you get

Run `extract_fesc_data.py` (see below) to produce two CSVs:

- **`matched_pairs.csv`** — one row per (HII region, matched cluster): the ID linkage
  plus region RA/DEC/L_Hα/size/class and cluster RA/DEC/age/mass/E(B−V)/class.
- **`fesc_per_region.csv`** — one row per HII region in the final sample: L_Hα, Q(H⁰)
  with percentiles, and the escape fraction f_esc with percentiles.

If you only want the already-computed values, the CSVs already shipped in this folder
are self-sufficient — you don't need the `.npy` files or to run anything.

## What you need to run it yourself

Put both files in one folder and point `--dat` at it:

| File | Size | What it is |
|------|------|-----------|
| `combined_catalogue.npy` | ~1.4 MB | HII region ↔ star-cluster matching. |
| `GenKroupc1234LHa_QH0_noCons.npy` | ~222 MB | Per-region Monte-Carlo Q(H⁰)/f_esc PDFs (final sample). |

**Where to get them**

- Easiest — request both from the author, **Jia Wei Teh** (`jiaweiteh.astro@gmail.com`).
- If you have the repository, you can rebuild `combined_catalogue.npy` yourself for free
  (it only needs the bundled LEGUS + SIGNALS catalogues). From the repo root:
  ```
  python -c "from src.tools.create_combined_table import create_combined_table as f; f()"
  ```
  This writes `src/dat/combined_catalogue.npy`.
- `GenKroupc1234LHa_QH0_noCons.npy` is built from the SLUG2 libraries (not in the public
  repo), so it can't be regenerated from the repo alone — request it.

## Run

```
python extract_fesc_data.py                    # both .npy in the current folder
python extract_fesc_data.py --dat /path/to/npy # or wherever they live
```

## How the escape fraction is defined

Per HII region (Teh et al. 2023):

```
Q_Hα  = L_Hα × 7.31e11             # ionising-photon rate implied by observed Hα (Kennicutt)
f_esc = (Q(H⁰) − Q_Hα) / Q(H⁰)     # = 1 − Q_Hα / Q(H⁰)
```

`Q(H⁰)` is the **summed** intrinsic ionising-photon budget of all clusters matched to the
region (10⁵ Monte-Carlo realisations of each cluster's cluster_slug Q(H⁰) posterior, added
across clusters). The `fesc_median_recomputed` column applies the formula above to the
stored Q(H⁰) samples and reproduces the pipeline's `fesc_median` exactly — so you can
change the L_Hα→Q_Hα conversion or the f_esc definition and recompute directly.

## Important: use the quality cut

`fesc_per_region.csv` is the **pre-cleaning** sample. Many regions have poorly constrained
Q(H⁰) and carry wild (large-negative) f_esc. The paper keeps only regions whose 1σ
log-Q(H⁰) width is < 0.5 dex — column **`passes_QH0_cut`** (1 = keep). **Filter on
`passes_QH0_cut == 1` before averaging or histogramming f_esc.**

## Column reference — `fesc_per_region.csv`

| Column | Meaning |
|--------|---------|
| `h2_ID` | SIGNALS HII-region ID |
| `n_clusters`, `matched_cluster_IDs` | number and LEGUS IDs of matched clusters |
| `logLHa`, `LHa_erg_s` | extinction-corrected Hα luminosity |
| `logQH0_{2.3,15.9,50,84.1,97.7}pct` | log Q(H⁰) percentiles (2σ, 1σ, median) |
| `logQH0_1sigma_width` | 84.1pct − 15.9pct, the quantity the cut is on |
| `fesc_{15.9,median,84.1}pct` | f_esc percentiles from the pipeline |
| `fesc_median_recomputed` | f_esc median recomputed here from Q(H⁰) & L_Hα |
| `log_clusterMass_median_Msun` | log of the region's summed cluster mass |
| `passes_QH0_cut` | 1 if `logQH0_1sigma_width < 0.5` (paper's selection) |

## Column reference — `matched_pairs.csv`

`h2_ID, h2_RA, h2_DEC, h2_rGalac_kpc, h2_LHa_erg_s, h2_logLHa, h2_size_pc, h2_class,
sc_ID, sc_RA, sc_DEC, separation_arcsec, sc_age_yr, sc_mass_Msun, sc_EBV, sc_class`

This table is the spatial matching only; Q(H⁰) is a per-*region* summed quantity (see
`fesc_per_region.csv`), not a per-cluster column here.

---
Contact: Jia Wei Teh · jiaweiteh.astro@gmail.com
