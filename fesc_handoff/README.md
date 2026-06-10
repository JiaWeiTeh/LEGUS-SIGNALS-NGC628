# NGC 628 — HII region × star-cluster data for escape-fraction work

A self-contained bundle to extract the matched HII-region/star-cluster catalogue and
the per-region escape-fraction inputs from the NGC 628 LEGUS × SIGNALS pipeline
(Teh et al. 2023, MNRAS 524, 1191). You do **not** need the full repository, the SLUG2
model libraries, or any of the multi-GB data to use this.

## What you need

Two pipeline arrays (ask Jia Wei for them — together ~224 MB):

| File | Size | What it is |
|------|------|-----------|
| `combined_catalogue.npy` | ~1.4 MB | Every HII region matched to the star cluster(s) inside its radius, with both catalogues' columns. |
| `GenKroupc1234LHa_QH0_noCons.npy` | ~222 MB | Per-region Monte-Carlo Q(H⁰) and f_esc PDFs for the final sample (Geneva tracks, Kroupa IMF, classes 1–4). |

Plus `extract_fesc_data.py` (this bundle) and Python 3 with `numpy`.

If you only want the **already-computed** f_esc values and don't need to re-sample the
PDFs, the small CSVs in this bundle (`matched_pairs.csv`, `fesc_per_region.csv`) may be
all you need — no `.npy` required.

## Run

```
python extract_fesc_data.py --dat /folder/with/the/npy/files --out .
```

Produces:

- **`matched_pairs.csv`** — one row per (HII region, matched cluster). The ID linkage
  plus the key associated values: region RA/DEC/L_Hα/size/class and cluster
  RA/DEC/age/mass/E(B−V)/class. 1031 pairs across 493 regions / 897 clusters.
- **`fesc_per_region.csv`** — one row per HII region in the final f_esc sample (139
  regions). L_Hα, Q(H⁰) percentiles, f_esc percentiles, and a recomputed f_esc.

## How the escape fraction is defined

Per HII region (Teh et al. 2023):

```
Q_Hα  = L_Hα × 7.31e11             # ionising-photon rate implied by observed Hα (Kennicutt)
f_esc = (Q(H⁰) − Q_Hα) / Q(H⁰)     # = 1 − Q_Hα / Q(H⁰)
```

`Q(H⁰)` is the **summed** intrinsic ionising-photon budget of all clusters matched to
the region, built by Monte-Carlo sampling each cluster's cluster_slug Q(H⁰) posterior
(10⁵ realisations) and adding across clusters. `f_esc_median_recomputed` in the CSV is
exactly this formula applied to the stored Q(H⁰) samples — it reproduces the pipeline's
`fesc_median` bit-for-bit, so you can change the L_Hα→Q_Hα conversion or the f_esc
definition and recompute directly.

## Important: use the quality cut

The 139-region table is the **pre-cleaning** sample. Many regions have poorly constrained
Q(H⁰) and carry wild (large-negative) f_esc. The paper keeps only regions whose 1σ
log-Q(H⁰) width is < 0.5 dex — column **`passes_QH0_cut`** (1 = keep). 46 regions pass;
the paper drops a few more by manual inspection. **Filter on `passes_QH0_cut == 1`
before averaging or histogramming f_esc.** Over the passing regions the median f_esc ≈ 0.1.

## `fesc_per_region.csv` columns

| Column | Meaning |
|--------|---------|
| `h2_ID` | SIGNALS HII-region ID |
| `n_clusters`, `matched_cluster_IDs` | number and LEGUS IDs of matched clusters |
| `logLHa`, `LHa_erg_s` | extinction-corrected Hα luminosity |
| `logQH0_{2.3,15.9,50,84.1,97.7}pct` | log Q(H⁰) percentiles (2σ, 1σ, median) |
| `logQH0_1sigma_width` | 84.1pct − 15.9pct, the quantity the cut is on |
| `fesc_{15.9,median,84.1}pct` | f_esc percentiles from the pipeline |
| `fesc_median_recomputed` | f_esc median recomputed here from Q(H⁰) & L_Hα (auditable) |
| `log_clusterMass_median_Msun` | log of the region's summed cluster mass |
| `passes_QH0_cut` | 1 if `logQH0_1sigma_width < 0.5` (paper's selection) |

## `matched_pairs.csv` columns

`h2_ID, h2_RA, h2_DEC, h2_rGalac_kpc, h2_LHa_erg_s, h2_logLHa, h2_size_pc, h2_class,
sc_ID, sc_RA, sc_DEC, separation_arcsec, sc_age_yr, sc_mass_Msun, sc_EBV, sc_class`

Note: this table is the spatial matching only. Q(H⁰) is a per-*region* summed quantity
(see `fesc_per_region.csv`), not a per-cluster column here.

---
Contact: Jia Wei Teh · jiaweiteh.astro@gmail.com
