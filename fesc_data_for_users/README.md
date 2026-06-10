# NGC 628 escape-fraction data

Matched H II-region ↔ star-cluster catalogue and escape fractions from the LEGUS × SIGNALS
pilot study of NGC 628 (Teh et al. 2023, MNRAS 524, 1191). Only `numpy` is needed — no
pipeline, no SLUG2 libraries, no multi-GB data.

**The CSVs in this folder are already complete — use them directly.** Re-run the script only
to regenerate them.

- `matched_pairs.csv` — one row per (H II region, matched cluster): the ID linkage plus each
  object's position and properties.
- `fesc_per_region.csv` — one row per region in the final sample: L_Hα, Q(H⁰), and f_esc,
  each with percentiles.

## Escape fraction

```
Q_Hα  = L_Hα × 7.31e11           # ionising-photon rate from observed Hα (Kennicutt)
f_esc = (Q(H⁰) − Q_Hα) / Q(H⁰)   # Q(H⁰) = summed budget of the matched clusters
```

Q(H⁰) is summed over all clusters matched to a region (10⁵ Monte-Carlo draws from each
cluster's posterior). `fesc_median_recomputed` applies this formula to the stored samples and
reproduces `fesc_median` exactly, so you can swap in your own assumptions.

## ⚠ Use the quality cut

`fesc_per_region.csv` is the **pre-cleaning** sample; poorly-constrained regions carry wild
negative f_esc. The paper keeps only regions with 1σ log-Q(H⁰) width < 0.5 dex.
**Filter on `passes_QH0_cut == 1` before averaging or histogramming f_esc.**

## Regenerate the CSVs

Put `combined_catalogue.npy` (~1.4 MB) and `GenKroupc1234LHa_QH0_noCons.npy` (~222 MB) in one
folder, then:

```
python extract_fesc_data.py --dat /path/to/npy --out .
```

Getting the inputs: request both from the author (below). `combined_catalogue.npy` you can
rebuild for free from the repo root (needs only the bundled catalogues):

```
python -c "from src.tools.create_combined_table import create_combined_table as f; f()"
```

`GenKroupc1234LHa_QH0_noCons.npy` is built from the SLUG2 libraries (not public) — request it.

## Columns — `fesc_per_region.csv`

| Column | Meaning |
|--------|---------|
| `h2_ID` | SIGNALS H II-region ID |
| `n_clusters`, `matched_cluster_IDs` | number and LEGUS IDs of matched clusters |
| `logLHa`, `LHa_erg_s` | extinction-corrected Hα luminosity |
| `logQH0_{2.3,15.9,50,84.1,97.7}pct` | log Q(H⁰) percentiles (2σ, 1σ, median) |
| `logQH0_1sigma_width` | 84.1pct − 15.9pct; the quantity the cut uses |
| `fesc_{15.9,median,84.1}pct` | f_esc percentiles from the pipeline |
| `fesc_median_recomputed` | f_esc median recomputed here from Q(H⁰) & L_Hα |
| `log_clusterMass_median_Msun` | log of the region's summed cluster mass |
| `passes_QH0_cut` | 1 if `logQH0_1sigma_width < 0.5` (paper's cut) |

## Columns — `matched_pairs.csv`

`h2_ID, h2_RA, h2_DEC` (region ID + position), `h2_rGalac_kpc` (galactocentric radius),
`h2_LHa_erg_s, h2_logLHa` (Hα luminosity), `h2_size_pc, h2_class`; then per matched cluster
`sc_ID, sc_RA, sc_DEC`, `separation_arcsec`, `sc_age_yr, sc_mass_Msun, sc_EBV, sc_class`.
Q(H⁰) is a per-*region* sum, not a per-cluster column here — see `fesc_per_region.csv`.

---
Contact: Jia Wei Teh · jiaweiteh.astro@gmail.com
