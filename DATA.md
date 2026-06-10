# Data files

This repo keeps only the small input files in Git. The big files — the SLUG2 model
libraries and the processed data the figures read — are too large for GitHub, so they
live only on your own machine and are git-ignored. This page says what each file is and
where it goes.

## In Git (small inputs)

- `lib/LEGUS/Catalog.csv` — the LEGUS star-cluster catalogue for NGC 628
  (columns explained in `lib/LEGUS/LEGUS_description.txt`).
- `lib/SIGNALS/NGC628_catalog_WCS_corr.fits` — the SIGNALS Hα H II-region catalogue
  (columns explained in `lib/SIGNALS/SIGNALS_description.txt`).

## Not in Git (kept on your machine)

- `lib/SLUG2/` — the SLUG2 stellar-population model libraries the analysis is built on:
  - `cluster_slug/` — lookup tables that turn a cluster's photometry into estimates of its
    ionising-photon output, mass, age, and reddening (one set for Geneva tracks, one for
    Padova).
  - `cluster_lib/` — the underlying SLUG2 model clusters (used by figures 4 and 9).

  Get these from the SLUG2 project or your own SLUG2 runs and place them under `lib/SLUG2/`.
- `src/dat/*.npy` — the processed data the figures read: the matched catalogue, the
  Monte-Carlo ionising-photon results, the lookup tables in array form, and a few summary
  numbers. The code writes these; they aren't versioned.
- `src/figs/` and other `*.pdf` — the output figures.

## Rebuilding the processed data

`combined_catalogue.npy` (the matched catalogue) needs only the two bundled catalogues —
no SLUG2 data. From the repo root:

```
python -c "from src.tools.create_combined_table import create_combined_table as f; f()"
```

Everything heavier — the lookup tables and the large `*LHa_QH0_noCons*.npy` Monte-Carlo
results — needs the `lib/SLUG2/` libraries. The Monte-Carlo results are produced by the
figure scripts themselves (set `reRun = True` in, e.g., `fg5_LHa_vs_QH0.py`), not by the
tools in `src/tools/`. Everything is written into `src/dat/`.

## Getting the data

Don't want to rebuild it? The `src/dat/` files and the SLUG2 libraries can be sent on
request — contact Jia Wei Teh at jiaweiteh.astro@gmail.com.

For the escape-fraction CSVs and what every column means, see
[`fesc_data_for_users/README.md`](fesc_data_for_users/README.md).
