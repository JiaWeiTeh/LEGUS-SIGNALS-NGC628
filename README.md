# LEGUS–SIGNALS–NGC 628

Analysis code and figures for **Teh et al. (2023)**, "Constraining the LyC escape fraction
from LEGUS star clusters with SIGNALS H II region observations: a pilot study of NGC 628,"
*MNRAS* **524**, 1191.

[![ADS](https://img.shields.io/badge/ADS-2023MNRAS.524.1191T-1f6feb.svg)](https://ui.adsabs.harvard.edu/abs/2023MNRAS.524.1191T/abstract)
[![arXiv](https://img.shields.io/badge/arXiv-2306.05457-b31b1b.svg)](https://arxiv.org/abs/2306.05457)
[![DOI](https://img.shields.io/badge/DOI-10.1093%2Fmnras%2Fstad1780-blue.svg)](https://doi.org/10.1093/mnras/stad1780)

## Install

Python 3.10+ and a LaTeX install (figures use `usetex`):

```
python -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
```

## Reproduce the figures

Run from the repo root:

```
python make_figures.py            # all figures → src/figs/
python make_figures.py fg5 fg7    # only matching scripts
python make_figures.py --list     # list, don't run
```

`make_figures.py` runs each `src/plots/fg*.py` in its own subprocess, skips past failures,
and prints a summary. You can also run one directly (`python src/plots/fg5_LHa_vs_QH0.py`).
Scripts add the repo root to `sys.path` (Spyder works too) and read cached arrays from
`src/dat/`, so the data must be in place — see [Data](#data).

## Just want the escape-fraction numbers?

`fesc_data_for_users/` ships the results as CSVs (`matched_pairs.csv`, `fesc_per_region.csv`)
— no pipeline or large data needed. To regenerate them from the `.npy` arrays:

```
cd fesc_data_for_users
python extract_fesc_data.py --dat ../src/dat --out .
```

`fesc_per_region.csv` gives, per region, L_Hα, Q(H⁰), and
`f_esc = (Q(H⁰) − L_Hα·7.31e11) / Q(H⁰)` with percentiles; filter on `passes_QH0_cut == 1`
before averaging. See [`fesc_data_for_users/README.md`](fesc_data_for_users/README.md).

## Layout

```
make_figures.py        run all figure scripts
lib/                   input data
  LEGUS/  SIGNALS/        bundled catalogues
  SLUG2/                  model libraries (not in Git — see Data)
src/
  paths.py               all paths, relative to the repo root
  tools/                 catalogue builders + analysis helpers
  reg_files/             export DS9 .reg files for visualisation
  plots/                 fg*.py, one per figure (+ alternative_params/, dark/, extras/)
  dat/                   post-processed .npy the plots read (not in Git — see Data)
  figs/                  output figure PDFs
fesc_data_for_users/   standalone f_esc CSV extractor
```

## Data

The plots read post-processed `.npy` arrays from `src/dat/` — mostly the Monte-Carlo Q(H⁰)
distributions (~850 MB) and clusterslug tables (~180 MB). These, plus the raw SLUG2 libraries
(`lib/SLUG2/`), are too large for Git and aren't included; the LEGUS and SIGNALS catalogues
are bundled.

Only `combined_catalogue.npy` rebuilds for free from the bundled catalogues (no SLUG2 needed):

```
python -c "from src.tools.create_combined_table import create_combined_table as f; f()"
```

The clusterslug tables and the heavy Q(H⁰) arrays need the raw SLUG2 libraries — and the
Q(H⁰) arrays are produced by the plot scripts themselves (`reRun = True` in
`fg5_LHa_vs_QH0.py`), not by `src/tools/`.

**To run the full pipeline, request the `src/dat/` bundle** (and SLUG2 libraries) from Jia Wei
Teh at **jiaweiteh.astro@gmail.com**, unpack into `src/dat/`, and run. Full inventory in
[`DATA.md`](DATA.md).
