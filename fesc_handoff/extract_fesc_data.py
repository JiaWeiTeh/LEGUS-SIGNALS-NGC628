#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
extract_fesc_data.py

Self-contained extractor for the NGC 628 LEGUS x SIGNALS escape-fraction data.
Reads the two pipeline arrays and writes plain CSVs. Needs only numpy -- it does
NOT import the `src` package, so it runs anywhere the two .npy files are present.

Inputs (place both in the same folder as this script, or pass --dat):
  combined_catalogue.npy          full HII-region <-> star-cluster spatial matching
  GenKroupc1234LHa_QH0_noCons.npy  per-region Monte-Carlo Q(H0) / f_esc PDFs (final sample)

Outputs:
  matched_pairs.csv     one row per (HII region, matched cluster): IDs + key values
  fesc_per_region.csv   one row per HII region in the final sample: L_Ha, Q(H0), f_esc

Escape fraction (per region; Teh et al. 2023):
  Q_Ha   = L_Ha * 7.31e11           # ionising photon rate implied by observed Halpha
  f_esc  = (Q(H0) - Q_Ha) / Q(H0)   # Q(H0) = summed intrinsic budget of matched clusters
fesc_per_region.csv reports f_esc both as the pipeline's stored percentiles and as a
value recomputed from the columns in the CSV, so the formula is transparent and the
colleague can change the L_Ha->Q_Ha conversion or the f_esc definition herself.
"""

import argparse
import os
import numpy as np

# Column order of each catalogue, copied from src/tools/read_catalogue.py so this
# script stays standalone. Index = position in a catalogue row.
H2_COLS = ['ID', 'RA', 'DEC', 'rGalac', 'LHa', 'HaDIG', 'class',
           'I0', 'Amp', 'sig', 'alpha', 'R2', 'rad', 'ext', 'ext_err',
           'N2Ha', 'N2Ha_err', 'N2Ha_SNR', 'S2Ha', 'S2Ha_err', 'S2Ha_SNR',
           'S2N2', 'S2N2_err', 'S2N2_SNR', 'O3Hb', 'O3Hb_err', 'O3Hb_SNR',
           'O2Hb', 'O2Hb_err', 'O2Hb_SNR', 'O23Hb', 'O23Hb_err', 'O23Hb_SNR',
           'O3O2', 'O3O2_err', 'O3O2_SNR', 'O3N2', 'O3N2_err', 'O3N2_SNR',
           'O2N2', 'O2N2_err', 'O2N2_SNR', 'S2S2', 'S2S2_err', 'S2S2_SNR',
           'dist', 'assoc']

SC_COLS = ['ID', 'x', 'y', 'RA', 'DEC', 'UV', 'dUV', 'U', 'dU', 'B', 'dB',
           'V', 'dV', 'I', 'dI', 'CI', 'age', 'maxage', 'minage',
           'mass', 'minmass', 'maxmass', 'ext', 'maxext', 'minext',
           'chi2UV', 'chi2U', 'chi2B', 'chi2V', 'chi2I', 'reduchi', 'Q',
           'nfil', 'class', 'mclass', 'dist', 'assoc']

H2 = {n: i for i, n in enumerate(H2_COLS)}
SC = {n: i for i, n in enumerate(SC_COLS)}

QHA_PER_LHA = 7.31e11   # Kennicutt: Q(Halpha) = L_Ha * this  (photon s^-1 per erg s^-1)


def build_matched_pairs(sc_cat, h2_cat):
    """One row per (HII region, matched cluster). Cluster Q(H0) is *not* here --
    Q(H0) is a per-region summed Monte-Carlo quantity (see fesc_per_region.csv)."""
    sc_by_id = {float(r[SC['ID']]): r for r in sc_cat}
    header = ['h2_ID', 'h2_RA', 'h2_DEC', 'h2_rGalac_kpc', 'h2_LHa_erg_s',
              'h2_logLHa', 'h2_size_pc', 'h2_class',
              'sc_ID', 'sc_RA', 'sc_DEC', 'separation_arcsec',
              'sc_age_yr', 'sc_mass_Msun', 'sc_EBV', 'sc_class']
    rows = []
    for h in h2_cat:
        assoc = np.atleast_1d(h[H2['assoc']])
        dists = np.atleast_1d(h[H2['dist']])
        if assoc.size == 0:
            continue
        LHa = float(h[H2['LHa']])
        logLHa = np.log10(LHa) if LHa > 0 else np.nan
        for scID, sep in zip(assoc, dists):
            s = sc_by_id.get(float(scID))
            if s is None:
                continue
            rows.append([
                int(h[H2['ID']]), float(h[H2['RA']]), float(h[H2['DEC']]),
                float(h[H2['rGalac']]), LHa, logLHa,
                float(h[H2['rad']]), int(h[H2['class']]),
                int(float(scID)), float(s[SC['RA']]), float(s[SC['DEC']]), float(sep),
                float(s[SC['age']]), float(s[SC['mass']]), float(s[SC['ext']]),
                int(s[SC['class']]),
            ])
    return header, rows


def build_fesc_per_region(clusterdata, h2_cat):
    """One row per HII region in the final f_esc sample (GenKroup...noCons.npy)."""
    # map h2 ID -> matched cluster IDs, for reference in the output
    assoc_by_id = {float(h[H2['ID']]): np.atleast_1d(h[H2['assoc']]) for h in h2_cat}

    header = ['h2_ID', 'n_clusters', 'matched_cluster_IDs',
              'logLHa', 'LHa_erg_s',
              'logQH0_2.3pct', 'logQH0_15.9pct', 'logQH0_median',
              'logQH0_84.1pct', 'logQH0_97.7pct', 'logQH0_1sigma_width',
              'fesc_15.9pct', 'fesc_median', 'fesc_84.1pct',
              'fesc_median_recomputed', 'log_clusterMass_median_Msun',
              'passes_QH0_cut']
    rows = []
    for entry in clusterdata:
        f_esc_pdf, f_esc_pct, qh0_pdf, qh0_pct, LHa_log, h2ID, med_mass = entry
        LHa = 10.0 ** float(LHa_log)
        ids = assoc_by_id.get(float(h2ID), np.array([]))
        # Recompute f_esc straight from the stored Q(H0) samples and L_Ha, so the
        # formula is auditable and the colleague can swap in her own assumptions.
        qHa = LHa * QHA_PER_LHA
        fesc_recomputed_median = float(np.median((np.asarray(qh0_pdf) - qHa) / np.asarray(qh0_pdf)))
        # Paper's selection (fg5/fg7): keep only well-constrained Q(H0), i.e. the
        # 1-sigma log-Q(H0) width (84.1pct - 15.9pct) < 0.5 dex. Regions failing this
        # carry the wild negative f_esc values and are dropped before quoting results.
        qh0_1sig_width = float(qh0_pct[3] - qh0_pct[1])
        passes = qh0_1sig_width < 0.5
        rows.append([
            int(h2ID), int(np.size(ids)),
            ';'.join(str(int(float(x))) for x in np.atleast_1d(ids)),
            float(LHa_log), LHa,
            *[float(v) for v in qh0_pct],          # 5 log-Q(H0) percentiles
            qh0_1sig_width,
            *[float(v) for v in f_esc_pct],         # f_esc 15.9 / 50 / 84.1
            fesc_recomputed_median,
            float(med_mass),
            int(passes),
        ])
    rows.sort(key=lambda r: r[0])
    return header, rows


def write_csv(path, header, rows):
    with open(path, 'w') as f:
        f.write(','.join(header) + '\n')
        for r in rows:
            f.write(','.join('' if v is None else
                             (f'{v:.6g}' if isinstance(v, float) else str(v))
                             for v in r) + '\n')


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--dat', default='.', help='folder holding the two .npy files (default: .)')
    p.add_argument('--out', default='.', help='output folder for the CSVs (default: .)')
    args = p.parse_args()

    combined = os.path.join(args.dat, 'combined_catalogue.npy')
    noCons = os.path.join(args.dat, 'GenKroupc1234LHa_QH0_noCons.npy')

    sc_cat, h2_cat = np.load(combined, allow_pickle=True)

    hdr, rows = build_matched_pairs(sc_cat, h2_cat)
    out1 = os.path.join(args.out, 'matched_pairs.csv')
    write_csv(out1, hdr, rows)
    print(f'matched_pairs.csv      : {len(rows):5d} (region, cluster) pairs  ->  {out1}')

    if os.path.exists(noCons):
        clusterdata = np.load(noCons, allow_pickle=True)
        hdr, rows = build_fesc_per_region(clusterdata, h2_cat)
        out2 = os.path.join(args.out, 'fesc_per_region.csv')
        write_csv(out2, hdr, rows)
        print(f'fesc_per_region.csv    : {len(rows):5d} HII regions (final sample) ->  {out2}')
    else:
        print(f'NOTE: {noCons} not found -- skipped fesc_per_region.csv.')
        print('      Send it too if she wants the computed Q(H0)/f_esc values.')


if __name__ == '__main__':
    main()
