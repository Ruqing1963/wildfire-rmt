# wildfire-rmt

**Level Repulsion in Wildfire Regimes: Fuel as the Reservoir of Isolated Forests**

A Random Matrix Theory Test on Fire-Scar Chronologies

Ruqing Chen, GUT Geoservice Inc., Montreal, Canada. June 2026.

## Summary

This repository tests whether wildfire regimes exhibit Wigner-Dyson level repulsion in the timing of fires, as predicted by the reservoir criterion of an ongoing random matrix theory (RMT) research program. The criterion holds that GOE-class level repulsion appears when, and only when, a single autonomous reservoir charges slowly and discharges substantially in discrete, datable events. A single isolated forest is a natural test on the biosphere: forest fuel accumulates slowly over years to decades (recharge), a fire consumes it (release), and the stand must rebuild fuel before the next fire. Wildfire is the program's first biosphere entry.

## Result

We analyze thirteen single-site fire-scar chronologies from the International Multiproxy Paleofire Database (IMPD), extract fire events by an objective Grissino-Mayer filter (years in which a fire scar is recorded on at least 10% of recording samples), and restrict to the pre-1900 natural fire regime to avoid twentieth-century fire suppression.

Of the thirteen sites, nine have adequate samples (n >= 30 events). All nine show GOE-class level repulsion, with spacing ratios from 0.515 to 0.622 (mean 0.562 +/- 0.050), each +2.6 to +6.3 sigma above the Poisson value (0.386), clustered between the GOE (0.531) and GUE (0.603) values. The four non-significant sites are precisely the four smallest samples (n = 15 to 24), as the criterion predicts. The result is robust to threshold, season-code definition, and unfolding. A control: including post-1900 fires inflates the coefficient of variation toward clustering, the fingerprint of a disrupted reservoir.

## Repository contents

- paper/ : manuscript (paper.tex, paper.pdf) and figure
- code/ : analysis script (fire_rmt.py)
- data/ : thirteen FHX2 fire-scar chronologies
- figures/ : figure in PDF and PNG
- results/ : locked per-site statistics (fire_results.csv)

## Reproducing

    cd code
    python fire_rmt.py

Requires numpy and scipy.

## Reference values

Poisson 0.386, GOE 0.531, GUE 0.603.

## Data source

International Multiproxy Paleofire Database (NOAA NCEI):
https://www.ncei.noaa.gov/products/paleoclimatology/fire-history

## Part of the RMT reservoir-criterion program

Companion: ENSO (https://zenodo.org/records/20797691), mantle plumes, pulsar glitches.
