#!/usr/bin/env python3
"""
Wildfire RMT analysis: level repulsion in fire-interval spacing.
Reservoir criterion test on single-forest fire-scar chronologies (FHX2 format).

Reads .fhx files from ../data/, extracts fire events by the Grissino-Mayer
10% filter, restricts to the pre-1900 natural fire regime, and computes the
spacing ratio <r>, the coefficient of variation, and significance against a
size-matched Poisson null ensemble. Writes a summary CSV to ../results/.

Reference values: Poisson 0.386, GOE 0.531, GUE 0.603.
"""
import re
import glob
import csv
import numpy as np
from scipy.interpolate import UnivariateSpline
import warnings
warnings.filterwarnings("ignore")

# FHX2 uppercase fire-scar season codes (A=early earlywood, D=dormant,
# E=earlywood, L=latewood, M=middle earlywood, U=undetermined season).
SCAR = set("ADELMU")


def parse_fhx(path):
    """Parse an FHX2 file. Return (rows, site_name).
    rows = list of (year, n_scars, n_recording) per data line."""
    raw = open(path, encoding="latin-1").read().replace("\r", "")
    lines = raw.split("\n")
    idx = [k for k, l in enumerate(lines) if l.strip().upper().startswith("FHX2")]
    if not idx:
        return None, ""
    i = idx[0]
    hdr = lines[i + 1].split()
    try:
        startyr, nsamp, namelen = int(hdr[0]), int(hdr[1]), int(hdr[2])
    except (ValueError, IndexError):
        return None, ""
    data = lines[i + 2 + namelen:]
    name = ""
    for l in lines[:2]:
        if ":" in l:
            name = l.split(":", 1)[1].strip()
            break
    rows = []
    for ln in data:
        m = re.match(r"^([.\|{}\[\]A-Za-z]+)\s+(\d{3,4})\s*$", ln)
        if not m:
            continue
        codes, yr = m.group(1), int(m.group(2))
        n_scar = sum(1 for c in codes if c in SCAR)
        n_rec = sum(1 for c in codes if c not in ".{}[]")
        rows.append((yr, n_scar, n_rec))
    return rows, name


def events(rows, frac=0.10):
    """Grissino-Mayer filter: fire years with scars on >= frac of recorders."""
    yrs = [yr for yr, nsc, nrec in rows
           if nsc > 0 and nrec > 0 and nsc >= frac * nrec]
    return np.array(sorted(set(yrs)))


def spacing_ratio(intervals):
    s = np.asarray(intervals, float)
    s = s / np.mean(s)
    r = np.minimum(s[:-1], s[1:]) / np.maximum(s[:-1], s[1:])
    return float(np.mean(r))


def cv(intervals):
    s = np.asarray(intervals, float)
    return float(np.std(s) / np.mean(s))


def unfold_r(ev):
    """Spacing ratio after removing the slow trend in fire frequency."""
    if len(ev) < 6:
        return float("nan")
    idx = np.arange(1, len(ev) + 1)
    spl = UnivariateSpline(ev, idx, s=len(ev) / 3.0)
    d = np.diff(spl(ev))
    d = d[d > 0]
    return spacing_ratio(d)


def poisson_significance(r_obs, n, seed=3, nsim=4000):
    """Sigma of r_obs above the size-matched Poisson null."""
    if n < 5:
        return float("nan")
    rng = np.random.default_rng(seed)
    rs = np.array([spacing_ratio(rng.exponential(1, n - 1)) for _ in range(nsim)])
    return float((r_obs - rs.mean()) / rs.std())


def analyze_site(path, frac=0.10, cutoff=1900):
    rows, name = parse_fhx(path)
    if rows is None:
        return None
    ev = events(rows, frac)
    ev = ev[ev < cutoff]            # pre-1900 natural regime
    if len(ev) < 8:
        return None
    iv = np.diff(ev).astype(float)
    return dict(
        name=name, n=len(ev),
        r=round(spacing_ratio(iv), 3),
        cv=round(cv(iv), 2),
        fri=round(float(np.mean(iv)), 1),
        unf_r=round(unfold_r(ev), 3),
        sig=round(poisson_significance(spacing_ratio(iv), len(ev)), 1),
        span=(int(ev.min()), int(ev.max())),
    )


def main():
    paths = sorted(glob.glob("../data/*.fhx"))
    results = {}
    for p in paths:
        code = re.search(r"(us[\w-]+?)\.fhx", p).group(1)
        res = analyze_site(p)
        if res:
            results[code] = res
            print("%-10s %-15s n=%3d <r>=%.3f CV=%.2f sig=%+.1f"
                  % (code, res["name"][:15], res["n"], res["r"], res["cv"], res["sig"]))

    adequate = {k: v for k, v in results.items() if v["n"] >= 30 and v["sig"] >= 2}
    rs = [v["r"] for v in adequate.values()]
    print("\nAdequate sites (n>=30, >+2 sigma): %d" % len(adequate))
    print("<r> mean=%.3f median=%.3f std=%.3f"
          % (np.mean(rs), np.median(rs), np.std(rs)))
    print("Reference: Poisson 0.386, GOE 0.531, GUE 0.603")

    # write CSV
    with open("../results/fire_results.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["code", "site", "n_events", "spacing_ratio", "cv",
                    "mean_fri_yr", "unfolded_r", "poisson_sigma",
                    "span_start", "span_end", "adequate"])
        for k, v in sorted(results.items()):
            w.writerow([k, v["name"], v["n"], v["r"], v["cv"], v["fri"],
                        v["unf_r"], v["sig"], v["span"][0], v["span"][1],
                        "yes" if (v["n"] >= 30 and v["sig"] >= 2) else "no"])
    print("\nWrote ../results/fire_results.csv")


if __name__ == "__main__":
    main()
