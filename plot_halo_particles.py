#!/usr/bin/env python3
"""Read FoF _strct_NNNNN files and plot a subsample of halo particles.

For each halo, every Nth particle (default N=10) is drawn. Each halo gets a
single fixed colour; the colour mapping is random-shuffled so spatially-
adjacent halos in the catalog don't end up with near-identical hues.

Three projections (xy, xz, yz) are plotted side by side. Positions are in
normalized box units ([0, 1]^3).

Usage:
  python3 plot_halo_particles.py --catalog outputs/
  python3 plot_halo_particles.py --catalog outputs/ --stride 10 \\
      --min-mass 500 --out halo_particles.png

The _strct_ file format (per MPI rank, Fortran sequential unformatted, see
modio.f90):
  record 1                            int32   ns       (number of halos kept)
  for each halo (ns times):
    record                            int32   mass
    record                            float32 (mass, 3)  positions
    record                            float32 (mass, 3)  velocities
    record                            int32   (mass,)    original /id values
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.io import FortranFile


def read_strct_file(path: Path):
    """Return a list of (mass, positions [N,3], ids [N]) for one rank file."""
    halos = []
    with FortranFile(str(path), "r") as f:
        ns = int(f.read_ints(np.int32)[0])
        for _ in range(ns):
            mass = int(f.read_ints(np.int32)[0])
            pos = f.read_reals(np.float32).reshape((mass, 3))
            _vel = f.read_reals(np.float32)   # not used here
            ids = f.read_ints(np.int32)
            halos.append((mass, pos, ids))
    return halos


def read_strct_dir(d: Path, root: str = ""):
    pattern = f"{root}_strct_*" if root else "*_strct_*"
    files = sorted(d.glob(pattern))
    if not files:
        raise FileNotFoundError(f"No _strct_ files matching {pattern!r} in {d}")
    halos = []
    for fp in files:
        halos.extend(read_strct_file(fp))
    return halos


def main() -> int:
    ap = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                 description=__doc__)
    ap.add_argument("--catalog", type=Path, required=True,
                    help="directory holding *_strct_NNNNN files")
    ap.add_argument("--stride", type=int, default=10,
                    help="plot every stride-th particle of each halo (default 10)")
    ap.add_argument("--min-mass", type=int, default=0,
                    help="only plot halos with at least this many particles (default 0)")
    ap.add_argument("--out", type=Path, default=Path("halo_particles.png"),
                    help="output PNG path")
    ap.add_argument("--seed", type=int, default=0,
                    help="RNG seed for the colour shuffle (default 0)")
    args = ap.parse_args()

    halos = read_strct_dir(args.catalog)
    print(f"Read {len(halos)} halos from {args.catalog}")

    # Filter by mass if requested
    if args.min_mass > 0:
        halos = [h for h in halos if h[0] >= args.min_mass]
        print(f"After --min-mass {args.min_mass}: {len(halos)} halos")

    # Stable, shuffled colours: one per halo.
    rng = np.random.default_rng(args.seed)
    cmap = plt.get_cmap("gist_rainbow")
    colors = cmap(rng.uniform(0.0, 1.0, size=len(halos)))

    fig, axes = plt.subplots(1, 3, figsize=(15, 5.4), sharex=True, sharey=True)
    proj = [(0, 1, "x", "y"),
            (0, 2, "x", "z"),
            (1, 2, "y", "z")]

    total_plotted = 0
    for ih, (mass, pos, _ids) in enumerate(halos):
        sub = pos[::args.stride]
        if sub.size == 0:
            continue
        total_plotted += sub.shape[0]
        col = colors[ih:ih + 1]    # array of shape (1,4) — broadcasts
        for ax, (i, j, _li, _lj) in zip(axes, proj):
            ax.scatter(sub[:, i], sub[:, j], s=1.0, c=col, alpha=0.85,
                       edgecolors="none")

    for ax, (_i, _j, li, lj) in zip(axes, proj):
        ax.set_xlim(0.0, 1.0)
        ax.set_ylim(0.0, 1.0)
        ax.set_aspect("equal")
        ax.set_xlabel(li)
        ax.set_ylabel(lj)
        ax.set_facecolor("black")

    fig.suptitle(f"FoF halo particles (1/{args.stride} subsample) — "
                 f"{len(halos)} halos, {total_plotted} points plotted",
                 fontsize=11)
    fig.tight_layout()
    fig.savefig(args.out, dpi=140, facecolor="white")
    print(f"Saved {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
