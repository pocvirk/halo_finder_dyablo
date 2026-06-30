#!/usr/bin/env python3
"""Plot the FoF halo mass function and compare to theory (Sheth-Tormen +
Tinker08) for a dyablo snapshot.

Reads:
  - the halo catalog written by `fof` (the _masst_NNNNN files under <catalog_dir>)
  - the dyablo restart .ini       -> cosmology (Om, Ob, h0)
  - the dyablo state .h5          -> current scale factor (scalar_data/aexp)
  - the .ini [mesh]/xmax          -> box size (SI metres)
  - the .ini [particle_grid]/total_mass + nx**3 -> DM particle mass

Output: hmf_compare.png next to this script (or wherever -o points).

Usage:
  python3 plot_hmf.py \\
      --catalog outputs/ \\
      --ini     /path/to/restart_star_NNNNNN.ini \\
      --h5      /path/to/star_iterNNNNNNN.h5

  or, reusing the same inputs you already gave to `fof`:
  python3 plot_hmf.py --fofin fof.in
  (--catalog, --ini, --h5, --particle-h5 and --b are taken from fof.in;
   any of them given explicitly on the command line override the fof.in value.)

Optional:
  --sigma8 0.81   (Planck-ish, not in the .ini)
  --ns     0.96
  --hmf-model SMT  (any model name accepted by `hmf`)
  --out    hmf_compare.png

Dependencies: numpy, scipy, h5py, matplotlib, hmf, astropy.
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import h5py
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.io import FortranFile


# ---------------------------------------------------------------------------
# physical constants
MSUN_KG = 1.98892e30
MPC_M   = 3.0856775814913673e22


# ---------------------------------------------------------------------------
# halo catalog reading (Fortran sequential unformatted; matches modio.f90).
#
# _masst_NNNNN per-halo record: int64 id, int32 mass, 3*float32 (x, y, z).
HALO_DTYPE = np.dtype([
    ("id",   "<i8"),
    ("mass", "<i4"),
    ("x",    "<f4"),
    ("y",    "<f4"),
    ("z",    "<f4"),
])


def read_masst_file(path: Path) -> np.ndarray:
    halos = []
    with FortranFile(str(path), "r") as f:
        _ = f.read_ints(np.int32)              # header: ns
        while True:
            try:
                rec = f.read_record(HALO_DTYPE)
            except Exception:
                break
            if rec.size == 0:
                break
            halos.append(rec)
    if not halos:
        return np.empty(0, dtype=HALO_DTYPE)
    return np.concatenate(halos).astype(HALO_DTYPE)


def read_catalog(catalog_dir: Path, root: str = "") -> np.ndarray:
    if root:
        pattern = f"{root}_masst_*"
    else:
        pattern = "*_masst_*"
    files = sorted(catalog_dir.glob(pattern))
    if not files:
        raise FileNotFoundError(
            f"No _masst_ files matching {pattern!r} under {catalog_dir}"
        )
    parts = [read_masst_file(p) for p in files]
    return np.concatenate(parts).astype(HALO_DTYPE)


# ---------------------------------------------------------------------------
# tiny .ini parser (just the values we need)
def parse_ini(path: Path) -> dict[str, dict[str, str]]:
    section = None
    out: dict[str, dict[str, str]] = {}
    with open(path) as f:
        for raw in f:
            line = raw.split(";", 1)[0].strip()
            if not line:
                continue
            m = re.match(r"\[([^]]+)\]", line)
            if m:
                section = m.group(1).strip()
                out.setdefault(section, {})
                continue
            if "=" not in line or section is None:
                continue
            k, v = line.split("=", 1)
            out[section][k.strip()] = v.strip()
    return out


def _f(s: str) -> float:
    return float(s.split()[0])


# ---------------------------------------------------------------------------
# fof.in parser. Field order is fixed and must match fof.f90:58-74.
#   1 root        2 code_index   3 pathinput   4 namepart   5 nameinfo
#   6 grpsize     7 perco        8 Mmin        9 Mmax
#   10-16 star/metal/outcube/dofof/readfromcube/dotimings/invent_ids
#   17 outdir
# String fields (path/name/dir) are read whole-line in Fortran (A512); the
# numeric `perco` is list-directed, so we take its first whitespace token.
def parse_fofin(path: Path) -> dict[str, str]:
    with open(path) as fh:
        lines = [raw.rstrip("\n") for raw in fh]

    def whole(i: int) -> str:
        return lines[i].strip() if i < len(lines) else ""

    def token(i: int) -> str:
        return lines[i].split()[0] if i < len(lines) and lines[i].split() else ""

    return {
        "root":      whole(0),
        "pathinput": whole(2),
        "namepart":  whole(3),
        "nameinfo":  whole(4),
        "perco":     token(6),
        "outdir":    whole(16),
    }


# ---------------------------------------------------------------------------
def warn_if_inputs_inconsistent(sim_files: dict) -> None:
    """Warn (very visibly) if the SIM-SOURCE files do not all live in the same
    directory, or if any is missing. Does NOT abort — it prints and returns.

    `sim_files` maps label -> Path for the files that MUST come from the same
    simulation: the state dump (.h5, gives the redshift), its config (.ini, gives
    the box size) and the particle file (gives the DM particle mass). The halo
    CATALOG is deliberately NOT passed in here: it is routinely kept in a separate
    directory. Feeding files from different sims (e.g. data from an 8 Mpc box with
    the .ini of a 4 Mpc box) silently corrupts the box size and rescales the whole
    mass function by (L_wrong / L_right)^3.
    """
    items   = {lab: p for lab, p in sim_files.items() if p is not None}
    missing = {lab: p for lab, p in items.items() if not Path(p).exists()}
    dirs    = {lab: str(Path(p).resolve().parent) for lab, p in items.items()}
    if len(set(dirs.values())) == 1 and not missing:
        print(f"[consistency OK] sim data + .ini co-located in: {next(iter(dirs.values()))}")
        return
    bar = "#" * 78
    print("\n" + bar)
    print("##  W A R N I N G   —   POSSIBLY INCONSISTENT SIMULATION INPUTS")
    print(bar)
    if missing:
        print("##  Missing file(s):")
        for lab, p in missing.items():
            print(f"##    - {lab}: {p}")
    if len(set(dirs.values())) > 1:
        print("##  These sim-source files do NOT all live in the same directory:")
        for lab, d in dirs.items():
            print(f"##    - {lab:<24s} dir: {d}")
        print("##  => They may come from DIFFERENT simulations. If so, the BOX SIZE read")
        print("##     from the .ini will not match the data, and the mass function will be")
        print("##     rescaled by (L_wrong / L_right)^3. CHECK THIS before trusting the plot.")
    print(f"{bar}\n##  Proceeding anyway — but treat the result as suspect.")
    print(bar + "\n")


def main() -> int:
    ap = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                 description=__doc__)
    ap.add_argument("--fofin", type=Path, default=None,
                    help="a fof.in input file; --catalog, --ini, --h5, "
                         "--particle-h5 and --b are derived from it unless given "
                         "explicitly on the command line")
    ap.add_argument("--catalog", type=Path, default=None,
                    help="directory holding *_masst_NNNNN files "
                         "(or taken from fof.in outdir)")
    ap.add_argument("--ini", type=Path, default=None,
                    help="dyablo restart .ini (or taken from fof.in)")
    ap.add_argument("--h5",  type=Path, default=None,
                    help="dyablo snapshot .h5 (the one with scalar_data/aexp); "
                         "with --fofin, derived from the particle file name")
    ap.add_argument("--omegam", type=float, default=None,
                    help="total matter density Omega_m (read from .ini if present, "
                         "else default 0.3089 / Planck15)")
    ap.add_argument("--h", type=float, default=None,
                    help="dimensionless Hubble constant h = H0/100 (not in the .ini; "
                         "default 0.6774 / Planck15)")
    ap.add_argument("--sigma8", type=float, default=0.81,
                    help="sigma_8 (default 0.81; not in the .ini)")
    ap.add_argument("--ns", type=float, default=0.96,
                    help="primordial spectral index n_s (default 0.96)")
    ap.add_argument("--hmf-model", default="SMT",
                    help="hmf fitting function (default SMT = Sheth-Tormen)")
    ap.add_argument("--out", type=Path, default=Path("hmf_compare.png"),
                    help="output PNG path")
    ap.add_argument("--b", type=float, default=None,
                    help="linking length (percolation parameter); "
                         "auto-detected from <catalog>/*.inp if omitted")
    ap.add_argument("--particle-h5", type=Path, default=None,
                    help="particle .h5 file (needed to auto-detect DM mass when "
                         "[particle_grid] section is absent from the .ini)")
    args = ap.parse_args()

    # ---- fill unset inputs from fof.in (CLI always wins) ----
    if args.fofin is not None:
        fin = parse_fofin(args.fofin)
        base = args.fofin.resolve().parent          # fof runs from fof.in's dir
        pathinput = Path(fin["pathinput"]) if fin["pathinput"] else base

        if args.catalog is None:
            outp = Path(fin["outdir"] or ".")
            args.catalog = outp if outp.is_absolute() else base / outp
        if args.ini is None and fin["nameinfo"]:
            args.ini = pathinput / fin["nameinfo"]
        if args.particle_h5 is None and fin["namepart"]:
            args.particle_h5 = pathinput / fin["namepart"]
        if args.h5 is None and fin["namepart"]:
            # fof.in names the *particle* file; the state file carrying
            # scalar_data/aexp is the same name minus dyablo's
            # '_particles_particles' infix (e.g. star_particles_particles_iterN
            # -> star_iterN).
            state_name = fin["namepart"].replace("_particles_particles", "")
            args.h5 = pathinput / state_name
        if args.b is None and fin["perco"]:
            try:
                args.b = float(fin["perco"])
            except ValueError:
                pass

    missing = [name for name, val in (("--catalog", args.catalog),
                                      ("--ini", args.ini),
                                      ("--h5", args.h5)) if val is None]
    if missing:
        ap.error("missing required input(s): " + ", ".join(missing)
                 + " — provide them explicitly or pass --fofin")
    if not args.h5.exists():
        ap.error(f"snapshot state file not found: {args.h5}\n"
                 "  (derived from fof.in particle name; pass --h5 explicitly if "
                 "your state dump is named differently)")

    # *** consistency guard ***  the .ini (box size), the state .h5 (redshift) and
    # the particle file must all come from the SAME simulation; warn loudly if not.
    # The halo catalog is intentionally excluded (routinely kept in a separate dir).
    warn_if_inputs_inconsistent({"state data (.h5)":   args.h5,
                                 "config (.ini, box)":  args.ini,
                                 "particle file":       args.particle_h5})

    # Linking length: try CLI, else parse the .inp echo file in the catalog dir.
    b_perco = args.b
    if b_perco is None:
        for inp in args.catalog.glob("*.inp"):
            with open(inp) as fh:
                for line in fh:
                    try:
                        v = float(line.strip())
                    except ValueError:
                        continue
                    if 0.0 < v < 1.0:
                        b_perco = v
                        break
            if b_perco is not None:
                break

    # ---- cosmology ----
    ini = parse_ini(args.ini)
    cos = ini.get("cosmology", {})
    mesh = ini["mesh"]
    # omegam/h: not stored in dyablo .ini; use CLI or Planck15 defaults.
    Om0 = args.omegam if args.omegam is not None else _f(cos["omegam"]) if "omegam" in cos else 0.3089
    Ob0 = _f(cos["omegab"]) if "omegab" in cos else 0.049
    h = args.h if args.h is not None else (_f(cos["h0"]) * MPC_M / 1000.0 / 100.0 if "h0" in cos else 0.6774)
    H0 = h * 100.0
    print(f"Cosmology: Om0={Om0:.4f}, Ob0={Ob0:.4f}, H0={H0:.2f} km/s/Mpc (h={h:.4f})")
    print(f"           sigma_8={args.sigma8}, n_s={args.ns} (user-supplied)")

    # ---- redshift ----
    with h5py.File(args.h5, "r") as f:
        aexp = float(f["scalar_data"].attrs["aexp"])
    z = 1.0 / aexp - 1.0
    print(f"Snapshot: aexp={aexp:.6f}  ->  z={z:.4f}")

    # ---- box ----
    Lbox_m = _f(mesh["xmax"]) - _f(mesh.get("xmin", "0"))
    Lbox_Mpc = Lbox_m / MPC_M
    Vbox_Mpc3 = Lbox_Mpc ** 3
    print(f"Box: L={Lbox_Mpc:.4f} Mpc (= {Lbox_Mpc*h:.3f} Mpc/h)")

    # ---- particle mass ----
    pg = ini.get("particle_grid", {})
    if "total_mass" in pg and "nx" in pg:
        m_p_kg = _f(pg["total_mass"]) / int(_f(pg["nx"])) ** 3
    else:
        # [particle_grid] absent (older dyablo restarts): read median mass from HDF5.
        pfile = args.particle_h5 or args.h5
        with h5py.File(pfile, "r") as f:
            if "mass" not in f:
                raise SystemExit(f"No 'mass' dataset in {pfile}; pass --particle-h5 "
                                 "pointing to the particle HDF5 file.")
            masses = f["mass"][:]
        m_p_kg = float(np.median(masses))
        print(f"[particle_grid] absent; DM particle mass from HDF5 median: {m_p_kg:.4e} kg")
    m_p_Msun = m_p_kg / MSUN_KG
    print(f"DM particle mass: {m_p_kg:.4e} kg  =  {m_p_Msun:.4e} Msun")

    # ---- catalog ----
    cat = read_catalog(args.catalog)
    M_Msun = cat["mass"].astype(np.float64) * m_p_Msun
    print(f"Catalog: {cat.size} halos, M in [{M_Msun.min():.3e}, {M_Msun.max():.3e}] Msun")

    # ---- measured dn/dlogM ----
    nbins = 16
    lgM_edges = np.linspace(np.log10(M_Msun.min()),
                            np.log10(M_Msun.max()) * 1.001,
                            nbins + 1)
    counts, _ = np.histogram(np.log10(M_Msun), bins=lgM_edges)
    dlogM = np.diff(lgM_edges)
    lgM_centres = 0.5 * (lgM_edges[1:] + lgM_edges[:-1])
    dn_dlogM = counts / (Vbox_Mpc3 * dlogM)
    err = np.sqrt(counts) / (Vbox_Mpc3 * dlogM)

    # ---- theoretical HMF (hmf + astropy + Eisenstein-Hu transfer) ----
    from hmf import MassFunction
    from astropy.cosmology import FlatLambdaCDM
    cosmo = FlatLambdaCDM(H0=H0, Om0=Om0, Ob0=Ob0, Tcmb0=2.725)
    common = dict(z=z, cosmo_model=cosmo, sigma_8=args.sigma8, n=args.ns,
                  Mmin=np.log10(M_Msun.min() * h) - 0.5,
                  Mmax=np.log10(M_Msun.max() * h) + 0.5,
                  dlog10m=0.02, transfer_model="EH")
    mf = MassFunction(hmf_model=args.hmf_model, **common)
    mf_t = MassFunction(hmf_model="Tinker08", **common)

    M_theory_Msun = mf.m / h
    dn_th_st = mf.dndlog10m * h**3              # (Mpc/h)^-3 -> Mpc^-3
    dn_th_t  = mf_t.dndlog10m * h**3

    # ---- plot ----
    fig, ax = plt.subplots(figsize=(7.5, 5.5))

    ax.plot(np.log10(M_theory_Msun), dn_th_st, "-", color="C0", lw=1.5,
            label=f"{args.hmf_model}")
    ax.plot(np.log10(M_theory_Msun), dn_th_t,  "--", color="C3", lw=1.5,
            label="Tinker08")
    keep = counts > 0
    ax.errorbar(lgM_centres[keep], dn_dlogM[keep], yerr=err[keep],
                fmt="o", ms=5, capsize=3, color="black",
                label=f"FoF (N={cat.size})")

    ax.set_yscale("log")
    ax.set_xlabel(r"$\log_{10}(M / M_\odot)$")
    ax.set_ylabel(r"$dn/d\log M \;\; [\mathrm{Mpc}^{-3}]$")
    b_str = f"b={b_perco:.3g}" if b_perco is not None else "b=?"
    ax.set_title(f"Halo mass function — dyablo FoF ({b_str}), z={z:.3f}\n"
                 f"L={Lbox_Mpc:.2f} Mpc, $m_p$={m_p_Msun:.2e} $M_\\odot$")
    ax.legend()
    ax.grid(True, which="both", alpha=0.3)
    fig.tight_layout()
    fig.savefig(args.out, dpi=130)
    print(f"\nSaved {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
