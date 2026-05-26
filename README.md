# fof_dyablo

A parallel Friends-of-Friends (FoF) halo finder for [dyablo](https://github.com/Dyablo-HPC/Dyablo)
HDF5 snapshots, adapted from the RAMSES pFoF code. It reads dyablo particle
dumps, links particles into halos in parallel with MPI, and writes a halo
catalog plus per-halo observables.

## Build

The finder is MPI + HDF5 Fortran. You need an MPI compiler and the HDF5 Fortran
library (`h5pfc` or `mpif90` against an HDF5 install).

```bash
make                       # uses h5pfc by default (see makefile)
# or override the toolchain:
make MPIFC=mpif90 FC=mpif90 HDF5_DIR=/path/to/hdf5
```

`-DLONGINT` (on by default) selects 8-byte particle IDs, needed when
`N_particles > 2^31`.

## Run

Edit `fof.in` (one parameter per line — see the user guide for the field
layout), then launch with MPI:

```bash
mpirun -np 64 ./fof
# if you have fewer physical cores than ranks, allow hardware threads:
mpirun --use-hwthread-cpus -np 64 ./fof
```

Outputs are written to the `outdir` named in `fof.in`.

## Analysis scripts

- **`plot_hmf.py`** — reads the FoF catalog and overplots analytic halo mass
  functions (Sheth–Tormen, Tinker08) at the snapshot redshift. The quickest
  invocation just reuses your input file:
  ```bash
  python3 plot_hmf.py --fofin fof.in
  ```
  All paths (catalog, .ini, snapshot .h5, particle .h5, linking length) are
  derived from `fof.in`; any can be overridden on the command line. See the
  user guide for details.
- **`plot_halo_particles.py`** — visualizes the particle members of individual
  halos.

Python dependencies: `numpy scipy h5py matplotlib hmf astropy`.

## Documentation

The full user guide lives in [`doc/`](doc/) (`user_guide.tex`, with a prebuilt
`user_guide.pdf`). Rebuild the PDF with `pdflatex user_guide.tex` (requires the
`listings` and `xcolor` LaTeX packages).
