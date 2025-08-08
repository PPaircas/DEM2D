# DEM2D (SUDEMOD2) Documentation

### Overview
Discrete-Element Modeling (2D) for acoustic wave propagation with contact mechanics. The code assembles mass/stiffness for linear triangular elements, constructs contact interfaces with springs, applies a Ricker-wavelet source, and outputs velocity snapshots.

- Main executable: `bin/DEMEXE`
- Language: C (OpenMP)
- Libraries: Seismic Unix (SU) `su`, `par`, `cwp`; internal lib `lib/libmyfun.a`
- Key headers: `include/DemRun.h`, `include/Element.h`, `include/basefun.h`

### Quick Start
- **Build library**: `make -C lib`
- **Build executable**: set `SU_HOME` (and optionally `SU_HOME1`) in `SRC/Makefile`, then `make -C SRC`
- **Run**: pass parameters as `key=value` pairs (SU getpar style)

Example:

```bash
./bin/DEMEXE \
  meshfnm=data/mesh.txt iModel=1 fpeak=25 srcx=500 srcz=500 sstrength=1e6 \
  density=2500 vp=3000 vs=1732 nx=101 nz=101 dx=10 dz=10 xmin=0 zmin=0 \
  tmax=0.5 dt=0.001 mt=50 snapt=0.2
```

- Outputs:
  - `source.txt`: time vs source amplitude
  - `VX_<it>.txt`, `VZ_<it>.txt`: velocity fields at snapshot index `<it> = snapt/dt`

### What's Inside
- `docs/build.md`: dependencies and build details
- `docs/usage.md`: runtime parameters and output formats
- `docs/api/*`: API reference for functions and types
- `docs/examples/*`: mesh format and runnable examples

For API details, start with `docs/api/DemRun.md`.