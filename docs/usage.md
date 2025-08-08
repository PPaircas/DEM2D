## Usage

Run the executable with SU-style `key=value` parameters.

### Required Parameters
- **meshfnm**: path to mesh file
- **tmax**: simulation duration (s)

### Common Parameters
- **iModel**: model type (1 = linear)
- **fpeak**: Ricker peak frequency (Hz)
- **srcx, srcz**: source coordinates (same units as mesh)
- **sstrength**: source amplitude scale
- **dt**: time step (s); if not set, defaults to 1e-3
- **nt**: number of steps; if not set, computed from `tmax/dt + 1`
- **mt**: print/monitor interval (default 100)
- **density, vp, vs**: material properties; `young` and `poisson` derived internally
- **nx, nz**: output grid samples in X/Z
- **dx, dz**: output grid spacing in X/Z
- **xmin, zmin**: origin of output grid
- **snapt**: time (s) to snapshot velocity fields (outputs `VX_*.txt`, `VZ_*.txt`)

Example:
```bash
./bin/DEMEXE meshfnm=data/mesh.txt iModel=1 fpeak=25 srcx=500 srcz=500 sstrength=1e6 \
  density=2500 vp=3000 vs=1732 nx=101 nz=101 dx=10 dz=10 xmin=0 zmin=0 \
  tmax=0.5 dt=0.001 mt=50 snapt=0.2
```

### Mesh File Format
Text file parsed by `GetMesh()`:
1. Line 1: `nNodes nElems`
2. Next `nNodes` lines: `NodeId NodeX NodeZ <ignored>`
3. Next `nElems` lines: `ElemId <ignored> <ignored> VertexId1 VertexId2 VertexId3`

See `docs/examples/minimal-mesh.md` for a concrete example.

### Outputs
- `source.txt`: two columns, time and source value per step
- `VX_<it>.txt`, `VZ_<it>.txt`: grid of `nx` rows by `nz` columns, ASCII floats per line at snapshot `it = snapt/dt`
  - Files are written from `OutputFile()` in `SRC/DemRun.c`.

### Execution Flow
1) `LoadParameters()` reads CLI parameters
2) `GetSize()` computes counts from mesh header
3) `InitializeVar()` allocates and initializes state
4) `GetMesh()` loads nodes/elements and builds discrete nodes
5) `SetElementArguments()`, `SetContactArguments()` choose models/materials
6) `Preprocess()` builds mass/stiffness and contact spring constants
7) `Calculation()` loops in time calling `CalcOneStep()` and writes outputs

For API details, see `docs/api/DemRun.md`.