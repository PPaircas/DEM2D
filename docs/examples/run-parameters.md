## Runtime Parameters Reference

- meshfnm: string (required)
- iModel: int, element/contact model
  - 1 = linear (default if provided in examples)
- fpeak: double, Ricker peak frequency (Hz)
- srcx, srcz: double, source coordinates
- sstrength: double, source amplitude scale
- density: double (kg/m^3)
- vp, vs: double (m/s)
- tmax: double, total time (s) [required]
- dt: double, time step (s). Default: 1e-3 if not set
- nt: int, total steps. Default: computed from tmax/dt if not set
- mt: int, monitor interval. Default: 100
- nx, nz: int, output grid samples
- dx, dz: double, output grid spacing
- xmin, zmin: double, output grid origin. Default: 0 if not set
- snapt: double, time (s) for snapshot output

Notes:
- Internally, `young` and `poisson` are derived from `density, vp, vs`:
  - `poisson = ( (vp/vs)^2 - 2 ) / ( 2*(vp/vs)^2 - 2 )`
  - `young = 2 * density * vs^2 * (1 + poisson)`
- Snapshot index: `it = snapt / dt`