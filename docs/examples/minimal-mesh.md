## Minimal Mesh Example

A tiny 2-element mesh illustrating the expected file format.

```
# nNodes nElems
4 2
# NodeId  X     Z   (extra ignored)
1  0.0   0.0   0
2  10.0  0.0   0
3  0.0   10.0  0
4  10.0  10.0  0
# ElemId  ig1 ig2  v1 v2 v3  (ig* ignored)
1  0   0   1  2  3
2  0   0   2  4  3
```

Save as `data/mesh.txt` and run:

```bash
./bin/DEMEXE \
  meshfnm=data/mesh.txt iModel=1 fpeak=25 srcx=5 srcz=5 sstrength=1e6 \
  density=2500 vp=3000 vs=1732 nx=21 nz=21 dx=0.5 dz=0.5 xmin=0 zmin=0 \
  tmax=0.2 dt=0.001 snapt=0.1
```

Outputs:
- `source.txt`
- `VX_100.txt` and `VZ_100.txt` (since `snapt/dt = 0.1/0.001 = 100`)