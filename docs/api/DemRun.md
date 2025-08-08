## DemRun.h API

Header: `include/DemRun.h`

### Globals (selected)
- Simulation/grid: `nx, nz, dx, dz, xmin, xmax, zmin, zmax, pX[], pZ[]`
- Time stepping: `t, tmax, dt, it, nt, mt, snapt`
- Source: `fpeak, srcx, srcz, sstrength, source[]`
- Material (bulk): `density, vp, vs, young, poisson, m_oBulkMaterial`
- Mesh sizes: `nElems, nNodes, dNodes, m_nDns, m_nDps`
- Coordinates: `afX, afZ` (original), `dfX, dfZ` (discrete)
- Kinematics: `fU,fV,fA` current; `fUU,fVV` next; `fUX,fUZ,fVX,fVZ`
- Stress/strain (at DPs and nodes): `fsxx, fszz, fsxz, fexx, fezz, fexz`, `fNsxx...`, `fNexx...`
- Forces: `m_fExtForce, m_fUnbalForce`
- Source location: `sElemId, sNodeId`

### Types (selected)
- `Point { float x,z; int GlobalNodeNo; }`
- `NodeInfo { int NodeId; float NodeX, NodeZ; }`
- `ElemInfo { int ElemId; int VertexId1, VertexId2, VertexId3; int VertexID[3]; }`
- `CbNode`
  - Nodal mass `m_fMass`, fixed flags `m_abFixU[2]`, delta displacement `m_afDeltaDisplace[2]`, external/internal force `m_afExtForce[2]`, `m_afForce[2]`, principal stress/strain and directions
- `CbNodePos { int m_iBlock; int m_iNode; }`
- `CbNodeGroup { CbNodePos *m_aNodePos; }`
- `CbForcePos { int m_iBlock; int m_iSpring, m_iForce; }`
- `CbNodeAux { int m_aiGlobalNodeNo; CbNodePos *m_aIdenticalNodes; }`
- `CbFaceAux { bool isContact; int m_aiGlobalNodeNo[2]; CbNodePos *m_aIdenticalFaces; }`
- `CbDemAux`
  - Global counts, bounding box, global faces/edges, contact counts, mapping arrays
- `CbBulkMaterial { double fDensity, fYoung, fPoisson; }`
- `CbSpring`
  - Node linkage, used/broken flags, area, normal/shear stiffness, length, displacement, force vectors
- `CbSpringForceGroup`
  - Node position bound to multiple spring forces
- `CbDem`
  - Holds `CbDemAux`, groups, and contact/touch stiffness factors. Global instance `BegDem`.
- `CbContact`
  - Two elements forming a contact face, face index, local coord system (tangent/normal), material and springs. Global array `MyContact`.

### Initialization Helpers (static)
- `Point_Init`, `CbNode_Init`, `CbNodePos_Init`, `CbDemAux_Init`, `CbFaceAux_Init`, `CbBulkMaterial_Init`, `CbSpringForceGroup_Init`, `CbDem_Init`, `CbContact_Init`

### Material Helpers (static)
- `int SetMaterial(CbContact *pContact, CbBulkMaterial oMaterial)`
- `CbBulkMaterial GetMaterial(CbContact pContact)`

### Orchestration Functions
- `void LoadParameters()`
  - Read CLI parameters via SU getpar, compute `nt` if omitted, print banner
- `void GetSize()`
  - Read first mesh line to set `nNodes`, `nElems`, set `dNodes = nElems*3`
- `void InitializeVar()`
  - Initialize materials (`young`, `poisson`), allocate all global arrays, zero state
- `void GetMesh()`
  - Load nodes/elements, populate `MyTri3` mappings, build discrete nodes `dfX/dfZ`, set counts in `BegDem.m_pAux`, call `FindIdenticalNodes/FindIdenticalFaces`, `SetAllContacts()`, and compute geometry extents
- `void SetElementArguments(int iModel, CbBulkMaterial m_oBulkMaterial)`
  - Assign element model/material to each element
- `void SetContactArguments(int iModel, CbBulkMaterial m_oBulkMaterial)`
  - Assign contact models and materials (homogeneous vs derived)
- `void Preprocess()`
  - Per-element: `GenerateMass`, `GenerateStiff`; per-contact: `Build()` spring stiffness; locate source node
- `void LocateSourceByCoord(double srcx, double srcz)`
  - Nearest node to `(srcx,srcz)` across all elements -> `sElemId`, `sNodeId`
- `void ApplySource(double t, double fpeak, double *source, double sstrength)`
  - Ricker wavelet delayed by `1/fpeak`; apply as Z-force at `(sElemId,sNodeId)`
- `void Calculation()`
  - Time loop: apply source, `CalcOneStep()`, store nodal fields, write snapshot at `snapt`, write `source.txt`
- `void CalcOneStep()`
  - Parallel over elements: `DoBlock_B_Inc` (kinematics, stresses, forces)
  - Parallel over contacts: `DoContact_Linear_Inc`
  - Parallel over spring-force groups: `SyncSpringGroup`
- `void FindAdjacentPoints(double pX, double pZ, int radius)`
  - Gather discrete nodes within a radius and compute distance weights
- `void GetWavefield(...)` and `void OutputFile(...)`
  - Interpolate discrete nodal fields to regular `nx x nz` grid and write ASCII files
- `void SetVelocityBySel()`, `void UpdateNodeStressStrain()`, `void UpdateNodePrinStressStrain()`
  - Boundary condition placeholder and stress post-processing
- `void FreeVar()`
  - Release allocated resources (see implementation)

### Kernel Functions
- `int DoNode_Inc(int threadIdx)`, `int DoBlock_B_Inc(int threadIdx)` in `ElementKernel.c`
- `int DoContact_Linear_Inc(int threadIdx)`, `int SyncSpringGroup(int threadIdx)` in `ContactKernel.c`

### Contact Setup
- `int SetSpringNodesGroup()`
- `int SetAllContacts()`
- `int Build(int iContact, double fContactStiffFactor)`

### Example Initialization Flow
```c
LoadParameters();
GetSize();
InitializeVar();
GetMesh();
SetElementArguments(iModel, m_oBulkMaterial);
SetContactArguments(iModel, m_oBulkMaterial);
Preprocess();
Calculation();
FreeVar();
```