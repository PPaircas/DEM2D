## Contact API

Header prototypes: `include/DemRun.h`
Implementation: `SRC/Contact.c`

### Overview
Contact faces are automatically detected between elements that share an edge. Each contact face has two springs (at the two edge nodes) with normal and shear stiffness. A local (tangent/normal) coordinate system is built per contact.

### Functions
- `int SetSpringNodesGroup()`
  - Build `BegDem.m_aSpringNodeGroup` mapping identical global node positions to lists of `(block, localNode)` pairs across neighboring elements.
- `int SetAllContacts()`
  - Allocate and initialize `MyContact` for all detected contact faces
  - For each contact:
    - Identify paired nodes across the two adjacent elements (by matching global node ids)
    - Initialize springs, local coordinate system (tangent, normal) based on edge geometry
  - Build `BegDem.m_aSpringForceGroup` to map per-node spring forces for accumulation
- `int Build(int iContact, double fContactStiffFactor)`
  - For contact `iContact`, compute normal/shear stiffness for each spring:
    - `kn = Area * young / L`
    - `ks = Area * young / (L * (1 + poisson))`
    - where `L` derived from total area and `fContactStiffFactor`
  - Reset spring force/displacement state

### Kernels (see `docs/api/kernels.md`)
- `int DoContact_Linear_Inc(int threadIdx)`
  - Compute spring forces in contact-local coordinates from displacement increments; map to global forces
- `int SyncSpringGroup(int threadIdx)`
  - Accumulate mapped spring forces into nodal internal forces

### Data Structures
- `CbContact`
  - `CbNodePos m_iElem1, m_iElem2;`
  - `int m_iGlobalFaceNo;`
  - `double LocalCoordSystem[2][2];`
  - `CbBulkMaterial m_oBulkMaterial;`
  - `bool m_bUsed, m_homogeneous;`
  - `int m_iContactModel;`
  - `double fContactStiffFactor;`
  - `CbSpring m_aoSpring[2];`

- `CbSpring`
  - `Brother1NodeNo`, `Brother2NodeNo`
  - `Brother1GlobalNodeNo`, `Brother2GlobalNodeNo`
  - `Brother1SpringNodeNo`, `Brother2SpringNodeNo`
  - `bool m_bUsed, m_bBroken, m_bHistoryBroken`
  - `double Area, NormalStiff, ShearStiff, m_fSpringLength`
  - `double SpringForce[2], Displace[2], rela_dis_local[2]`
  - `double m_afEForce[2][2]` (mapped forces per side)

### Notes
- Contact detection and grouping rely on `FindIdenticalFaces()` and `BegDem.m_pAux` prepared in `GetMesh()`.
- The model id `iModel` controls the contact constitutive model; current kernels implement linear (id = 1).