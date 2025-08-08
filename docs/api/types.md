## Core Types

Defined in `include/DemRun.h` and `include/Element.h`.

### Node and Element Topology
- `typedef struct NodeInfo { int NodeId; float NodeX, NodeZ; } NodeInfo;`
- `typedef struct ElemInfo { int ElemId; int VertexId1, VertexId2, VertexId3; int VertexID[3]; } ElemInfo;`
- `typedef struct Point { float x, z; int GlobalNodeNo; } Point;`

### DEM State
- `typedef struct CbNode`
  - `int m_iBrokenType;`
  - `double m_fMass;`
  - `bool m_abFixU[2];`
  - `double m_afDeltaDisplace[2];`
  - `double m_afExtForce[2];`
  - `double m_afForce[2];`
  - `double m_afPrinStress[2];`
  - `double m_afPrinStrain[2];`
  - `double m_afPrinStressDir[2][2];`
  - `double m_afPrinStrainDir[2][2];`

- `typedef struct CbNodePos { int m_iBlock; int m_iNode; } CbNodePos;`
- `typedef struct CbNodeGroup { CbNodePos *m_aNodePos; } CbNodeGroup;`
- `typedef struct CbForcePos { int m_iBlock; int m_iSpring, m_iForce; } CbForcePos;`

### DEM Aux/Topology
- `typedef struct CbNodeAux { int m_aiGlobalNodeNo; CbNodePos *m_aIdenticalNodes; } CbNodeAux;`
- `typedef struct CbFaceAux { bool isContact; int m_aiGlobalNodeNo[2]; CbNodePos *m_aIdenticalFaces; } CbFaceAux;`
- `typedef struct CbDemAux`
  - counts: `m_nNodes, m_nElems, m_nVarNodes, m_nInitialNodes, m_nInitialElems`
  - contacts/faces/edges: `m_nContacts, m_nGlobalFaces, m_nGlobalEdges, m_nBrokenFaces`
  - geometry: `m_fMax/Min/AveElemLength`, `m_fMax/MinCoordX/Z`
  - arrays: `m_aNodes, m_aVarNodes, m_aDps, m_aGlobalFaces, m_aAuxNodeNoFromVarNodeNo`

### Material and Contact
- `typedef struct CbBulkMaterial { double fDensity, fYoung, fPoisson; }`
- `typedef struct CbSpring`
  - node mapping, flags, `Area`, `NormalStiff`, `ShearStiff`, `m_fSpringLength`
  - `double rela_dis_local[2];`
  - `double SpringForce[2];`
  - `double Displace[2];`
  - `double m_afEForce[2][2]; // [side][x/z]`

- `typedef struct CbSpringForceGroup`
  - `int m_nSprings;`
  - `CbNodePos m_oNodePos;`
  - `CbForcePos m_aoSpringForces[2];`

- `typedef struct CbDem`
  - `CbDemAux m_pAux;`
  - `CbNodeGroup *m_aSpringNodeGroup;`
  - `CbSpringForceGroup *m_aSpringForceGroup;`
  - `double m_fContactStiffFactor, m_fTouchStiffFactor;`

- `typedef struct CbContact`
  - `CbNodePos m_iElem1, m_iElem2;`
  - `int m_iGlobalFaceNo;`
  - `double LocalCoordSystem[2][2]; // tangent/normal`
  - `CbBulkMaterial m_oBulkMaterial;`
  - `bool m_bUsed, m_homogeneous;`
  - `int m_iContactModel;`
  - `double fContactStiffFactor;`
  - `CbSpring m_aoSpring[2];`

### Element State
- `typedef struct MyElem`
  - See `docs/api/Element.md` for detailed fields.

Global instances:
- `MyElem *MyTri3;`
- `CbDem BegDem;`
- `CbContact *MyContact;`