## Element.h API

Header: `include/Element.h`

### Types
- `MyElem`
  - Element state for linear triangle with 3 nodes
  - Key fields:
    - `CbNode aNodes[3]`: nodal physical state (forces, displacement increments, etc.)
    - `double Stiff[6][6]`, `MatrixB[3][6]`, `MatrixD[3][3]`, `DeltJacobi`
    - `double DisplaceVector[6]`, `NodalForce[6]`
    - `double DeltaStress[3]`, `DeltaStrain[3]`
    - `double AverStress[3]`, `AverStrain[3]`
    - `double PrinStress[2]`, `PrinStrain[2]`, directions
    - `double VertexCoord[3][2]`, `CenterCoord[2]`
    - `int m_iElemModel` (1 = linear)
    - Node/edge mappings: `m_aiNodes[3]`, `m_bgNodes[3]`, `m_agNodes[3]`, `m_aiEdges[3]`, `m_aiEdgeNodes[3][2]`
    - Discrete indices: `m_iDnPos`, `m_iDpPos`
- `MyElem *MyTri3` global array of size `nElems`

### Functions
- `void GetElemsCoord()`
  - Load node coordinates into element-wise arrays `X1..Z3`
- `void GetTriNodes()`
  - Populate local node/edge indices and centers for each element
- `void BuildDps()`
  - Build discrete points (Gauss/integration points) and allocate stress/strain arrays
- `void GetMax_MinElemLength()`
  - Compute min/max/avg edge length across mesh
- `void GetMax_MinNodeCoord()`
  - Compute min/max coordinates across nodes
- `int CalAreaTri3(int threadIdx)`
  - Compute area for element `threadIdx`
- `int MakeMatrixD(int threadIdx)`
  - Build constitutive matrix (linear elastic)
- `int GenerateMass(int threadIdx)`
  - Assemble lumped mass for nodal DOFs
- `int GenerateStiff(int threadIdx)`
  - Assemble element stiffness
- `int CalculateDeltaStress(int threadIdx)`
  - Compute stress/strain increments from displacement increments
- `int CalculateNodeForce(int threadIdx)`
  - Internal nodal force from stress state
- `int CalculateAverStressStrain(int threadIdx)`
  - Average element stress/strain
- `int CalculatePrinStressStrain(int threadIdx)`
  - Principal stress/strain and directions

### Notes
- Many functions accept an element index `threadIdx` to enable OpenMP parallelism.
- See `SRC/Element.c` for implementation details.