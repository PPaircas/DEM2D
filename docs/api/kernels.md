## Kernel Functions

### Element Kernels (`SRC/ElementKernel.c`)
- `int DoNode_Inc(int threadIdx)`
  - For element `threadIdx`, update nodal acceleration, velocity, displacement increments using external/internal forces and `dt`.
  - Writes to global arrays `fA`, `fVV`, `fUU`, `fUX/fUZ`, `fVX/fVZ`.
- `int DoBlock_B_Inc(int threadIdx)`
  - For element `threadIdx`, runs element pipeline for model 1 (linear):
    - `DoNode_Inc -> CalculateDeltaStress -> CalculateNodeForce`

### Contact Kernels (`SRC/ContactKernel.c`)
- `int DoContact_Linear_Inc(int threadIdx)`
  - For contact `threadIdx`, compute spring forces in local (tangent/normal) system from nodal displacement increments and map to global forces on the two elements
- `int SyncSpringGroup(int threadIdx)`
  - For spring force group `threadIdx`, accumulate spring forces to associated nodal internal forces `m_afForce`

These kernels are called in `CalcOneStep()` in the order: elements -> contacts -> spring-force sync.