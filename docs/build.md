## Build and Dependencies

### Prerequisites
- GCC with OpenMP (`-fopenmp`)
- Seismic Unix (SU) installed and built
  - Required libs: `-lsu -lpar -lcwp -lm`
  - Set `SU_HOME` to SU install root containing `include/` and `lib/`
- POSIX environment (Linux)

### Directory Layout
- `lib/`: builds `libmyfun.a` from `basefun.c`
- `SRC/`: builds the main executable `bin/DEMEXE`

### Build Steps
1) Build internal library
```bash
make -C lib
```
This produces `lib/libmyfun.a`.

2) Configure SU paths
Edit `SRC/Makefile` or export env vars:
- `SU_HOME=/path/to/su` (contains `include/` and `lib/`)
- Optionally `SU_HOME1` if you have extra Complex modules

3) Build executable
```bash
make -C SRC
```
This produces `bin/DEMEXE`.

4) Clean
```bash
make -C SRC clean
make -C lib clean
```

### Makefile Flags (SRC/Makefile)
- Includes: `-I$(SU_HOME)/include -I$(JOB_HOME)/include`
- Libs: `-L$(SU_HOME)/lib -L$(JOB_HOME)/lib -lmyfun -lsu -lpar -lcwp -lm`

Ensure `bin/` exists or is created by the build; output target is `../bin/DEMEXE` from `SRC`.