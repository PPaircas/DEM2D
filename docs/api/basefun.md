## basefun.h API

Header: `include/basefun.h`
Library: built into `lib/libmyfun.a`

### Functions
- `int SetPositionByLine(FILE *fp, int nLine)`
  - Rewinds and advances file pointer to the beginning of `nLine` (1-based)
- `void norm(double w[2])`
  - In-place normalize a 2D vector
- `void normalize_matrix(double **matrix, int rows, int cols)`
  - Min-max normalize matrix values to [-1, 1]
- `double TwoPointDistance(double point1[2], double point2[2])`
  - Euclidean distance between (x,z) points
- `double Cal_Principal_Stress(double stress_strain[3], double principal[2], double lm[2][2])`
  - 2D principal stress magnitudes and direction cosines from stress components
- `double ricker(double t, double fpeak)`
  - Ricker wavelet value at time `t` and peak frequency `fpeak`

### Usage Example
```c
#include "basefun.h"

double p1[2] = {0,0}, p2[2] = {3,4};
double d = TwoPointDistance(p1, p2); // 5.0

double v[2] = {3,4};
norm(v); // v -> {0.6, 0.8}

double s[3] = {sx, sz, sxz};
double princ[2], lm[2][2];
Cal_Principal_Stress(s, princ, lm);

double a = ricker(0.01, 25.0);
```