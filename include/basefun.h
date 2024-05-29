#include <stdio.h>
#ifndef __BASEFUN_H__
#define __BASEFUN_H__

int SetPositionByLine(FILE *fp, int nLine);
void norm(double w[2]);
void normalize_matrix(double **matrix, int rows, int cols);
double TwoPointDistance(double point1[2], double point2[2]);
double Cal_Principal_Stress(double stress_strain[3], double principal[2], double lm[2][2]);
double ricker(double t, double fpeak);

#endif