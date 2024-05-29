#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include "../include/basefun.h"

#define BUFSIZE 1024
#define PI 3.14159265

int SetPositionByLine(FILE *fp, int nLine)
{
    int i = 0;
    char buffer[BUFSIZE + 1];
    fpos_t pos;

    rewind(fp);
    for (; i < nLine; i++)
        fgets(buffer, BUFSIZE, fp);
    fgetpos(fp, &pos);
    return 0;
}

void norm(double w[2])  ///////////////////向量单位化
{
    int i;
    double size;

    size = (double)sqrt(w[0] * w[0] + w[1] * w[1]);
    for (i = 0; i < 2; i++)
    {
        w[i] = w[i] / size;
    }
}

//////////////////计算两点距离
double TwoPointDistance(double point1[2], double point2[2])
{
    double dx = point2[0] - point1[0];
    double dz = point2[1] - point1[1];
    double dist = sqrt(dx * dx + dz * dz);
    return (dist);
}

// 计算二维平面的主应力及方向余弦
double Cal_Principal_Stress(double stress_strain[3], double principal[2], double lm[2][2])
{
    double temp1 = (stress_strain[0] + stress_strain[1]) / 2.0;
    double temp2 = (stress_strain[0] - stress_strain[1]) / 2.0;
    double temp31 = pow(temp2, 2.0);
    double temp32 = pow(stress_strain[3], 2.0);

    principal[0] = temp1 + sqrt(temp31 + temp32);
    principal[1] = temp1 - sqrt(temp31 + temp32);

    double mvl[2] = {0.0};

    for (int i = 0; i < 2; i++)
    {
        mvl[i] = (principal[i] - (stress_strain[0] + stress_strain[2])) / (stress_strain[1] + stress_strain[2] - principal[i]);

        lm[0][i] = sqrt(1.0 / (1.0 + pow(mvl[i], 2.0)));
        lm[1][i] = sqrt(1.0 / (1.0 + pow(1.0 / mvl[i], 2.0)));
    }

    return 0;
}

void normalize_matrix(double **matrix, int rows, int cols)
{
    double min = matrix[0][0];
    double max = matrix[0][0];

    // 找到数组中的最小值和最大值
    for (int i = 0; i < rows; i++)
    {
        for (int j = 0; j < cols; j++)
        {
            if (matrix[i][j] < min)
            {
                min = matrix[i][j];
            }
            if (matrix[i][j] > max)
            {
                max = matrix[i][j];
            }
        }
    }

    // 对数组中的每个元素进行归一化
    for (int i = 0; i < rows; i++)
    {
        for (int j = 0; j < cols; j++)
        {
            matrix[i][j] = -1.0 + 2.0 * (matrix[i][j] - min) / (max - min);
        }
    }
}

// 雷克子波
double ricker(double t, double fpeak)
{
    double x, xx;
    x = PI * fpeak * t;
    xx = x * x;
    return exp(-xx) * (1.0 - 2.0 * xx);
}

// 写文本文件


// 写二进制文件
