#include "Element.h"
#include "DemRun.h"
#include "basefun.h"
#include "su.h"

void GetElemsCoord() // 获取所有单元顶点坐标信息
{
    afX = alloc1double(nNodes);
    afZ = alloc1double(nNodes);
    X1 = alloc1double(nElems);
    X2 = alloc1double(nElems);
    X3 = alloc1double(nElems);
    Z1 = alloc1double(nElems);
    Z2 = alloc1double(nElems);
    Z3 = alloc1double(nElems);

    for (int iNode = 0; iNode < nNodes; iNode++)
    {
        afX[iNode] = Node[iNode].NodeX;
        afZ[iNode] = Node[iNode].NodeZ;
    }

    for (int i = 0; i < nElems; i++)
    {
        X1[i] = afX[Elem[i].VertexId1 - 1]; // 节点编号从1开始，存储从0开始，所以要减一
        X2[i] = afX[Elem[i].VertexId2 - 1];
        X3[i] = afX[Elem[i].VertexId3 - 1];
        Z1[i] = afZ[Elem[i].VertexId1 - 1];
        Z2[i] = afZ[Elem[i].VertexId2 - 1];
        Z3[i] = afZ[Elem[i].VertexId3 - 1];
    }
    printf("\n");
    printf("Vertex Coordinate of First Element:\n");
    printf("ElemsVer1NodeX=%f, ElemsVer1NodeZ=%f\n", X1[0], Z1[0]);
    printf("ElemsVer2NodeX=%f, ElemsVer2NodeZ=%f\n", X2[0], Z2[0]);
    printf("ElemsVer3NodeX=%f, ElemsVer3NodeZ=%f\n", X3[0], Z3[0]);
    printf("\n");
    printf("===================GOT Elems Coord=================");
    printf("\n");
}

void GetTriNodes()
{
    for (int ii = 0; ii < nElems; ii++)
    {
        MyTri3[ii].VertexCoord[0][0] = X1[ii];
        MyTri3[ii].VertexCoord[0][1] = Z1[ii];
        MyTri3[ii].VertexCoord[1][0] = X2[ii];
        MyTri3[ii].VertexCoord[1][1] = Z2[ii];
        MyTri3[ii].VertexCoord[2][0] = X3[ii];
        MyTri3[ii].VertexCoord[2][1] = Z3[ii];

        MyTri3[ii].CenterCoord[0] = (X1[ii] + X2[ii] + X3[ii]) / 3;
        MyTri3[ii].CenterCoord[1] = (Z1[ii] + Z2[ii] + Z3[ii]) / 3;

        for (int i2 = 0; i2 < 3; i2++)
        {
            MyTri3[ii].m_aiNodes[i2] = i2;                    // 节点在单元内局部编号
            MyTri3[ii].m_aiEdges[i2] = i2;                    // 每个单元的三条边赋予局部编号
            MyTri3[ii].m_bgNodes[i2] = Elem[ii].VertexID[i2]; // 节点在离散前的全局编号
        }

        // 每条边有两个节点, 给这些节点逆时针局部编号
        MyTri3[ii].m_aiEdgeNodes[0][0] = 0;
        MyTri3[ii].m_aiEdgeNodes[0][1] = 1;
        MyTri3[ii].m_aiEdgeNodes[1][0] = 1;
        MyTri3[ii].m_aiEdgeNodes[1][1] = 2;
        MyTri3[ii].m_aiEdgeNodes[2][0] = 2;
        MyTri3[ii].m_aiEdgeNodes[2][1] = 0;
    }

    printf("LocalNodesNo=%d\n", MyTri3[0].m_aiNodes[1]);
    printf("GobalNodesNo=%d\n", MyTri3[0].m_bgNodes[1]);
    printf("LocalEdgesNo=%d\n", MyTri3[0].m_aiEdges[1]);
    printf("LocalEdgesNodesNo=%d\n", MyTri3[0].m_aiEdgeNodes[0][1]);
    printf("===================GOT Tri Nodes=================\n");
    printf("\n");
}

void BuildDps()
{
    m_nDns = 0;
    m_nDps = 0;
    int count = 0;

    for (int ii = 0; ii < nElems; ii++)
    {
        MyTri3[ii].m_iDnPos = m_nDns;
        MyTri3[ii].m_iDpPos = m_nDps;

        m_nDns++;
        m_nDps++;
    }

    fsxx = alloc1double(m_nDps);
    fszz = alloc1double(m_nDps);
    fsxz = alloc1double(m_nDps);
    fexx = alloc1double(m_nDps);
    fezz = alloc1double(m_nDps);
    fexz = alloc1double(m_nDps);

    for (int iDp = 0; iDp < m_nDps; iDp++)
    {
        fsxx[iDp] = fszz[iDp] = fsxz[iDp] = 0.0;
        fexx[iDp] = fezz[iDp] = fexz[iDp] = 0.0;
    }
}

void GetMax_MinElemLength()
{
    double Coord1[2], Coord2[2], TempLength, minLength, maxLength;

    for (int iElem = 0; iElem < nElems; iElem++)
    {
        int NodesCount = 3;
        for (int i = 0, j = NodesCount - 1; i < NodesCount; j = i++)
        {
            int iGlobalNodeNo1 = MyTri3[iElem].m_bgNodes[i];
            int iGlobalNodeNo2 = MyTri3[iElem].m_bgNodes[j];

            Coord1[0] = afX[iGlobalNodeNo1 - 1];
            Coord1[1] = afZ[iGlobalNodeNo1 - 1];
            Coord2[0] = afX[iGlobalNodeNo2 - 1];
            Coord2[1] = afZ[iGlobalNodeNo2 - 1];

            TempLength = TwoPointDistance(Coord1, Coord2);
            minLength = BegDem.m_pAux.m_fMinElemLength;
            maxLength = BegDem.m_pAux.m_fMaxElemLength;

            if (TempLength < minLength)
            {
                minLength = TempLength;
            }
            if (TempLength > maxLength)
            {
                maxLength = TempLength;
            }

            BegDem.m_pAux.m_fMinElemLength = minLength;
            BegDem.m_pAux.m_fMaxElemLength = maxLength;
        }
    }

    // 单元平均长度
    BegDem.m_pAux.m_fAveElemLength = (minLength + maxLength) / 2.0;

    printf("MinElemLength=%f\n", BegDem.m_pAux.m_fMinElemLength);
    printf("MaxElemLength=%f\n", BegDem.m_pAux.m_fMaxElemLength);
    printf("AveElemLength=%f\n", BegDem.m_pAux.m_fAveElemLength);
    printf("\n");
    printf("==================GOT AveElemLength==================\n");
    printf("\n");
}

void GetMax_MinNodeCoord()
{
    int NodesNum = BegDem.m_pAux.m_nNodes;
    for (int inode = 0; inode < NodesNum; inode++)
    {
        int ib = BegDem.m_pAux.m_aNodes[inode].m_aIdenticalNodes[0].m_iBlock;
        int in = BegDem.m_pAux.m_aNodes[inode].m_aIdenticalNodes[0].m_iNode;
        int iVarNodeNo = MyTri3[ib].m_aiNodes[in];

        if (BegDem.m_pAux.m_fMaxCoordX < afX[iVarNodeNo])
        {
            BegDem.m_pAux.m_fMaxCoordX = afX[iVarNodeNo];
        }
        if (BegDem.m_pAux.m_fMinCoordX > afX[iVarNodeNo])
        {
            BegDem.m_pAux.m_fMinCoordX = afX[iVarNodeNo];
        }
        if (BegDem.m_pAux.m_fMaxCoordZ < afZ[iVarNodeNo])
        {
            BegDem.m_pAux.m_fMaxCoordZ = afZ[iVarNodeNo];
        }
        if (BegDem.m_pAux.m_fMinCoordZ > afZ[iVarNodeNo])
        {
            BegDem.m_pAux.m_fMinCoordZ = afZ[iVarNodeNo];
        }
    }
    printf("NodeMinCoordX=%f\n", BegDem.m_pAux.m_fMinCoordX);
    printf("NodeMinCoordZ=%f\n", BegDem.m_pAux.m_fMinCoordZ);
    printf("NodeMaxCoordX=%f\n", BegDem.m_pAux.m_fMaxCoordX);
    printf("NodeMaxCoordZ=%f\n", BegDem.m_pAux.m_fMaxCoordZ);
    printf("\n");
    printf("==================GOT Max_MinNodeCoord==================\n");
}

int CalAreaTri3(int threadIdx) // 计算单元面积
{
    int i = threadIdx;
    float det;

    det = X1[i] * Z2[i] + X3[i] * Z1[i] + X2[i] * Z3[i] - X3[i] * Z2[i] - X2[i] * Z1[i] - X1[i] * Z3[i];
    ElemsArea[i] = 0.5 * det;

    return 0;
}

int MakeMatrixD(int threadIdx)
{
    MyElem *pElem = &MyTri3[threadIdx];
    for (int p = 0; p < 3; p++)
    {
        for (int q = 0; q < 3; q++)
        {
            pElem->MatrixD[p][q] = 0.0;
        }
    }

    double lamda, miu; // 拉梅系数

    lamda = young * poisson / ((1.0 + poisson) * (1.0 - 2.0 * poisson));
    miu = young / (2.0 * (1.0 + poisson));

    // 均匀各向同性完全弹性介质,设定弹性系数矩阵Cij(平面问题矩阵为3*3)
    pElem->MatrixD[0][0] = lamda + 2 * miu;
    pElem->MatrixD[1][1] = lamda + 2 * miu;
    pElem->MatrixD[2][2] = miu;

    pElem->MatrixD[0][1] = lamda;
    pElem->MatrixD[0][2] = 0;
    pElem->MatrixD[1][0] = lamda;
    pElem->MatrixD[1][2] = 0;
    pElem->MatrixD[2][0] = 0;
    pElem->MatrixD[2][1] = 0;

    return 0;
}

int GenerateMass(int threadIdx)
{
    CalAreaTri3(threadIdx);

    MyElem *pElem = &MyTri3[threadIdx];

    int i, j;
    double M[3][3] = {0.0}; // 质量矩阵

    // 质量矩阵生成: 运用三角形单元内基函数积分计算式（Ni*Nj=A/12）

    for (i = 0; i < 3; i++)
    {
        for (j = 0; j < 3; j++)
        {
            if (i == j)
            {
                M[i][j] = density * ElemsArea[threadIdx] / 6;
            }
            else
            {
                M[i][j] = density * ElemsArea[threadIdx] / 12;
            }
        }
    }

    /// 将质量矩阵各元素向对角线上聚集
    for (i = 0; i < 3; i++)
    {
        double value = 0.0;
        for (j = 0; j < 3; j++)
        {
            value += M[i][j];
            M[i][j] = 0.0;
        }
        // 为防止由于块体单元局部坐标系与整体坐标系反号所造成的体积负值，特加绝对值
        M[i][i] = fabs(value);
    }

    /// 将值传递给节点质量
    for (i = 0; i < 3; i++)
    {
        pElem->aNodes[i].m_fMass = M[i][i];
    }

    return 0;
}

int GenerateStiff(int threadIdx)
{
    MakeMatrixD(threadIdx);

    MyElem *pElem = &MyTri3[threadIdx];
    int ii = threadIdx;
    int i, j, k, l, inode;
    double x, z;
    double pNx, pNz; // 每个形函数的偏导数

    double w = 1.0 / 3.0;

    double N[3];             // 形函数
    double a[3], b[3], c[3]; // 坐标变换的辅助参数

    double pN[2][3] = {0.0};    // 形函数偏导矩阵
    double aNode[3][2] = {0.0}; // 单元坐标矩阵
    double B[3][6] = {0.0};     // 应变矩阵
    double B_T[6][3] = {0.0};   // 应变矩阵的转置

    a[0] = X2[ii] * Z3[ii] - X3[ii] * Z2[ii];
    a[1] = X3[ii] * Z1[ii] - X1[ii] * Z3[ii];
    a[2] = X1[ii] * Z2[ii] - X2[ii] * Z1[ii];

    b[0] = Z2[ii] - Z3[ii];
    b[1] = Z3[ii] - Z1[ii];
    b[2] = Z1[ii] - Z2[ii];

    c[0] = X3[ii] - X2[ii];
    c[1] = X1[ii] - X3[ii];
    c[2] = X2[ii] - X1[ii];

    for (k = 0; k < 3; k++)
    {
        // N[k] = (a[k] + b[k] * x + c[k] * z) / (2.0 * ElemsArea[ii]);    // 形函数

        pN[0][k] = b[k] / (2.0 * ElemsArea[ii]); // x偏导
        pN[1][k] = c[k] / (2.0 * ElemsArea[ii]); // z偏导
    }

    for (inode = 0; inode < 3; inode++)
    {
        aNode[inode][0] = afX[Elem[ii].VertexID[inode] - 1];
        aNode[inode][1] = afZ[Elem[ii].VertexID[inode] - 1];
    }

    double J[2][2] = {0.0};     // 雅可比矩阵
    double inv_J[2][2] = {0.0}; // 雅可比矩阵逆矩阵
    double det_J = 0.0;         // 雅可比矩阵行列式

    // 矩阵相乘，求雅克比矩阵
    for (i = 0; i < 2; i++)
    {
        for (l = 0; l < 2; l++)
        {
            for (j = 0; j < 3; j++)
            {
                J[i][l] += pN[i][j] * aNode[j][l];
            }
        }
    }

    /// 求雅可比矩阵行列式
    det_J = J[0][0] * J[1][1] - J[1][0] * J[0][1];

    pElem->DeltJacobi = det_J;

    /// 雅可比矩阵求逆
    inv_J[0][0] = J[1][1] / det_J;
    inv_J[0][1] = -J[0][1] / det_J;
    inv_J[1][0] = -J[1][0] / det_J;
    inv_J[1][1] = J[0][0] / det_J;

    // B矩阵生成
    for (k = 0; k < 3; k++)
    {
        pNx = pN[0][k];
        pNz = pN[1][k];

        B[0][0 + 2 * k] = pNx;
        B[0][1 + 2 * k] = 0.0;
        B[1][0 + 2 * k] = 0.0;
        B[1][1 + 2 * k] = pNz;
        B[2][0 + 2 * k] = pNz;
        B[2][1 + 2 * k] = pNx;
    }

    // 计算B矩阵的转置
    for (i = 0; i < 3; i++)
    {
        for (j = 0; j < 6; j++)
        {
            B_T[j][i] = B[i][j];
        }
    }

    for (i = 0; i < 3; i++)
    {
        for (j = 0; j < 6; j++)
        {
            pElem->MatrixB[i][j] = B[i][j];
        }
    }

    // 刚度矩阵生成
    for (i = 0; i < 6; i++)
    {
        for (j = 0; j < 6; j++)
        {
            pElem->Stiff[i][j] = 0.0; // 刚度矩阵清零

            for (k = 0; k < 3; k++)
            {
                for (l = 0; l < 3; l++)
                {
                    pElem->Stiff[i][j] += B_T[i][k] * pElem->MatrixD[k][l] * B[l][j] * ElemsArea[ii];
                }
            }
        }
    }

    /////////////////////////////////////////////////////////////////////////
    return 0;
}

int CalculateDeltaStress(int threadIdx) // 根据单元节点位移计算单元应变应力
{
    int i, j, k, inode;
    MyElem *pElem = &(MyTri3[threadIdx]);

    for (inode = 0; inode < 3; inode++) // 单元节点位移向量
    {
        CbNode *pNode = &(MyTri3[threadIdx].aNodes[inode]);
        pElem->DisplaceVector[inode * 2 + 0] = pNode->m_afDeltaDisplace[0];
        pElem->DisplaceVector[inode * 2 + 1] = pNode->m_afDeltaDisplace[1];
    }

    // 单元节点位移表示的单元应变
    for (i = 0; i < 3; i++)
    {
        for (j = 0; j < 6; j++)
        {
            pElem->DeltaStrain[i] += pElem->MatrixB[i][j] * pElem->DisplaceVector[j];
        }
    }

    // 单元应力
    for (i = 0; i < 3; i++)
    {
        for (k = 0; k < 3; k++)
        {
            pElem->DeltaStress[i] += pElem->MatrixD[i][k] * pElem->DeltaStrain[k];
        }
    }
    // printf("====================CalculateDeltaStress========================\n");
    return 0;
}

int CalculateNodeForce(int threadIdx) // 根据单元节点位移计算单元节点力
{
    MyElem *pElem = &(MyTri3[threadIdx]);
    // 单元节点力向量
    for (int i = 0; i < 6; i++)
    {
        for (int j = 0; j < 6; j++)
        {
            pElem->NodalForce[i] += pElem->Stiff[i][j] * pElem->DisplaceVector[j];
        }
    }

    // 将节点力向量赋值给单元节点
    for (int iNode = 0; iNode < 3; iNode++)
    {
        CbNode *pNode = &(MyTri3[threadIdx].aNodes[iNode]);
        pNode->m_afForce[0] = pElem->NodalForce[iNode * 2];     // X方向
        pNode->m_afForce[1] = pElem->NodalForce[iNode * 2 + 1]; // Z方向
    }

    for (int i = 0; i < 6; i++)
    {
        if (pElem->NodalForce[i] != 0)
            printf("threadIdx=%d, i=%d, NodalForce=%e\n", threadIdx, i, pElem->NodalForce[i]);
    }
    // printf("====================CalculateNodeForce========================\n");
    return 0;
}

int CalculateAverStressStrain(int threadIdx)
{
    MyElem *pElem = &MyTri3[threadIdx];
    double afAverStress[3] = {0.0};
    double afAverStrain[3] = {0.0};

    double **afStress, **afStrain;
    afStress = alloc2double(3, m_nDps);
    afStrain = alloc2double(3, m_nDps);

    for (int i2 = 0; i2 < m_nDps; i2++)
    {
        afStress[i2][0] = fsxx[i2];
        afStress[i2][1] = fszz[i2];
        afStress[i2][2] = fsxz[i2];

        afStrain[i2][0] = fexx[i2];
        afStrain[i2][1] = fezz[i2];
        afStrain[i2][2] = fexz[i2];
    }

    int nDps = 1;
    for (int i = 0; i < nDps; i++)
    {
        int iDpPos = pElem->m_iDpPos;
        for (int j = 0; j < 3; j++)
        {
            afAverStress[j] += afStress[iDpPos][j];
            afAverStrain[j] += afStrain[iDpPos][j];
        }
    }

    for (int i = 0; i < 3; i++)
    {
        pElem->AverStress[i] = afAverStress[i];
        pElem->AverStrain[i] = afAverStrain[i];
    }

    free2double(afStress);
    free2double(afStrain);

    // printf("====================CalculateAverStressStrain========================\n");
    return 0;
}

int CalculatePrinStressStrain(int threadIdx)
{
    MyElem *pElem = &MyTri3[threadIdx];
    double PrinStress[2] = {0.0};
    double PrinStrain[2] = {0.0};
    double lmStress[2][2] = {0.0};
    double lmStrain[2][2] = {0.0}; // 方向余弦

    Cal_Principal_Stress(pElem->AverStress, PrinStress, lmStress);
    Cal_Principal_Stress(pElem->AverStrain, PrinStrain, lmStrain);

    for (int i = 0; i < 2; i++)
    {
        pElem->PrinStress[i] = PrinStress[i];
        pElem->PrinStrain[i] = PrinStrain[i];

        for (int j = 0; j < 2; j++)
        {
            pElem->PrinStrainDir[i][j] = lmStress[i][j];
            pElem->PrinStrainDir[i][j] = lmStrain[i][j];
        }
    }
    // printf("====================CalculatePrinStressStrain========================\n");
    return 0;
}