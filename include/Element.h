#ifndef __ELEMENT_H__
#define __ELEMENT_H__

#include "DemRun.h"

typedef struct MyElemStruct
{
    CbNode aNodes[3];               // 存放单元节点的物理量（应力应变、内力外力等）
    double m_afFaceNodeAreas[3][2]; // 单元每条边上弹簧面积矩阵

    int m_aiEdges[3];        // 单元上每条边的局部编号
    int m_aiEdgeNodes[3][2]; // 单元每条边上的局部节点索引

    double MatrixB[3][6];
    double DeltJacobi;
    double MatrixD[3][3];
    double Stiff[6][6];

    double DisplaceVector[6]; // 单元所有节点位移向量
    double NodalForce[6];     // 单元上的节点力向量

    double DeltaStress[3]; // 单元应力
    double DeltaStrain[3]; // 单元应变
    double AverStress[3];
    double AverStrain[3];

    double PrinStress[2]; // PrinStress[0]-max; PrinStress[1]-min
    double PrinStrain[2];
    double PrinStressDir[2][2];
    double PrinStrainDir[2][2];

    double VertexCoord[3][2]; // 顶点坐标
    double CenterCoord[2];    // 单元中心点坐标 CenterCoord[0]-x; CenterCoord[1]-z

    int m_iElemModel; // 0-none 1-linear 2-MC
    CbBulkMaterial m_oMaterial;

    int m_aiNodes[3]; // 存储单元节点局部编号
    int m_bgNodes[3]; // 存储单元离散前的节点全局编号
    int m_agNodes[3]; // 存储单元离散后的节点全局编号
    int m_iDnPos;     // Discrete Nodes. Each element has a number of nodes which are not shared by other elements.
    int m_iDpPos;     // Discrete points. Each element has a given number of points that can be unequal to the number of nodes.

} MyElem;

MyElem *MyTri3;

//===========================================================
void GetElemsCoord();
void GetTriNodes();
void BuildDps();
void GetMax_MinElemLength();
void GetMax_MinNodeCoord();
int CalAreaTri3(int threadIdx);
int MakeMatrixD(int threadIdx);
int GenerateMass(int threadIdx);
int GenerateStiff(int threadIdx);
int CalculateDeltaStress(int threadIdx);
int CalculateNodeForce(int threadIdx);
int CalculateAverStressStrain(int threadIdx);
int CalculatePrinStressStrain(int threadIdx);

#endif