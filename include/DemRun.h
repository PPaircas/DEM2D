#ifndef __DEMRUN_H__
#define __DEMRUN_H__

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdbool.h>

#define BUFSIZE 1024
#define MAX_SHARE_ELEMENTS 10

//========================================================parameters required

char buf[BUFSIZE];
FILE *fp;
char *meshfnm;
int iModel;
double fpeak;	   // 子波峰值频率
double srcx, srcz; // 震源坐标
double t;
double tmax; // 计算的最大时长
double dt;	 // 时间步长
double *source;
double sstrength; // 震源振幅
int it, nt, mt;
double snapt; // 输出波场的时刻
double vp, vs;
double density; // 密度
double young;	// 杨氏模量
double poisson; // 泊松比

int nx, nz;
double dx, dz;
double xmin, xmax, zmin, zmax;
double *pX, *pZ; // 模型网格点坐标

//========================================================define variables

int nElems; //	网格单元数量
int nNodes; //	单元离散前网格节点数
int dNodes; //  单元离散后网格节点数
int m_nDns; //	节点积分点数量
int m_nDps; //	高斯积分点数量

double *afX, *X1, *X2, *X3; // 所有单元节点的X坐标
double *afZ, *Z1, *Z2, *Z3; // 所有单元节点的Z坐标
double *dfX, *dfZ;			// 所有离散单元节点的坐标

double *ElemsArea; // 所有单元面积

double **fU, **fV, **fA;	   // 当前时刻的位移、速度、加速度
double **fUU, **fVV;		   // 下一时刻的位移、速度、加速度
double *fUX, *fUZ;			   // 离散单元节点位移分量
double *fVX, *fVZ;			   // 离散单元节点速度分量
double *fsxx, *fszz, *fsxz;	   // 高斯点应力
double *fexx, *fezz, *fexz;	   // 高斯点应变
double *fNsxx, *fNszz, *fNsxz; // 离散单元节点应力
double *fNexx, *fNezz, *fNexz; // 离散单元节点应变

double m_fExtForce;	  // 所有单元的外力绝对值之和
double m_fUnbalForce; // 所有单元的不平衡力绝对值之和

int sElemId, sNodeId; // 施加震源的点，其所属的单元号及节点局部编号

//=======================================================define structs

typedef struct PointStruct
{
	float x;
	float z;
	int GlobalNodeNo;
} Point;

static void Point_Init(Point *pPoint)
{
	pPoint->x = -1;
	pPoint->z = -1;
	pPoint->GlobalNodeNo = -1;
}

Point *adjacentPoints;
double *weights;
int numAdjacentPoints;

typedef struct NodeInfoStruct
{
	int NodeId;
	float NodeX;
	float NodeZ;
} NodeInfo;

NodeInfo *Node;

typedef struct ElemInfoStruct
{
	int ElemId;
	int VertexId1;
	int VertexId2;
	int VertexId3;
	int VertexID[3];
} ElemInfo;

ElemInfo *Elem;

typedef struct CbNodeStruct
{
	int m_iBrokenType; // 0-未破坏 1-拉破坏 2-剪破坏
	double m_fMass;
	bool m_abFixU[2]; // 节点是否固定，施加边界条件时，实现将选中的点进行固定

	double m_afDeltaDisplace[2]; // [0]-X方向位移；[1]-Z方向位移
	double m_afExtForce[2];		 // [0]-X方向合力；[1]-Z方向合力
	double m_afForce[2];		 // [0]-X方向合力；[1]-Z方向合力

	double m_afPrinStress[2];
	double m_afPrinStrain[2];
	double m_afPrinStressDir[2][2];
	double m_afPrinStrainDir[2][2];
} CbNode;

static void CbNode_Init(CbNode *pNode)
{
	pNode->m_iBrokenType = 0;
	pNode->m_fMass = 0;
	pNode->m_abFixU[0] = false;
	pNode->m_abFixU[1] = false;
	pNode->m_afDeltaDisplace[0] = 0;
	pNode->m_afDeltaDisplace[1] = 0;
	pNode->m_afExtForce[0] = 0;
	pNode->m_afExtForce[1] = 0;
	pNode->m_afForce[0] = 0;
	pNode->m_afForce[1] = 0;
}

typedef struct CbNodePosStruct
{
	int m_iBlock;
	int m_iNode;
} CbNodePos;

static void CbNodePos_Init(CbNodePos *pNode)
{
	pNode->m_iBlock = -1;
	pNode->m_iNode = -1;
}

typedef struct CbNodeGroupStruct
{
	CbNodePos *m_aNodePos;
} CbNodeGroup;

typedef struct CbForcePosStruct
{
	int m_iBlock;
	int m_iSpring, m_iForce;
} CbForcePos;

//=======================================================Aux 拓扑类

typedef struct CbNodeAuxStruct
{
	int m_aiGlobalNodeNo;
	CbNodePos *m_aIdenticalNodes; // block_id, node_id
} CbNodeAux;

typedef struct CbFaceAuxStruct
{
	bool isContact;
	int m_aiGlobalNodeNo[2];
	CbNodePos *m_aIdenticalFaces; // block_id, edge_id
} CbFaceAux;

static void CbFaceAux_Init(CbFaceAux *pFace)
{
	pFace->isContact = false;
	pFace->m_aiGlobalNodeNo[0] = 0;
	pFace->m_aiGlobalNodeNo[1] = 0;
}

static bool isEqual(CbFaceAux *left, CbFaceAux *right)
{
	int a0 = left->m_aiGlobalNodeNo[0];
	int a1 = left->m_aiGlobalNodeNo[1];

	int b0 = right->m_aiGlobalNodeNo[0];
	int b1 = right->m_aiGlobalNodeNo[1];

	return (a0 == b0 && a1 == b1);
}

typedef struct CbDemAuxStruct // 拓扑类 例如实现：查找打散的节点编号，分属哪个单元
{
	int m_nNodes;
	int m_nElems;
	int m_nVarNodes;
	int m_nInitialNodes;
	int m_nInitialElems;

	int m_nContacts; // 接触数量
	int m_nGlobalFaces;
	int m_nGlobalEdges;
	int m_nBrokenFaces;

	double m_fMaxElemLength;
	double m_fMinElemLength;
	double m_fAveElemLength;

	double m_fMaxCoordX;
	double m_fMinCoordX;
	double m_fMaxCoordZ;
	double m_fMinCoordZ;

	CbNodeAux *m_aNodes;
	CbNodePos *m_aVarNodes;
	CbNodePos *m_aDps;
	CbFaceAux *m_aGlobalFaces;
	int *m_aAuxNodeNoFromVarNodeNo;

} CbDemAux;

static void CbDemAux_Init(CbDemAux *pAux)
{
	pAux->m_nNodes = 0;
	pAux->m_nElems = 0;
	pAux->m_nVarNodes = 0;
	pAux->m_nInitialNodes = 0;
	pAux->m_nInitialElems = 0;

	pAux->m_nContacts = 0;
	pAux->m_nBrokenFaces = 0;

	pAux->m_fMaxElemLength = -1.0;
	pAux->m_fMinElemLength = 1e20;
	pAux->m_fAveElemLength = -1.0;

	pAux->m_fMaxCoordX = -1.0;
	pAux->m_fMinCoordX = 1e20;
	pAux->m_fMaxCoordZ = -1.0;
	pAux->m_fMinCoordZ = 1e20;
}

typedef struct CbBulkMaterialStruct // 材料参数
{
	double fDensity;
	double fYoung;
	double fPoisson;
} CbBulkMaterial;

static void CbBulkMaterial_Init(CbBulkMaterial *pBulk)
{
	pBulk->fDensity = 2500; // unit: kg/m3
	pBulk->fYoung = 25e9;	// Pa
	pBulk->fPoisson = 0.25;
}

CbBulkMaterial m_oBulkMaterial;

typedef struct CbSpringStruct
{
	int Brother1NodeNo; // 接触面一侧块体的节点局部序号
	int Brother2NodeNo; // 接触面另一侧块体的节点局部序号
	int Brother1GlobalNodeNo;
	int Brother2GlobalNodeNo;
	int Brother1SpringNodeNo;
	int Brother2SpringNodeNo;

	bool m_bUsed;
	bool m_bBroken; // 0-未破坏 1-拉破坏 2-剪破坏
	bool m_bHistoryBroken;

	double Area;		// 弹簧对应的表面积
	double NormalStiff; // 弹簧法向刚度
	double ShearStiff;	// 弹簧切向刚度

	double m_fSpringLength;
	double rela_dis_local[2]; // 局部坐标系下
	double SpringForce[2];	  // 界面弹簧力
	double Displace[2];
	double m_afEForce[2][2];
} CbSpring;

#define SPRINGS_PER_FORCE_GRP 2
typedef struct CbSpringForceGroupStruct
{
	int m_nSprings;
	CbNodePos m_oNodePos; // bind
	CbForcePos m_aoSpringForces[SPRINGS_PER_FORCE_GRP];

} CbSpringForceGroup;

static void CbSpringForceGroup_Init(CbSpringForceGroup *pSpring)
{
	pSpring->m_nSprings = 0;
	pSpring->m_oNodePos.m_iBlock = -1;
	pSpring->m_oNodePos.m_iNode = -1;
}

typedef struct CbDemStruct
{
	CbDemAux m_pAux;
	CbNodeGroup *m_aSpringNodeGroup;
	CbSpringForceGroup *m_aSpringForceGroup;

	double m_fContactStiffFactor;
	double m_fTouchStiffFactor;

} CbDem;

CbDem BegDem;

static void CbDem_Init(CbDem *pDem)
{
	pDem->m_fContactStiffFactor = 10;
	pDem->m_fTouchStiffFactor = 10;
}

typedef struct CbContactStruct
{
	CbNodePos m_iElem1;	 // One element in contact face
	CbNodePos m_iElem2;	 // Another element in contact face
	int m_iGlobalFaceNo; // The global face number of contact face

	double LocalCoordSystem[2][2]; // Loacl coordinate system in contact face

	CbBulkMaterial m_oBulkMaterial;

	bool m_bUsed;
	bool m_homogeneous;
	int m_iContactModel;		// The contact face constitutive model: 1-linear, 2-brittleMC
	double fContactStiffFactor; // Initial the physical information of contact face

	CbSpring m_aoSpring[2]; // <CbSpring>	The springs in contact face

} CbContact;

CbContact *MyContact; // 全局接触

static void CbContact_Init(CbContact *pContact)
{
	pContact->m_iGlobalFaceNo = -1;
	pContact->m_bUsed = true;
	pContact->m_homogeneous = true;
	pContact->m_iContactModel = 1;
	pContact->fContactStiffFactor = 10;
}

/// Set material of contact face
static int SetMaterial(CbContact *pContact, CbBulkMaterial oMaterial)
{
	pContact->m_oBulkMaterial = oMaterial;
	return 0;
}

/// Obtain material of contatc face
static CbBulkMaterial GetMaterial(CbContact pContact)
{
	return pContact.m_oBulkMaterial;
}

//======================================================================DemRun.c

void LoadParameters();
void InitializeVar();
void GetSize();
void FindIdenticalNodes();
void FindIdenticalFaces();
void GetMesh();
void SetElementArguments(int iModel, CbBulkMaterial m_oBulkMaterial);
void SetContactArguments(int iModel, CbBulkMaterial m_oBulkMaterial);
void Preprocess();
void LocateSourceByCoord(double srcx, double srcz);
void ApplySource(double t, double fpeak, double *source, double sstrength);
void SetVelocityBySel();
void Calculation();
void CalcOneStep();
void FindAdjacentPoints(double pX, double pZ, int radius);
void GetWavefield(double **fux, double **fuz, double **fvx, double **fvz, int it);
void OutputFile(double **UX, double **UZ, double **VX, double **VZ, int nx, int nz, int it);
void UpdateNodeStressStrain();
void UpdateNodePrinStressStrain();
void FreeVar();

//=======================================================================Initialization Function

static void Point_Init(Point *pPoint);
static void CbNode_Init(CbNode *pNode);
static void CbNodePos_Init(CbNodePos *pNode);
static void CbDemAux_Init(CbDemAux *pAux);
static void CbFaceAux_Init(CbFaceAux *pFace);
static void CbBulkMaterial_Init(CbBulkMaterial *pBulk);
static void CbNodePos_Init(CbNodePos *pNode);
static void CbSpringForceGroup_Init(CbSpringForceGroup *pSpring);
static void CbDem_Init(CbDem *pDem);
static void CbContact_Init(CbContact *pContact);

//=====================================================basic function

static bool isEqual(CbFaceAux *left, CbFaceAux *right);
static int SetMaterial(CbContact *pContact, CbBulkMaterial oMaterial);
static CbBulkMaterial GetMaterial(CbContact pContact);

//=====================================================ElementKernel.c

int DoNode_Inc(int threadIdx);
int DoBlock_B_Inc(int threadIdx);

//=====================================================ContactKernel.c

int DoContact_Linear_Inc(int threadIdx);
int SyncSpringGroup(int threadIdx);

//=====================================================Contact.c

int SetSpringNodesGroup();
int SetAllContacts();
int Build(int iContact, double fContactStiffFactor);

#endif
