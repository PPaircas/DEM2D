#include "su.h"
#include "basefun.h"
#include "DemRun.h"
#include "Element.h"
#include <math.h>
#include <stdio.h>
#include <fcntl.h>
#include <unistd.h>
#include <omp.h>

void LoadParameters()
{
	getparstring("meshfnm", &meshfnm);
	getparint("iModel", &iModel);
	getpardouble("fpeak", &fpeak);
	getpardouble("srcx", &srcx);
	getpardouble("srcz", &srcz);
	getpardouble("sstrength", &sstrength);
	getparint("mt", &mt);
	getpardouble("density", &density);
	getpardouble("vp", &vp);
	getpardouble("vs", &vs);

	getparint("nx", &nx);
	getparint("nz", &nz);
	getpardouble("dx", &dx);
	getpardouble("dz", &dz);
	getpardouble("snapt", &snapt);

	// set calculation arguments
	if (!getpardouble("tmax", &tmax))
		err("must specify tmax!\n");
	if (!getpardouble("dt", &dt))
		dt = 1e-3;
	if (!getparint("nt", &nt))
		nt = 1 + tmax / dt;
	if (!getparint("mt", &mt))
		mt = 100;
	if (!getpardouble("xmin", &xmin))
		xmin = 0;
	if (!getpardouble("zmin", &zmin))
		zmin = 0;

	printf("==================PARAMETERS LOADED BEGIN==================\n");
	printf("fpeak=%f\n", fpeak);
	printf("srcx=%f, srcz=%f\n", srcx, srcz);
	printf("tmax=%f\n", tmax);
	printf("dt=%f\n", dt);
	printf("nt=%d\n", nt);
	printf("mt=%d\n", mt);
	printf("snapt=%f\n", snapt);
	printf("density=%f, vp=%f, vs=%f\n", density, vp, vs);
	printf("==================PARAMETERS LOADED END====================\n");
	printf("\n");
}

void GetSize() // 获取网格单元、节点数量
{
	int rowNum = 0;
	char c;
	fp = fopen(meshfnm, "r");
	if (fp == NULL)
	{
		printf("Cannot open file!\n");
	}

	while (!feof(fp)) // 计算文件中数据的行数
	{
		c = fgetc(fp);
		if (c == '\n')
			rowNum++;
	}
	printf("row_number=%d\n", rowNum);

	SetPositionByLine(fp, 1);
	fgets(buf, BUFSIZE, fp);
	sscanf(buf, "%d %d", &(nNodes), &(nElems));
	fclose(fp);
	printf("NodeNum=%d ElemNum=%d\n", nNodes, nElems);

	dNodes = nElems * 3; // 离散单元节点数量

	printf("==================GOT MESH SIZE==================\n");
	printf("\n");
}

void InitializeVar()
{
	// 初始化结构体数组
	CbDem_Init(&BegDem);
	CbDemAux_Init(&BegDem.m_pAux);
	CbBulkMaterial_Init(&m_oBulkMaterial);

	// 计算杨氏模量和泊松比
	double vRatio = vp / vs;
	poisson = (vRatio * vRatio - 2.0) / (2.0 * vRatio * vRatio - 2.0);
	young = 2.0 * density * vs * vs * (1.0 + poisson);
	printf("young=%e, poisson=%f\n", young, poisson);

	pX = alloc1double(nx);
	pZ = alloc1double(nz);
	for (int i = 0; i < nx; i++)
		pX[i] = xmin + dx * i;
	for (int j = 0; j < nz; j++)
		pZ[j] = zmin + dz * j;

	m_fExtForce = 0.0;
	m_fUnbalForce = 0.0;

	// 初始化变量集
	ElemsArea = alloc1double(nElems);
	fU = alloc2double(2, dNodes);
	fV = alloc2double(2, dNodes);
	fA = alloc2double(2, dNodes);
	fUU = alloc2double(2, dNodes);
	fVV = alloc2double(2, dNodes);
	fUX = alloc1double(dNodes);
	fUZ = alloc1double(dNodes);
	fVX = alloc1double(dNodes);
	fVZ = alloc1double(dNodes);
	fNsxx = alloc1double(dNodes);
	fNszz = alloc1double(dNodes);
	fNsxz = alloc1double(dNodes);
	fNexx = alloc1double(dNodes);
	fNezz = alloc1double(dNodes);
	fNexz = alloc1double(dNodes);

	memset((void *)fU[0], 0, (dNodes * 2) * DSIZE);
	memset((void *)fV[0], 0, (dNodes * 2) * DSIZE);
	memset((void *)fA[0], 0, (dNodes * 2) * DSIZE);
	memset((void *)fUU[0], 0, (dNodes * 2) * DSIZE);
	memset((void *)fVV[0], 0, (dNodes * 2) * DSIZE);
	memset(fUX, 0, dNodes * DSIZE);
	memset(fUZ, 0, dNodes * DSIZE);
	memset(fVX, 0, dNodes * DSIZE);
	memset(fVZ, 0, dNodes * DSIZE);
	memset(fNsxx, 0, dNodes * DSIZE);
	memset(fNszz, 0, dNodes * DSIZE);
	memset(fNsxz, 0, dNodes * DSIZE);
	memset(fNexx, 0, dNodes * DSIZE);
	memset(fNezz, 0, dNodes * DSIZE);
	memset(fNexz, 0, dNodes * DSIZE);

	printf("\n");
	printf("==================Initialized Variables==================\n");
	printf("\n");
}

void FindIdenticalNodes() // 找到离散单元节点
{
	int i1, i2;
	int NodeCount = 3;
	int DpCount = 1;
	int iNode;
	int NodesNum = BegDem.m_pAux.m_nNodes;

	BegDem.m_pAux.m_aNodes = (CbNodeAux *)malloc(sizeof(CbNodeAux) * BegDem.m_pAux.m_nNodes); // size: nElems*3
	BegDem.m_pAux.m_aVarNodes = (CbNodePos *)malloc(sizeof(CbNodePos) * BegDem.m_pAux.m_nVarNodes);
	BegDem.m_pAux.m_aAuxNodeNoFromVarNodeNo = (int *)malloc(sizeof(int) * BegDem.m_pAux.m_nVarNodes);
	BegDem.m_pAux.m_aDps = (CbNodePos *)malloc(sizeof(CbNodePos) * m_nDps);

	int *aiNodeElems;
	aiNodeElems = (int *)malloc(sizeof(int) * BegDem.m_pAux.m_nNodes);

	for (int i = 0; i < BegDem.m_pAux.m_nNodes; i++)
	{
		BegDem.m_pAux.m_aNodes[i].m_aIdenticalNodes = (CbNodePos *)malloc(sizeof(CbNodePos) * MAX_SHARE_ELEMENTS);
		CbNodePos_Init(BegDem.m_pAux.m_aNodes[i].m_aIdenticalNodes);
		aiNodeElems[i] = 0;
	}

	for (i1 = 0; i1 < nElems; i1++)
	{
		for (i2 = 0; i2 < NodeCount; i2++)
		{
			iNode = i1 * NodeCount + i2;
			int eleId = Elem[i1].VertexID[i2];	  // 离散前的单元节点全局编号
			int varId = MyTri3[i1].m_agNodes[i2]; // 离散后的单元节点全局编号

			CbNodePos oNodePos;
			oNodePos.m_iBlock = i1;
			oNodePos.m_iNode = i2;

			BegDem.m_pAux.m_aNodes[iNode].m_aiGlobalNodeNo = eleId;
			BegDem.m_pAux.m_aVarNodes[varId] = oNodePos;			// 通过varId索引离散节点从属单元
			BegDem.m_pAux.m_aAuxNodeNoFromVarNodeNo[varId] = eleId; // 通过varId索引eleId
		}

		for (i2 = 0; i2 < DpCount; i2++)
		{
			CbNodePos oNodePos;
			oNodePos.m_iBlock = i1;
			oNodePos.m_iNode = i2;

			int dpId = MyTri3[i1].m_iDpPos + i2;
			BegDem.m_pAux.m_aDps[dpId] = oNodePos; // 找到高斯积分点从属单元
		}
	}

	for (i1 = 0; i1 < nElems; i1++)
	{
		for (i2 = 0; i2 < NodeCount; i2++)
		{
			iNode = i1 * NodeCount + i2;

			CbNodePos oNodePos;
			oNodePos.m_iBlock = i1;
			oNodePos.m_iNode = i2;

			int eleId = Elem[i1].VertexID[i2]; // 离散前的单元节点全局编号
			BegDem.m_pAux.m_aNodes[iNode].m_aiGlobalNodeNo = eleId;

			// 找到每个离散节点所从属的单元号及局部编号
			for (iNode = 0; iNode < NodesNum; iNode++)
			{
				if (BegDem.m_pAux.m_aNodes[iNode].m_aiGlobalNodeNo == eleId)
				{
					BegDem.m_pAux.m_aNodes[iNode].m_aIdenticalNodes[aiNodeElems[iNode]] = oNodePos;
					aiNodeElems[iNode]++;
					break;
				}
			}
		}
	}

	printf("iBlock1=%d\n", BegDem.m_pAux.m_aNodes[0].m_aIdenticalNodes[0].m_iBlock);
	printf("iNode1=%d\n", BegDem.m_pAux.m_aNodes[0].m_aIdenticalNodes[0].m_iNode);

	printf("iBlock2=%d\n", BegDem.m_pAux.m_aNodes[0].m_aIdenticalNodes[1].m_iBlock);
	printf("iNode2=%d\n", BegDem.m_pAux.m_aNodes[0].m_aIdenticalNodes[1].m_iNode);

	printf("iNode=%d\n", iNode);
	printf("nNodes=%d\n", NodesNum);

	printf("\n");
	printf("==================Find Identical Nodes==================\n");
	printf("\n");
}

void FindIdenticalFaces() // 找到全局边数，边上全局节点号，共用该边的单元号，边在单元上的局部编号
{
	int i1, i2, i3, i4;
	int *aiNodeFaces; // 所有离散后节点, 具有相同全局编号所占边数

	int nFaces = BegDem.m_pAux.m_nNodes; // size: nElems*3
	int iFace;							 // face serial no
	aiNodeFaces = (int *)malloc(sizeof(int) * nFaces);
	BegDem.m_pAux.m_aGlobalFaces = (CbFaceAux *)malloc(sizeof(CbFaceAux) * nFaces);

	for (int i = 0; i < nFaces; i++)
	{
		CbFaceAux_Init(&(BegDem.m_pAux.m_aGlobalFaces[i]));
		BegDem.m_pAux.m_aGlobalFaces[i].m_aIdenticalFaces = (CbNodePos *)malloc(sizeof(CbNodePos *) * 2);
		CbNodePos_Init(&(BegDem.m_pAux.m_aGlobalFaces[i].m_aIdenticalFaces[0]));
		CbNodePos_Init(&(BegDem.m_pAux.m_aGlobalFaces[i].m_aIdenticalFaces[1]));
		aiNodeFaces[i] = 0;
	}

	for (i1 = 0; i1 < nElems; i1++) // 遍历每个三角形单元
	{
		MyElem pElem = MyTri3[i1];
		for (i2 = 0; i2 < 3; i2++) // 遍历三角单元的每条边
		{
			iFace = i1 * 3 + i2;
			for (i3 = 0; i3 < 2; i3++) // 遍历边上的两个顶点
			{
				int iElemNode = pElem.m_bgNodes[pElem.m_aiEdgeNodes[i2][i3]];
				BegDem.m_pAux.m_aGlobalFaces[iFace].m_aiGlobalNodeNo[i3] = iElemNode;
			} // m_aiGlobalNodeNo存储了所有边的所有顶点

			CbFaceAux bFace = BegDem.m_pAux.m_aGlobalFaces[iFace];

			if (bFace.m_aiGlobalNodeNo[0] > bFace.m_aiGlobalNodeNo[1]) // 保持节点编号从小到大存储
			{
				int temp = bFace.m_aiGlobalNodeNo[0];
				BegDem.m_pAux.m_aGlobalFaces[iFace].m_aiGlobalNodeNo[0] = bFace.m_aiGlobalNodeNo[1];
				BegDem.m_pAux.m_aGlobalFaces[iFace].m_aiGlobalNodeNo[1] = temp;
			}
		}
	}

	for (i1 = 0; i1 < nElems; i1++)
	{
		MyElem pElem = MyTri3[i1];
		for (i2 = 0; i2 < 3; i2++) // 遍历三角形单元的三条边
		{
			CbFaceAux oFace;
			CbFaceAux_Init(&oFace);
			oFace.m_aIdenticalFaces = (CbNodePos *)malloc(sizeof(CbNodePos *) * 2);
			CbNodePos_Init(&oFace.m_aIdenticalFaces[0]);
			CbNodePos_Init(&oFace.m_aIdenticalFaces[1]);

			for (i3 = 0; i3 < 2; i3++) // 遍历每条边上的单元节点
			{
				int iElemNode = pElem.m_bgNodes[pElem.m_aiEdgeNodes[i2][i3]]; // 离散前边上的节点全局编号
				oFace.m_aiGlobalNodeNo[i3] = iElemNode;
			}

			if (oFace.m_aiGlobalNodeNo[0] > oFace.m_aiGlobalNodeNo[1]) // 保持节点编号从小到大排列
			{
				int temp = oFace.m_aiGlobalNodeNo[0];
				oFace.m_aiGlobalNodeNo[0] = oFace.m_aiGlobalNodeNo[1];
				oFace.m_aiGlobalNodeNo[1] = temp;
			}

			CbNodePos oNodePos;
			oNodePos.m_iBlock = i1;
			oNodePos.m_iNode = i2; // 边的编号

			for (iFace = 0; iFace < nFaces; iFace++)
			{
				if (isEqual(&BegDem.m_pAux.m_aGlobalFaces[iFace], &oFace))
				{
					// printf("iFace=%d\n",iFace);
					BegDem.m_pAux.m_aGlobalFaces[iFace].m_aIdenticalFaces[aiNodeFaces[iFace]] = oNodePos;
					aiNodeFaces[iFace]++;
					break;
				}
			}
		}
	}

	printf("blockid1=%d\n", BegDem.m_pAux.m_aGlobalFaces[0].m_aIdenticalFaces[0].m_iBlock);
	printf("edgeid1=%d\n", BegDem.m_pAux.m_aGlobalFaces[0].m_aIdenticalFaces[0].m_iNode);
	printf("blockid2=%d\n", BegDem.m_pAux.m_aGlobalFaces[0].m_aIdenticalFaces[1].m_iBlock);
	printf("edgeid2=%d\n", BegDem.m_pAux.m_aGlobalFaces[0].m_aIdenticalFaces[1].m_iNode);

	int nContacts = 0;
	for (iFace = 0; iFace < nFaces; iFace++) // 找到单元之间接触的边
	{
		if (aiNodeFaces[iFace] == 2)
		{
			nContacts++;
			BegDem.m_pAux.m_aGlobalFaces[iFace].isContact = true;
		}
	}

	// int nContacts = nContactEdges / 2; 		// 所有接触的数量

	BegDem.m_pAux.m_nGlobalFaces = nFaces; // 所有的离散单元的边数
	BegDem.m_pAux.m_nContacts = nContacts;
	printf("iFace=%d\n", iFace);
	printf("nFaces=%d\n", nFaces);
	// printf("nContactEdges=%d\n", nContactEdges);
	printf("nContacts=%d\n", nContacts);

	printf("\n");
	printf("==================Find Identical Edges==================\n");
	printf("\n");
}

void GetMesh()
{
	Node = (NodeInfo *)malloc(sizeof(NodeInfo) * nNodes);
	Elem = (ElemInfo *)malloc(sizeof(ElemInfo) * nElems);
	MyTri3 = (MyElem *)malloc(sizeof(MyElem) * nElems);

	// 从文件里将网格信息读入到指定数据结构中存储
	int i, j;
	int beglineNum = 2;

	fp = fopen(meshfnm, "r");
	if (fp == NULL)
	{
		printf("Cannot open file!\n");
	}

	SetPositionByLine(fp, beglineNum); // 将指针移向指定行

	for (i = 0; i < nNodes; i++)
	{
		fgets(buf, BUFSIZE, fp);															// 一次读取一行
		sscanf(buf, "%d %f %f %*f", &(Node[i].NodeId), &(Node[i].NodeX), &(Node[i].NodeZ)); // 分别读取序号、坐标数据
	}

	for (j = 0; j < nElems; j++)
	{
		fgets(buf, BUFSIZE, fp); // 一次读取一行
		sscanf(buf, "%d %*d %*d %d %d %d", &(Elem[j].ElemId), &(Elem[j].VertexId1), &(Elem[j].VertexId2), &(Elem[j].VertexId3));

		Elem[j].VertexID[0] = Elem[j].VertexId1;
		Elem[j].VertexID[1] = Elem[j].VertexId2;
		Elem[j].VertexID[2] = Elem[j].VertexId3;
	}

	fclose(fp);

	printf("==================MESH READ BEGIN==================\n");
	printf("Node info read:\n");
	for (i = 0; i < 3; i++)
	{
		printf("NodeId=%d NodeX=%f NodeZ=%f\n", Node[i].NodeId, Node[i].NodeX, Node[i].NodeZ);
	}
	printf("\nElement info read:\n");
	for (j = 0; j < 3; j++)
	{
		printf("ElemId=%d VertexId1=%d VertexId2=%d VertexId3=%d\n", Elem[j].ElemId, Elem[j].VertexId1, Elem[j].VertexId2, Elem[j].VertexId3);
	}

	printf("==================MESH READ END====================\n");

	//========================================================================

	GetElemsCoord(); // 获取所有单元顶点坐标信息
	GetTriNodes();	 //	获取单元节点的局部序号、全局序号

	BuildDps(); // 创建高斯积分点

	dfX = alloc1double(dNodes); // 存储离散后的单元节点坐标
	dfZ = alloc1double(dNodes);

	// 把节点打散，存储在不同的单元
	int nVarNodes = 0; // 存放离散后节点数
	for (int i1 = 0; i1 < nElems; i1++)
	{
		for (int inode = 0; inode < 3; inode++)
		{
			int iGlobalNodeNo = MyTri3[i1].m_bgNodes[inode]; // 从网格中读取的节点全局序号

			dfX[nVarNodes] = afX[iGlobalNodeNo - 1];
			dfZ[nVarNodes] = afZ[iGlobalNodeNo - 1];

			MyTri3[i1].m_agNodes[inode] = nVarNodes; // 对节点按照单元编号的次序重新编号（全局序号）
			nVarNodes++;
		}
	}

	printf("nVarNodes=%d\n", nVarNodes);

	//////////////////////////////////////////////////// 导入网格结束
	BegDem.m_pAux.m_nElems = nElems;
	BegDem.m_pAux.m_nInitialElems = nElems;
	BegDem.m_pAux.m_nNodes = nVarNodes;		   // size: nElems*3
	BegDem.m_pAux.m_nInitialNodes = nVarNodes; // size: nElems*3
	BegDem.m_pAux.m_nVarNodes = nVarNodes;	   // size: nVarNodes = nElems*3

	FindIdenticalNodes(); // 找到单元离散后的节点
	FindIdenticalFaces(); // 找到单元之间接触的边

	SetAllContacts(); // 设置接触，关键一步，让单元离散添加弹簧

	GetMax_MinElemLength(); // 计算三角单元边长的最大最小值
	GetMax_MinNodeCoord();	// 计算单元节点坐标的最大最小值
}

void SetElementArguments(int iModel, CbBulkMaterial m_oBulkMaterial)
{
	for (int iElem = 0; iElem < nElems; iElem++)
	{
		MyTri3[iElem].m_iElemModel = iModel;		 // 对单元施加模型
		MyTri3[iElem].m_oMaterial = m_oBulkMaterial; // 设置单元材料参数
	}
	printf("\n");
	printf("==================Set Elem Arguments====================\n");
}

void SetContactArguments(int iModel, CbBulkMaterial m_oBulkMaterial)
{
	int nContacts = BegDem.m_pAux.m_nContacts;
	for (int i1 = 0; i1 < nContacts; i1++)
	{
		MyContact[i1].m_iContactModel = iModel; // 对接触面施加模型

		if (MyContact->m_homogeneous)
		{
			MyContact[i1].m_oBulkMaterial = m_oBulkMaterial; // 设置接触面材料参数
		}
		else
		{
			// 如果是非均匀介质，找到接触面两侧的单元材料参数，取最小的值赋给接触面
			CbNodePos iNodePos1 = MyContact[i1].m_iElem1;
			CbNodePos iNodePos2 = MyContact[i1].m_iElem2;

			CbBulkMaterial oMat1 = GetMaterial(MyContact[iNodePos1.m_iBlock]);
			CbBulkMaterial oMat2 = GetMaterial(MyContact[iNodePos2.m_iBlock]);

			CbBulkMaterial oMat;
			oMat.fDensity = fmin(oMat1.fDensity, oMat2.fDensity);
			oMat.fYoung = fmin(oMat1.fYoung, oMat2.fYoung);
			oMat.fPoisson = fmin(oMat1.fPoisson, oMat2.fPoisson);

			SetMaterial(&MyContact[i1], oMat);
		}
	}
	printf("\n");
	printf("==================Set Contact Arguments====================\n");
}

void Preprocess()
{
	for (int iElem = 0; iElem < nElems; iElem++)
	{
		GenerateMass(iElem);
		GenerateStiff(iElem);
	}
	printf("Mass0=%e\n", MyTri3[0].aNodes[0].m_fMass);
	printf("stiff0=%e\n", MyTri3[0].Stiff[0][0]);

	double fContactStiffFactor = BegDem.m_fContactStiffFactor;

	for (int iContact = 0; iContact < BegDem.m_pAux.m_nContacts; iContact++)
	{
		Build(iContact, fContactStiffFactor); // 计算接触面弹簧的法向刚度、切向刚度
	}
	printf("NormalStiff=%e\n", MyContact[0].m_aoSpring[0].NormalStiff);
	printf("ShearStiff=%e\n", MyContact[0].m_aoSpring[0].ShearStiff);

	// 找到离震源坐标最近的节点，将震源施加在该节点上
	LocateSourceByCoord(srcx, srcz);
	printf("sGlobalNo=%d, sElemId=%d, sNodeId=%d\n", MyTri3[sElemId].m_agNodes[sNodeId], sElemId, sNodeId);
	printf("srcX=%f, srcZ=%f\n", srcx, srcz);
	printf("sNodeCoordX=%f, sNodeCoordZ=%f\n", dfX[MyTri3[sElemId].m_agNodes[sNodeId]], dfZ[MyTri3[sElemId].m_agNodes[sNodeId]]);

	printf("\n");
	printf("==================Preprocess====================\n");
}

void LocateSourceByCoord(double srcx, double srcz)
{
	// 给定震源坐标，寻找距离震源坐标最近的点，返回该节点所在单元编号，及其再单元内局部编号
	int iElem, iNode;
	double GivenCoord[2] = {srcx, srcz}, QueryNodeCoord[2];
	double QueryDistance = 1.0e8, TempDistance;
	double AveElemLength = BegDem.m_pAux.m_fAveElemLength;

	for (iElem = 0; iElem < nElems; iElem++)
	{
		MyElem pMyElem = MyTri3[iElem];
		for (iNode = 0; iNode < 3; iNode++)
		{
			QueryNodeCoord[0] = afX[pMyElem.m_bgNodes[iNode] - 1];
			QueryNodeCoord[1] = afZ[pMyElem.m_bgNodes[iNode] - 1];

			TempDistance = TwoPointDistance(GivenCoord, QueryNodeCoord);
			if (TempDistance < QueryDistance)
			{
				QueryDistance = TempDistance;
				sElemId = iElem;
				sNodeId = iNode;
			}
		}
	}

	printf("querydistance=%f\n", QueryDistance);
	printf("\n");
	printf("====================Locate Source Position========================\n");
}

void ApplySource(double t, double fpeak, double *source, double sstrength) // 施加源项
{
	// 震源函数采用给定点上延迟的雷克子波作震源
	double tdelay, ts;
	int it = t / dt;

	tdelay = 1.0 / fpeak;
	if (t > 2.0 * tdelay)
		ts = 0;
	ts = ricker(t - tdelay, fpeak);
	source[it] = sstrength * density * ts;

	//MyTri3[sElemId].aNodes[sNodeId].m_afExtForce[0] = 0.5 * source[it]; // x方向施加震源
	MyTri3[sElemId].aNodes[sNodeId].m_afExtForce[1] = source[it]; // z方向施加震源
																  // printf("t=%f, source=%f\n", t, source);
}

void Calculation()
{
	int snapit = snapt / dt;
	double **aUx, **aUz;
	double **aVx, **aVz;
	aUx = alloc2double(dNodes, nt);
	aUz = alloc2double(dNodes, nt);
	aVx = alloc2double(dNodes, nt);
	aVz = alloc2double(dNodes, nt);
	source = alloc1double(nt);

	for (it = 0, t = 0.0; it < nt; it++, t += dt)
	{
		ApplySource(t, fpeak, source, sstrength); // 加载震源

		CalcOneStep(); // do one time step

		for (int iNode = 0; iNode < dNodes; iNode++)
		{
			aUx[it][iNode] = fUX[iNode];
			aUz[it][iNode] = fUZ[iNode];
			aVx[it][iNode] = fVX[iNode];
			aVz[it][iNode] = fVZ[iNode];
		}

		// if (it % mt == 0 && it <= 100)
		// {
		// 	printf("it=%d, time=%f\n", it, t);
		// 	GetWavefield(aUx, aUz, aVx, aVz, it);
		// }

		if (it == snapit)
		{
			printf("snapit=%d, snapt=%f s\n", it, snapit * dt);
			GetWavefield(aUx, aUz, aVx, aVz, snapit);
		}

		fU = fUU;
		fV = fVV;
	}

	// 输出震源
	FILE *file = fopen("source.txt", "w");
	for (it = 0, t = 0.0; it < nt; ++it, t += dt)
	{
		fprintf(file, "%f %f\n ", t, source[it]);
	}
	fclose(file);

	free2double(aUx);
	free2double(aUz);
	free2double(aVx);
	free2double(aVz);
	printf("====================Calculation Wavefield===================\n");
}

void CalcOneStep()
{
	/// Calculation the block elements.
#pragma omp parallel for num_threads(40)
	for (int iElem = 0; iElem < nElems; iElem++)
	{
		DoBlock_B_Inc(iElem);
	}

	// UpdateNodeStressStrain();
	// UpdateNodePrinStressStrain();

	/// Calculation the contact faces.
	int nContacts = BegDem.m_pAux.m_nContacts;
#pragma omp parallel for num_threads(40)
	for (int iContact = 0; iContact < nContacts; iContact++)
	{
		DoContact_Linear_Inc(iContact);
	}

	/// Integrating nodal forces.
	int nSpriGroup = BegDem.m_pAux.m_nContacts * 2;
#pragma omp parallel for num_threads(40)
	for (int iSpriGroup = 0; iSpriGroup < nSpriGroup; iSpriGroup++)
	{
		SyncSpringGroup(iSpriGroup);
	}
}

void FindAdjacentPoints(double pX, double pZ, int radius)
{
	int MaxPointNum = radius * radius * 2;
	numAdjacentPoints = 0;
	memset(weights, 0, MaxPointNum * DSIZE);
	for (int i = 0; i < MaxPointNum; i++)
		Point_Init(&adjacentPoints[i]);

#pragma omp parallel for num_threads(40)
	for (int i = 0; i < dNodes; i++)
	{
		// 判断点是否在圆内
		if (pow(dfX[i] - pX, 2) + pow(dfZ[i] - pZ, 2) <= pow(radius, 2))
		{
			// 计算距离加权值
			double weight = (radius - sqrt(pow(dfX[i] - pX, 2) + pow(dfZ[i] - pZ, 2))) / radius;
			Point adjacentPoint = {dfX[i], dfZ[i], i};
#pragma omp critical
			{
				adjacentPoints[numAdjacentPoints] = adjacentPoint;
				weights[numAdjacentPoints] = weight;
				numAdjacentPoints++;
			}
		}
	}
}

void GetWavefield(double **fux, double **fuz, double **fvx, double **fvz, int it)
{
	double **UX, **UZ; // 全局位移场
	double **VX, **VZ; // 全局速度场
	UX = alloc2double(nz, nx);
	UZ = alloc2double(nz, nx);
	VX = alloc2double(nz, nx);
	VZ = alloc2double(nz, nx);

	float factor = 1.0;
	int radius = round(BegDem.m_pAux.m_fAveElemLength * factor);
	int MaxPointNum = radius * radius * 2;
	printf("radius=%d, MaxPointNum=%d\n", radius, MaxPointNum);

	adjacentPoints = (Point *)malloc(sizeof(Point) * MaxPointNum);
	weights = alloc1double(MaxPointNum);

	for (int i = 0; i < nx; i++)
	{
		for (int j = 0; j < nz; j++)
		{
			FindAdjacentPoints(pX[i], pZ[j], radius);

			for (int k = 0; k < numAdjacentPoints; k++)
			{
				UX[i][j] += weights[k] * fux[it][adjacentPoints[k].GlobalNodeNo];
				UZ[i][j] += weights[k] * fuz[it][adjacentPoints[k].GlobalNodeNo];

				VX[i][j] += weights[k] * fvx[it][adjacentPoints[k].GlobalNodeNo];
				VZ[i][j] += weights[k] * fvz[it][adjacentPoints[k].GlobalNodeNo];
			}
		}
	}

	OutputFile(UX, UZ, VX, VZ, nx, nz, it);

	free2double(UX);
	free2double(UZ);
	free2double(VX);
	free2double(VZ);
}

void OutputFile(double **UX, double **UZ, double **VX, double **VZ, int nx, int nz, int it)
{
	FILE *file;

	char filename1[] = "UX_";
	char filename2[] = "UZ_";
	char filename3[] = "VX_";
	char filename4[] = "VZ_";
	char filetype[] = ".txt";
	char extension[20];

	sprintf(extension, "%d%s", it, filetype);
	printf("extension=%s\n", extension);

	char UxName[20], UzName[20], VxName[20], VzName[20];
	sprintf(UxName, "%s%s", filename1, extension);
	sprintf(UzName, "%s%s", filename2, extension);
	sprintf(VxName, "%s%s", filename3, extension);
	sprintf(VzName, "%s%s", filename4, extension);

	// file = fopen(UxName, "w"); // 打开文件，以写入模式打开（清空文件内容，如果文件已存在）
	// for (int i = 0; i < nx; i++)
	// {
	// 	for (int j = 0; j < nz; j++)
	// 	{
	// 		fprintf(file, "%e ", UX[i][j]);
	// 	}
	// 	fprintf(file, "\n");
	// }
	// fclose(file);

	// file = fopen(UzName, "w");
	// for (int i = 0; i < nx; i++)
	// {
	// 	for (int j = 0; j < nz; j++)
	// 	{
	// 		fprintf(file, "%e ", UZ[i][j]);
	// 	}
	// 	fprintf(file, "\n");
	// }
	// fclose(file);

	// file = fopen("UX.bin", "wb");
	// fwrite(UX[0], DSIZE, nx * nz, file);
	// fclose(file);

	// file = fopen("UZ.bin", "wb");
	// fwrite(UZ[0], DSIZE, nx * nz, file);
	// fclose(file);

	file = fopen(VxName, "w"); // 打开文件，以写入模式打开（清空文件内容，如果文件已存在）
	for (int i = 0; i < nx; i++)
	{
		for (int j = 0; j < nz; j++)
		{
			fprintf(file, "%e ", VX[i][j]);
		}
		fprintf(file, "\n");
	}
	fclose(file);

	file = fopen(VzName, "w");
	for (int i = 0; i < nx; i++)
	{
		for (int j = 0; j < nz; j++)
		{
			fprintf(file, "%e ", VZ[i][j]);
		}
		fprintf(file, "\n");
	}
	fclose(file);

	// file = fopen("VX.bin", "wb");
	// fwrite(VX[0], FSIZE, nx * nz, file);
	// fclose(file);

	// file = fopen("VZ.bin", "wb");
	// fwrite(VZ[0], FSIZE, nx * nz, file);
	// fclose(file);
}

void SetVelocityBySel() // 施加速度边界条件，默认施加在节点上
{
}

void UpdateNodeStressStrain()
{
	// 首先将高斯点应力作为节点力赋值给各个单元节点，然后将同一位置的未破坏节点求平均后再赋值给各个节点
	double **afS, **afE, **afNS, **afNE;

	afS = alloc2double(3, m_nDps);
	afE = alloc2double(3, m_nDps);
	afNS = alloc2double(3, dNodes);
	afNE = alloc2double(3, dNodes);

	for (int i2 = 0; i2 < m_nDps; i2++)
	{
		afS[i2][0] = fsxx[i2];
		afS[i2][1] = fszz[i2];
		afS[i2][2] = fsxz[i2];

		afE[i2][0] = fexx[i2];
		afE[i2][1] = fezz[i2];
		afE[i2][2] = fexz[i2];
	}

	for (int ii = 0; ii < dNodes; ii++)
	{
		afNS[ii][0] = fNsxx[ii];
		afNS[ii][1] = fNszz[ii];
		afNS[ii][2] = fNsxz[ii];

		afNE[ii][0] = fNexx[ii];
		afNE[ii][1] = fNezz[ii];
		afNE[ii][2] = fNexz[ii];
	}

	for (int ielem = 0; ielem < nElems; ielem++)
	{
		MyElem pMyElem = MyTri3[ielem];
		int nDps = 1;
		int iNodes = 3;
		// 将高斯点应力作为节点应力赋值给各个单元节点
		if (nDps == 1)
		{
			for (int inode = 0; inode < iNodes; inode++)
			{
				int iVarNodeNo = pMyElem.m_agNodes[inode];
				for (int i = 0; i < 3; i++)
				{
					afNS[iVarNodeNo][i] = afS[pMyElem.m_iDpPos][i];
					afNE[iVarNodeNo][i] = afE[pMyElem.m_iDpPos][i];
				}
			}
		}
		else
		{
			for (int inode = 0; inode < iNodes; inode++)
			{
				int iVarNodeNo = pMyElem.m_aiNodes[inode];

				for (int i = 0; i < 3; i++)
				{
					afNS[iVarNodeNo][i] = afS[pMyElem.m_iDpPos + inode][i];
					afNE[iVarNodeNo][i] = afE[pMyElem.m_iDpPos + inode][i];
				}
			}
		}
	}

	int nnNodes = BegDem.m_pAux.m_nNodes;
	for (int inode = 0; inode < nnNodes; inode++)
	{
		printf("inode=%d, nnNodes=%d\n", inode, nnNodes);
		int nNodeCount = 0;
		double afAverageStress[3] = {0.0};
		double afAverageStrain[3] = {0.0};

		CbNodePos *aoNodePos;
		aoNodePos = (CbNodePos *)malloc(sizeof(CbNodePos) * MAX_SHARE_ELEMENTS);
		for (int i = 0; i < MAX_SHARE_ELEMENTS; i++)
		{
			CbNodePos_Init(&aoNodePos[i]);
		}
		aoNodePos = BegDem.m_pAux.m_aNodes[inode].m_aIdenticalNodes;

		for (int i1 = 0; i1 < MAX_SHARE_ELEMENTS; i1++)
		{
			int ib = aoNodePos[i1].m_iBlock;
			int in = aoNodePos[i1].m_iNode;

			if (ib >= 0 && in >= 0)
			{
				if (MyTri3[ib].aNodes[in].m_iBrokenType == 0 && MyTri3[ib].m_iElemModel > 0)
				{
					nNodeCount++;
				}
			}
		}

		for (int i1 = 0; i1 < MAX_SHARE_ELEMENTS; i1++)
		{
			if (nNodeCount > 0)
			{
				int ib = aoNodePos[i1].m_iBlock;
				int in = aoNodePos[i1].m_iNode;

				if (ib >= 0 && in >= 0)
				{
					int iVarNodeNO = MyTri3[ib].m_aiNodes[in];

					// 将同一位置的未破坏节点应力求平均后再赋值给各个节点
					if (MyTri3[ib].aNodes[in].m_iBrokenType == 0 && MyTri3[ib].m_iElemModel > 0)
					{
						for (int i2 = 0; i2 < 3; i2++)
						{
							afAverageStress[i2] += afNS[iVarNodeNO][i2] / nNodeCount;
							afAverageStrain[i2] += afNE[iVarNodeNO][i2] / nNodeCount;
						}

						fNsxx[iVarNodeNO] = afAverageStress[0];
						fNszz[iVarNodeNO] = afAverageStress[1];
						fNsxz[iVarNodeNO] = afAverageStress[2];
						fNexx[iVarNodeNO] = afAverageStrain[0];
						fNezz[iVarNodeNO] = afAverageStrain[1];
						fNexz[iVarNodeNO] = afAverageStrain[2];
					}
				}
			}
		}
	}

	free2double(afS);
	free2double(afE);
	free2double(afNS);
	free2double(afNE);
	printf("====================Update Node Stress Strain=====================\n");
}

void UpdateNodePrinStressStrain()
{
	double **afSN, **afEN;
	afSN = alloc2double(3, dNodes);
	afEN = alloc2double(3, dNodes);

	for (int ii = 0; ii < dNodes; ii++)
	{
		afSN[ii][0] = fNsxx[ii];
		afSN[ii][1] = fNszz[ii];
		afSN[ii][2] = fNsxz[ii];

		afEN[ii][0] = fNexx[ii];
		afEN[ii][1] = fNezz[ii];
		afEN[ii][2] = fNexz[ii];
	}

	for (int inode = 0; inode < BegDem.m_pAux.m_nNodes; inode++)
	{
		printf("inode=%d, nnNodes=%d\n", inode, BegDem.m_pAux.m_nNodes);
		CbNodePos *aoNodePos;
		aoNodePos = (CbNodePos *)malloc(sizeof(CbNodePos) * MAX_SHARE_ELEMENTS);
		for (int i = 0; i < MAX_SHARE_ELEMENTS; i++)
		{
			CbNodePos_Init(&aoNodePos[i]);
		}
		aoNodePos = BegDem.m_pAux.m_aNodes[inode].m_aIdenticalNodes;

		for (int i1 = 0; i1 < MAX_SHARE_ELEMENTS; i1++)
		{
			int ib = aoNodePos[i1].m_iBlock;
			int in = aoNodePos[i1].m_iNode;

			if (ib >= 0 && in >= 0)
			{

				int iVarNodeNO = MyTri3[ib].m_aiNodes[in];

				double afStress[3] = {0.0};
				double afStrain[3] = {0.0};

				for (int i2 = 0; i2 < 3; i2++)
				{
					afStress[i2] = afSN[iVarNodeNO][i2];
					afStrain[i2] = afEN[iVarNodeNO][i2];
				}

				double PrinStress[2] = {0.0};
				double PrinStrain[2] = {0.0};
				double lmStress[2][2] = {0.0}; // 方向余弦
				double lmStrain[2][2] = {0.0};

				Cal_Principal_Stress(afStress, PrinStress, lmStress);
				Cal_Principal_Stress(afStrain, PrinStrain, lmStrain);

				for (int i = 0; i < 2; i++)
				{
					MyTri3[ib].aNodes[in].m_afPrinStress[i] = PrinStress[i];
					MyTri3[ib].aNodes[in].m_afPrinStrain[i] = PrinStrain[i];

					for (int j = 0; j < 2; j++)
					{
						MyTri3[ib].aNodes[in].m_afPrinStressDir[i][j] = lmStress[i][j];
						MyTri3[ib].aNodes[in].m_afPrinStrainDir[i][j] = lmStrain[i][j];
					}
				}
			}
		}
	}

	free2double(afSN);
	free2double(afEN);
	printf("====================Update Node Prin Stress Strain=====================\n");
}

void FreeVar()
{
	free(Node);
	free(Elem);
	free(MyTri3);
	free(MyContact);
	free1double(ElemsArea);
	free1double(source);

	free1double(pX);
	free1double(pZ);
	free(adjacentPoints);
	free1double(weights);

	free1double(afX);
	free1double(afZ);
	free1double(X1);
	free1double(X2);
	free1double(X3);
	free1double(Z1);
	free1double(Z2);
	free1double(Z3);
	free1double(dfX);
	free1double(dfZ);

	free2double(fU);
	free2double(fV);
	free2double(fA);
	// free2double(fUU);
	// free2double(fVV);

	free1double(fUX);
	free1double(fUZ);
	free1double(fVX);
	free1double(fVZ);
	free1double(fsxx);
	free1double(fszz);
	free1double(fsxz);
	free1double(fexx);
	free1double(fezz);
	free1double(fexz);

	free1double(fNsxx);
	free1double(fNszz);
	free1double(fNsxz);
	free1double(fNexx);
	free1double(fNezz);
	free1double(fNexz);

	printf("==================Variables Memory Have Freed===============\n");
}
