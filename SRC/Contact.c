#include <math.h>
#include "basefun.h"
#include "DemRun.h"
#include "Element.h"

int SetSpringNodesGroup() // 把节点编号信息赋给m_aSpringNodeGroup，构成弹簧节点编号信息
{
    int i1, i2;
    int NodeSize = BegDem.m_pAux.m_nNodes; // size: nElems * 3

    BegDem.m_aSpringNodeGroup = (CbNodeGroup *)calloc(NodeSize, sizeof(CbNodeGroup));

    for (i1 = 0; i1 < NodeSize; i1++)
    {
        CbNodePos *aoNodePos = BegDem.m_pAux.m_aNodes[i1].m_aIdenticalNodes;
        BegDem.m_aSpringNodeGroup[i1].m_aNodePos = (CbNodePos *)malloc(sizeof(CbNodePos) * MAX_SHARE_ELEMENTS);

        for (i2 = 0; i2 < MAX_SHARE_ELEMENTS; i2++)
        {
            CbNodePos oNodePos = aoNodePos[i2];
            BegDem.m_aSpringNodeGroup[i1].m_aNodePos[i2] = oNodePos;
        }
    }

    return 0;
}

int SetAllContacts()
{
    SetSpringNodesGroup(); // 把节点编号信息赋给m_aSpringNodeGroup，构成弹簧编号信息

    int nContacts = BegDem.m_pAux.m_nContacts;
    MyContact = (CbContact *)malloc(sizeof(CbContact) * nContacts);

    for (int i = 0; i < nContacts; i++)
        CbContact_Init(&MyContact[i]);

    CbContact *pContact; // 抽象接触面类：相互接触的两块体的编号，相互接触的两块体的材料

    int nSprings = nContacts * 2;
    int nSpringNodes = nContacts * 4;

    int iSpringNode = 0;
    int iContact = 0;

    for (int iGlobalFace = 0; iGlobalFace < BegDem.m_pAux.m_nGlobalFaces; iGlobalFace++)
    {
        CbFaceAux *piFace = &BegDem.m_pAux.m_aGlobalFaces[iGlobalFace];
        if (piFace->isContact)
        {
            if (iContact < nContacts)
            {
                pContact = &MyContact[iContact];
                iContact++;
            }

            // printf("iGlobalFace=%d\n", iGlobalFace);
            pContact->m_iGlobalFaceNo = iGlobalFace; // 所有接触面在所有面里的编号（全局编号）
            pContact->m_iElem1.m_iBlock = piFace->m_aIdenticalFaces[0].m_iBlock;
            pContact->m_iElem1.m_iNode = piFace->m_aIdenticalFaces[0].m_iNode; // 边
            pContact->m_iElem2.m_iBlock = piFace->m_aIdenticalFaces[1].m_iBlock;
            pContact->m_iElem2.m_iNode = piFace->m_aIdenticalFaces[1].m_iNode;
            pContact->m_bUsed = true;

            int *aiNodes1 = Elem[pContact->m_iElem1.m_iBlock].VertexID; // 节点全局编号
            int *aiNodes2 = Elem[pContact->m_iElem2.m_iBlock].VertexID;

            int *aiVarNodes1 = MyTri3[pContact->m_iElem1.m_iBlock].m_aiNodes; // 节点单元内局部编号
            int *aiVarNodes2 = MyTri3[pContact->m_iElem2.m_iBlock].m_aiNodes;

            int *aiEdgeNodes1 = MyTri3[pContact->m_iElem1.m_iBlock].m_aiEdgeNodes[pContact->m_iElem1.m_iNode]; // 返回单元里局部边编号下的边的行指针
            int *aiEdgeNodes2 = MyTri3[pContact->m_iElem2.m_iBlock].m_aiEdgeNodes[pContact->m_iElem2.m_iNode];

            int i1;
            int SpringCount = 2;
            int nEdgeNodes = 2; // 接触面涉及的连续节点数

            for (i1 = 0; i1 < SpringCount; i1++)
            {
                CbSpring *pSpring = &(pContact->m_aoSpring[i1]);
                pSpring->m_bUsed = true;
                pSpring->m_bBroken = false;
                pSpring->m_bHistoryBroken = false;
                pSpring->Brother1NodeNo = aiEdgeNodes1[i1]; // 单元的某个边的i1节点局部编号，其中aiEdgeNodes1为行首指针

                int iGlobalNode1 = aiNodes1[aiEdgeNodes1[i1]];
                for (int i2 = 0; i2 < nEdgeNodes; i2++)
                {
                    int iGlobalNode2 = aiNodes2[aiEdgeNodes2[i2]];
                    if (iGlobalNode2 == iGlobalNode1)
                    {
                        pSpring->Brother2NodeNo = aiEdgeNodes2[i2];
                        break;
                    }
                } // 寻找全局编号相同的另一个弹簧节点的局部编号和它所在行的标号（位置）

                pSpring->Area = 1; // 弹簧的截面积

                //==================================================aoSpringNum: 计算弹簧的数量？

                int *pNode1 = &aiNodes1[pSpring->Brother1NodeNo]; // 弹簧接触对中一个接触节点在单元节点编号数组中的指针
                int *pNode2 = &aiNodes2[pSpring->Brother2NodeNo];

                pSpring->Brother1GlobalNodeNo = *pNode1;
                pSpring->Brother2GlobalNodeNo = *pNode2;

                pSpring->Brother1SpringNodeNo = (iSpringNode + i1) * 2 + 1;
                pSpring->Brother2SpringNodeNo = (iSpringNode + i1) * 2 + 2;

                // printf("SpringNo1=%d, =====pNode1=%d\n", pSpring->Brother1SpringNodeNo, *pNode1);
                // printf("SpringNo2=%d, =====pNode2=%d\n", pSpring->Brother2SpringNodeNo, *pNode2);

            } // oMapSpringGroup里的迭代器：每个弹簧接触对的节点在单元局部编号下的指针，该弹簧接触对的节点在全局弹簧节点里的编号
            iSpringNode += 2;

            //===================================================================配置接触面的单位坐标系统

            double FaceCenter[2], Normal[2], Tangent[2];
            FaceCenter[0] = 0.0;
            FaceCenter[1] = 0.0;
            for (int iNode = 0; iNode < nEdgeNodes; iNode++)
            {
                FaceCenter[0] += afX[aiNodes1[aiVarNodes1[iNode]] - 1];
                FaceCenter[1] += afZ[aiNodes1[aiVarNodes1[iNode]] - 1];
            }

            FaceCenter[0] /= (double)nEdgeNodes;
            FaceCenter[1] /= (double)nEdgeNodes; // 中心点坐标

            Tangent[0] = afX[aiNodes1[aiVarNodes1[0]] - 1] - FaceCenter[0];
            Tangent[1] = afZ[aiNodes1[aiVarNodes1[0]] - 1] - FaceCenter[1]; // 切向向量

            Normal[0] = -Tangent[1];
            Normal[1] = Tangent[0]; // 依据相互垂直的两个向量，点积为0，计算出法向向量

            // 向量单位化
            norm(Tangent);
            norm(Normal);

            for (i1 = 0; i1 < 2; i1++)
            {
                pContact->LocalCoordSystem[0][i1] = Tangent[i1];
                pContact->LocalCoordSystem[1][i1] = Normal[i1];
            }
        }
        // printf("==============iter===========%d\n", iGlobalFace);
    } // 所有接触面的循环完毕

    //================================================================给所有有接触的节点进行编号 Generate ContactNodeGroup

    BegDem.m_aSpringForceGroup = (CbSpringForceGroup *)malloc(sizeof(CbSpringForceGroup) * nSprings);

    for (int i1 = 0; i1 < nSprings; i1++)
    {
        BegDem.m_aSpringForceGroup[i1].m_nSprings = 0;
    }

    int i1, i2;
    int SpringCount = 2;
    for (i1 = 0; i1 < nContacts; i1++)
    {
        pContact = &MyContact[i1];
        // printf("iContact=%d\n", i1);

        CbSpring *pSpring; // 弹簧类：弹簧物理量，表面积，刚度系数等
        for (i2 = 0; i2 < SpringCount; i2++)
        {
            pSpring = &(pContact->m_aoSpring[i2]);

            int iSForceGroup = 2 * i1 + i2; // 弹簧力组全局编号
            // printf("SpringForceGroup ID=%d\n", iSForceGroup);

            CbNodePos posNode1, posNode2;
            posNode1.m_iBlock = pContact->m_iElem1.m_iBlock;
            posNode2.m_iBlock = pContact->m_iElem2.m_iBlock;
            posNode1.m_iNode = pSpring->Brother1NodeNo;
            posNode2.m_iNode = pSpring->Brother2NodeNo;

            CbSpringForceGroup *pSForceGroup;
            pSForceGroup = &BegDem.m_aSpringForceGroup[iSForceGroup];

            pSForceGroup->m_aoSpringForces[pSForceGroup->m_nSprings].m_iBlock = i1;  // 某个节点的某个接触面上的弹簧力所属于的接触面在所有接触面内的编号
            pSForceGroup->m_aoSpringForces[pSForceGroup->m_nSprings].m_iSpring = i2; // 。。。弹簧力在所属接触面里的局部编号
            pSForceGroup->m_aoSpringForces[pSForceGroup->m_nSprings].m_iForce = 0;
            pSForceGroup->m_nSprings++;
            pSForceGroup->m_oNodePos = posNode1; // 某个节点的所属单元编号，该节点的单元内局部编号

            pSForceGroup->m_aoSpringForces[pSForceGroup->m_nSprings].m_iBlock = i1;
            pSForceGroup->m_aoSpringForces[pSForceGroup->m_nSprings].m_iSpring = i2;
            pSForceGroup->m_aoSpringForces[pSForceGroup->m_nSprings].m_iForce = 1;
            pSForceGroup->m_nSprings++;
            pSForceGroup->m_oNodePos = posNode2;

        } // m_aSpringForceGroup: 给所有有接触的节点进行编号—节点弹簧力组的编号和该组里的各个弹簧
    }

    printf("===================Set All Contacts===============\n");
    return 0;
}

int Build(int iContact, double fContactStiffFactor) // 计算接触面弹簧的法向刚度、切向刚度
{
    CbSpring *pSpring;
    int i1;
    int SpringCount = 2;

    double fArea = 0;
    for (i1 = 0; i1 < SpringCount; i1++)
    {
        pSpring = &(MyContact[iContact].m_aoSpring[i1]);

        fArea += pSpring->Area;
    }
    fArea = sqrt(fArea) / fContactStiffFactor;

    for (i1 = 0; i1 < SpringCount; i1++)
    {
        pSpring = &(MyContact[iContact].m_aoSpring[i1]);
        pSpring->NormalStiff = pSpring->Area * young / fArea; // kn = EA/L
        pSpring->ShearStiff = pSpring->Area * young / fArea / (1 + poisson);
        pSpring->m_fSpringLength = fArea;

        pSpring->m_afEForce[0][0] = 0.0;
        pSpring->m_afEForce[0][1] = 0.0;
        pSpring->m_afEForce[1][0] = 0.0;
        pSpring->m_afEForce[1][1] = 0.0;

        pSpring->SpringForce[0] = 0.0;
        pSpring->SpringForce[1] = 0.0;
        pSpring->Displace[0] = 0.0;
        pSpring->Displace[1] = 0.0;
        pSpring->rela_dis_local[0] = 0;
        pSpring->rela_dis_local[1] = 0;
    }

    return 0;
}