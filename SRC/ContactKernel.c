#include "DemRun.h"
#include "Element.h"

#define verbose 2

// 计算接触面上的弹簧力
int DoContact_Linear_Inc(int threadIdx)
{
    int i1 = threadIdx;
    int iSpring, nSpringCount;
    int ip1, ip2;
    double aa[2][2], NormalStiff, ShearStiff;

    MyElem *pElem1, *pElem2;

    CbContact *pContact = &MyContact[i1];

    pElem1 = &MyTri3[pContact->m_iElem1.m_iBlock];
    pElem2 = &MyTri3[pContact->m_iElem2.m_iBlock];

    // 接触面局部坐标系统
    aa[0][0] = pContact->LocalCoordSystem[0][0];
    aa[0][1] = pContact->LocalCoordSystem[0][1];
    aa[1][0] = pContact->LocalCoordSystem[1][0];
    aa[1][1] = pContact->LocalCoordSystem[1][1];

    nSpringCount = 2; // 接触面上弹簧个数

    for (iSpring = 0; iSpring < nSpringCount; iSpring++) // 接触面上弹簧组循环
    {
        ip1 = pContact->m_aoSpring[iSpring].Brother1NodeNo;
        ip2 = pContact->m_aoSpring[iSpring].Brother2NodeNo;
        NormalStiff = pContact->m_aoSpring[iSpring].NormalStiff; // 法向弹簧刚度系数
        ShearStiff = pContact->m_aoSpring[iSpring].ShearStiff;   // 切向弹簧刚度系数

        // 计算对应节点间的接触力
        double Cb[2], ca[2];
        Cb[0] = pElem1->aNodes[ip1].m_afDeltaDisplace[0] - pElem2->aNodes[ip2].m_afDeltaDisplace[0]; // 弹簧伸长量X分量
        Cb[1] = pElem1->aNodes[ip1].m_afDeltaDisplace[1] - pElem2->aNodes[ip2].m_afDeltaDisplace[1]; // 弹簧伸长量Z分量

        if (verbose == 2)
        {
            if (Cb[0] != 0 || Cb[1] != 0)
                printf("iContact=%d, Cb[0]=%e, Cb[1]=%e\n", i1, Cb[0], Cb[1]);
        }

        ca[0] = aa[0][0] * Cb[0] + aa[0][1] * Cb[1]; // 全局转到局部
        ca[1] = aa[1][0] * Cb[0] + aa[1][1] * Cb[1];

        pContact->m_aoSpring[iSpring].SpringForce[0] += ca[0] * ShearStiff; // 弹簧力
        pContact->m_aoSpring[iSpring].SpringForce[1] += ca[1] * NormalStiff;

        pContact->m_aoSpring[iSpring].m_bBroken = 0;

        ///////////////////////////////////把弹簧力转到全局坐标////////////弹簧力作用两侧的单元节点

        pContact->m_aoSpring[iSpring].m_afEForce[1][0] = aa[0][0] * pContact->m_aoSpring[iSpring].SpringForce[0] + aa[1][0] * pContact->m_aoSpring[iSpring].SpringForce[1];
        pContact->m_aoSpring[iSpring].m_afEForce[1][1] = aa[0][1] * pContact->m_aoSpring[iSpring].SpringForce[0] + aa[1][1] * pContact->m_aoSpring[iSpring].SpringForce[1];

        pContact->m_aoSpring[iSpring].m_afEForce[0][0] = -pContact->m_aoSpring[iSpring].m_afEForce[1][0];
        pContact->m_aoSpring[iSpring].m_afEForce[0][1] = -pContact->m_aoSpring[iSpring].m_afEForce[1][1];
    }

    return 0;
}

/// 将接触面上弹簧力作用到两侧单元节点
int SyncSpringGroup(int threadIdx)
{
    int i1 = threadIdx;
    CbSpringForceGroup rGroup = BegDem.m_aSpringForceGroup[i1];
    for (int i2 = 0; i2 < rGroup.m_nSprings; i2++)
    {
        CbNodePos rPos = rGroup.m_oNodePos;
        CbForcePos rForcePos = rGroup.m_aoSpringForces[i2];
        CbSpring pSpring = MyContact[rForcePos.m_iBlock].m_aoSpring[rForcePos.m_iSpring];

        double *rForce = pSpring.m_afEForce[rForcePos.m_iForce];
        CbNode *pNode = &(MyTri3[rPos.m_iBlock].aNodes[rPos.m_iNode]);

        pNode->m_afForce[0] -= rForce[0]; // 与节点力方向相反
        pNode->m_afForce[1] -= rForce[1];
    }

    return 0;
}
