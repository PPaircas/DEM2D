#include "DemRun.h"
#include "Element.h"
#include <stdbool.h>
#include <math.h>

#define DEBUG 4

////////////节点计算//////////
int DoNode_Inc(int threadIdx)
{
	int inode, ii;
	int Tri3Nodes = 3;

	for (inode = 0; inode < Tri3Nodes; inode++)
	{
		double ExtForce = 0.0, UnbalForce = 0.0;
		int iGlobalNo = MyTri3[threadIdx].m_agNodes[inode];
		CbNode *pNode = &(MyTri3[threadIdx].aNodes[inode]);

		for (ii = 0; ii < 2; ii++)
		{
			double UnbalanceForce = pNode->m_afExtForce[ii] - pNode->m_afForce[ii];

			if (DEBUG == 1 && UnbalanceForce != 0)
			{
				printf("ExtF=%e, InF=%e, unbalanceforce=%e\n", pNode->m_afExtForce[ii], pNode->m_afForce[ii], UnbalanceForce);
			}

			ExtForce += fabs(pNode->m_afExtForce[ii]);
			UnbalForce += fabs(UnbalanceForce);

			fA[iGlobalNo][ii] = UnbalanceForce / pNode->m_fMass;
			//fV[iGlobalNo][ii] += fA[iGlobalNo][ii] * dt;
			fVV[iGlobalNo][ii] = fV[iGlobalNo][ii] + fA[iGlobalNo][ii] * dt;

			pNode->m_afDeltaDisplace[ii] = fVV[iGlobalNo][ii] * dt;
			//fU[iGlobalNo][ii] += pNode->m_afDeltaDisplace[ii];
			fUU[iGlobalNo][ii] = fU[iGlobalNo][ii] + pNode->m_afDeltaDisplace[ii];

			if (DEBUG == 4 && fU[iGlobalNo][ii])
			{
				printf(" iGlobalNo=%d, CoordX=%f, CoordZ=%f, displace=%e\n", iGlobalNo, dfX[iGlobalNo], dfZ[iGlobalNo], fU[iGlobalNo][ii]);
			}

			// pNode->m_afForce[ii] = 0.0;

			// 更新单元节点的位移、速度、加速度
			if (ii == 0) // X分量
			{
				fUX[iGlobalNo] = fUU[iGlobalNo][ii];
				fVX[iGlobalNo] = fVV[iGlobalNo][ii];
			}
			else // Z分量
			{
				fUZ[iGlobalNo] = fUU[iGlobalNo][ii];
				fVZ[iGlobalNo] = fVV[iGlobalNo][ii];
			}
		}
		m_fExtForce += ExtForce;	 // 所有单元上所有节点的外力绝对值之和
		m_fUnbalForce += UnbalForce; // 所有单元上所有节点的不平衡力绝对值之和
	}

	//printf("threadId=%d, m_fExtForce=%e, m_fUnbalForce=%e\n", threadIdx, m_fExtForce, m_fUnbalForce);

	return 0;
}

////////////B计算////////////////
int DoBlock_B_Inc(int threadIdx)
{
	if (MyTri3[threadIdx].m_iElemModel == 1) // linear model
	{
		DoNode_Inc(threadIdx); // 计算单元节点加速度，速度，位移

		CalculateDeltaStress(threadIdx); // 根据单元节点位移增量计算单元应力增量和应变增量

		CalculateNodeForce(threadIdx); // 	根据单元节点位移计算单元节点力

		// CalculateAverStressStrain(threadIdx); // 计算单元平均应力应变

		// CalculatePrinStressStrain(threadIdx); // 计算单元主应力主应变及方向
	}

	return 0;
}