#include <stdio.h>
#include <stdlib.h>
#include <memory.h>
#include <math.h>

#include "cephes.h"
#include "rats.h"

#define LEN 1024


//void RunsDistribution_InitRunsArr(unsigned int *Runs0, unsigned int *Runs1)
//{
//	for (int i = 0; i < LEN; i++)
//	{
//		Runs0[i] = 0;
//		Runs1[i] = 0;
//	}
//}
//
//void RunsDistribution_Calek(int bits, unsigned int *k, double *e)
//{
//	int i;
//	for (i = 1;; i++)
//	{
//		e[i] = ((double)bits - i + 3) / pow((double)2, i + 2);
//		if (e[i] < 5)
//		{
//			break;
//		}
//	}
//	*k = i - 1;
//}

void RunsDistribution_CountRuns(unsigned char *data, int bits, unsigned int *Runs0, unsigned int *Runs1)
{
	int i, j;
	int count = 0;
	int currentserial = 0;

	for (i = 0; bits - i >= 8; i += 8)
	{
		unsigned char tmp = data[i >> 3];
		for (j = 0; j < 8; j++)
		{
			if (tmp & 0x80)		/* 读到1 */
			{
				if (currentserial == 1)	/* 且上一个也是1 */
				{
					count++;
				}
				else					/* 但上一个是0 */
				{
					Runs0[count]++;
					currentserial = 1;
					count = 1;
				}
			}
			else						/* 读到0 */
			{
				if (currentserial == 0)	/* 且上一个也是0 */
				{
					count++;
				}
				else					/* 但上一个是1 */
				{
					Runs1[count]++;
					currentserial = 0;
					count = 1;
				}
			}
			tmp <<= 1;
		}
	}
	/* 处理最后一个 */
	if (currentserial == 0)
	{
		Runs0[count]++;
	}
	else
	{
		Runs1[count]++;
	}
}

//void RunsDistribution_CalV(unsigned int k, double *e, unsigned int *Runs0, unsigned int *Runs1, double &V)
//{
//	for (int i = 1; i <= k; i++)
//	{
//		double v0 = (double)Runs0[i] - e[i];
//		double v1 = (double)Runs1[i] - e[i];
//		V += (v0 * v0 + v1 * v1) / e[i];
//	}
//}

int RunsDistribution(double alpha, unsigned char *data, int bits, RunsDistribution_V *value)
{
	int i, k;								//游程数取值: 1-k
	double e[32];								//从[1]开始存储
	unsigned int Runs0[LEN];					//从[1]开始存储
	unsigned int Runs1[LEN];					//从[1]开始存储
	double p_value, v_obs, v0, v1;

	//printf("bits: %d\n", bits);

	for (i = 1; i < 32; i++)
	{
		if ((e[i] = ((double)bits - i + 3) / pow((double)2, i + 2)) < 5)
		{
			k = i - 1;
			break;
		}
	}

	memset(Runs0, 0, sizeof(unsigned int)* (k + 1));
	memset(Runs1, 0, sizeof(unsigned int)* (k + 1));

	//RunsDistribution_Calek(bits, &k, e);
	//printf("k: %d\n\n", k);
	//for (int i = 1; i <= k + 1; i++)
	//{
	//	printf("%f\n", e[i]);
	//}
	//printf("\n");

	RunsDistribution_CountRuns(data, bits, Runs0, Runs1);

	//for (int i = 1; i <= 10; i++)
	//{
	//	printf("%d, %d\n", Runs0[i], Runs1[i]);
	//}
	//printf("last line not count\n\n");

	v_obs = 0;
	for (i = 1; i <= k; i++)
	{
		v0 = (double)Runs0[i] - e[i];
		v1 = (double)Runs1[i] - e[i];
		v_obs += (v0 * v0 + v1 * v1) / e[i];
	}
	p_value = cephes_igamc((double)(k - 1), v_obs / 2.0);

	//RunsDistribution_CalV(k, e, Runs0, Runs1, V);
	//printf("V: %f\n", V);

	
	//printf("P-value: %f\n", p_value);

	if (value)
	{
		value->k = k;
		value->v_obs = v_obs;
		value->p_value = p_value;
		memcpy(value->e, e, sizeof(double)* (k + 1));
		memcpy(value->Runs0, Runs0, sizeof(unsigned int)* (k + 1));
		memcpy(value->Runs1, Runs1, sizeof(unsigned int)* (k + 1));
	}

	if (p_value < alpha)
		return 0;//未通过
	return 1;//通过
}

void RunsDistribution_demo()
{
	//unsigned char a[] = { 0x5f, 0x79 };
	unsigned char a[] = { 0x5f, 0x79, 0x63, 0xd8, 0x87, 0x40, 0x8a, 0x9d, 0xc2, 0x4a, 0xa8, 0xe3, 0x94, 0x1f, 0x7f, 0x88, 0xca };
	//unsigned char a[] = "000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000";
	//unsigned char a[] = "erq823alsfjqo;wjfao;eiufalsfj;lsdfwr632/rwaw4rqur88 qwfo lofjqwo;-3=\"2wef]qwfqwr";

	RunsDistribution(0, a, (sizeof(a)) * 8, NULL);
}

//int main()
//{
//	RunsDistribution_demo();
//	system("pause");
//	return 0;
//}