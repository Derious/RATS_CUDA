/* Poker.c
 *
 * Copyright (C) 2016-2018 rats Inc.
 *
 * ITEM-3: Poker Test
 *
 * This file is part of rats.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "cephes.h"
#include "rats.h"

//Get the frequency of all possible non-overlapping m-bit patterns
void rats_get_poker_frequency(unsigned char *data, int n, int m, int *v)
{
	int i, N, v_size;

	N = n / m;

	v_size = 1 << m;
	memset(v, 0, v_size * sizeof(int));

	if (m == 2)
	{
		for (i = 0; i < N/4; i++)
		{
			v[(data[i] >> 6) & 0x03]++;
			v[(data[i] >> 4) & 0x03]++;
			v[(data[i] >> 2) & 0x03]++;
			v[(data[i]     ) & 0x03]++;
		}
		if (N % 4)
		{
			if (N % 4 == 3)
			{
				v[(data[i] >> 6) & 0x03]++;
				v[(data[i] >> 4) & 0x03]++;
				v[(data[i] >> 2) & 0x03]++;
			}
			else if (N % 4 == 2)
			{
				v[(data[i] >> 6) & 0x03]++;
				v[(data[i] >> 4) & 0x03]++;
			}
			else
			{
				v[(data[i] >> 6) & 0x03]++;
			}
		}
	}
	else if (m == 4)
	{
		for (i = 0; i < N/2; i++)
		{
			v[(data[i] >> 4) & 0x0F]++;
			v[(data[i]     ) & 0x0F]++;
		}
		if (N % 2)
		{
			v[(data[i] >> 4) & 0x0F]++;
		}
	}
	else if (m == 8)
	{
		for (i = 0; i < N; i++)
		{
			v[data[i]]++;
		}
	}
	else
		return;
}

int Poker(double alpha, unsigned char *data, int bits, int m, Poker_V *value)
{
	double p_value, v_obs, sum;
	unsigned int *ni, ni_size;
	int i, p, N;

#ifdef _DEBUG
	N = 0;
#else
	N = bits / m;
#endif // _DEBUG

	ni_size = 1 << m;
	if ((ni = (unsigned int*)calloc(ni_size, sizeof(unsigned int))) == NULL)
		return -1;

	rats_get_poker_frequency(data, bits, m, ni);
	p = 1 << m;
	for (i = 0, sum = 0; i < p; i++)
	{
		sum += (double)ni[i] * (double)ni[i];
#ifdef _DEBUG
		N += ni[i];
#endif // _DEBUG
	}
	free(ni);
#ifdef _DEBUG
	if (N != bits / m)
		return -1;//统计失败
#endif // _DEBUG

	v_obs = (double)p * sum / N - N;
	p_value = cephes_igamc((p - 1) / 2.0, v_obs / 2.0);

	if (value)
	{
		value->p_value = p_value;
		value->v_obs = v_obs;
		value->N = N;
	}

	if (p_value < alpha)
		return 0;//未通过
	return 1;//通过
}

// void Poker_demo()
// {
// 	unsigned char a[] = { 0x5f, 0x79, 0x63, 0xd8, 0x87, 0x40, 0x8a, 0x9d, 0xc2, 0x4a, 0xa8, 0xe3, 0x94, 0x1f, 0x7f, 0x88, 0xca };
// 	//unsigned char a[] = "000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000";
// 	//Poker_P param;
// 	// param.m = 4;
// 	// Poker(0, a, (sizeof(a)-1) * 8, 4, 0);
// }

// int main()
// {       
// 	Poker_demo();
// 	system("pause");
// 	return 0;
// }