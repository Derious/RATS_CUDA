/* BlockFrequency.c
 *
 * Copyright (C) 2016-2018 rats Inc.
 *
 * ITEM-2: BlockFrequency Test
 *
 * This file is part of rats.
 */
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "cephes.h"
#include "rats.h"

static const unsigned char ratsTable8_Hamming[256] = {
	0, 1, 1, 2, 1, 2, 2, 3, 1, 2, 2, 3, 2, 3, 3, 4,		//0000-XXXX
	1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5,		//0001-XXXX
	1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5,		//0010-XXXX
	2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,		//0011-XXXX
	1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5,		//0100-XXXX
	2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,		//0101-XXXX
	2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,		//0110-XXXX
	3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7,		//0111-XXXX
	1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5,		//1000-XXXX
	2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,		//1001-XXXX
	2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,		//1010-XXXX
	3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7,		//1011-XXXX
	2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,		//1100-XXXX
	3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7,		//1101-XXXX
	3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7,		//1110-XXXX
	4, 5, 5, 6, 5, 6, 6, 7, 5, 6, 6, 7, 6, 7, 7, 8		//1111-XXXX
};

int BlockFrequency(double alpha, unsigned char *data, int bits,
	int M, BlocktFrequency_V *value)
{
	double	p_value, v_obs, sum;
	int		i, n, N, blockSum;
	double	pi, v;
	int		bi, br, blockRes;
	int totalSum;

	//Step 1: Check input size
	n = bits;
	if (n < 100 || M < 20)
		return -1;//Error: Input size error.

	//Step 2: Get statistic
	N = bits / M; 		/* # OF SUBSTRING BLOCKS      */
	sum = 0.0;
	totalSum = 0;

	for (i = 0; i < N; i++)
	{
		blockSum = 0;
		blockRes = M;

		bi = (i * M) / 8;
		if (br = (i * M) % 8)
		{
			blockSum += ratsTable8_Hamming[(data[bi++] << br) & 0xFF];
			blockRes -= (8 - br);
		}
		br = bi + blockRes / 8;
		while (bi < br)
		{
			blockSum += ratsTable8_Hamming[data[bi++]];
		}
		if (br = blockRes % 8)
		{
			blockSum += ratsTable8_Hamming[(data[bi++] >> (8 - br)) & 0xFF];
		}

		pi = (double)blockSum / (double)M;
		v = pi - 0.5;
		sum += v * v;
		totalSum += blockSum;
	}
	v_obs = 4.0 * M * sum;
	p_value = cephes_igamc(N / 2.0, v_obs / 2.0);

	if (value)
	{
		value->p_value = p_value;
		value->v_obs = v_obs;
	}

	//Step 3: Decision
	if (p_value < alpha)
		return 0;
	return 1;
}
