/* MonobitFrequency.c
 *
 * Copyright (C) 2016-2018 rats Inc.
 *
 * ITEM-1: MonobitFrequency Test
 *
 * This file is part of rats.
 */
#include <math.h>

#include "cephes.h"
#include "rats.h"

#ifndef SQRT2
#define SQRT2		1.41421356237309504880
#endif // !SQRT2

//8-bits Hamming Weight Table
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

int MonobitFrequency(double alpha, unsigned char *data, int bits, 
	MonobitFrequency_V *value)
{
	double p_value, s_obs;
	int n, sum;
	int i, p;

	//Step 1: Check input size
	if (bits < 100)
		return -1;//Error: Input size error.
	n = bits;

	//Step 2: Get statistic
	sum = 0;
	p = n / 8;
	for (i = 0; i < p; i++)
		sum += ratsTable8_Hamming[data[i]];
	if (p = n % 8)
		sum += ratsTable8_Hamming[data[i] >> (8 - p)];
	sum = 2 * sum - n;

	s_obs = fabs((double)sum) / sqrt(n);
	p_value = erfc(s_obs / SQRT2);

	if (value)
	{
		value->p_value = p_value;
		value->s_obs = s_obs;
		value->sum = sum;
	}

	//Step 3: Decision
	if (p_value < alpha)
		return 0;
	return 1;
}

//void MonobitFrequency_demo()
//{
//	unsigned char a[] = { 0x5f, 0x78 };
//	//unsigned char a[] = { 0x5f,0x79,0x63,0xd8,0x87,0x40,0x8a,0x9d,0xc2,0x4a,0xa8,0xe3,0x94,0x1f, 0x7f, 0x88, 0xca };
//	//unsigned char a[] = "000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000";
//	//unsigned char a[] = "erq823alsfjqo;wjfao;eiufalsfj;lsdfwr632/rwaw4rqur88 qwfo lofjqwo;-3=\"2wef]qwfqwr";
//
//	MonobitFrequency(0, a, (sizeof(a)) * 8, nullptr);
//}

//int main()
//{
//	//Autocorrelation_Init8bitWeightTable();
//	MonobitFrequency_demo();
//	system("pause");
//	return 0;
//}