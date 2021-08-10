/* Serial.c
*
* Copyright (C) 2016-2018 rats Inc.
*
* ITEM-4: Serial Test
*
* This file is part of rats.
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "cephes.h"  
#include "rats.h"

#ifndef GET_EPSILON
#define GET_EPSILON(d, i)	(((d)[(i) / 8] >> (7 - (i) % 8)) & 1)
#endif // !GET_EPSILON

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

//Get the frequency of all possible overlapping m-bit patterns
void rats_get_serial_frequency(unsigned char *data, int n, int m, int *v)
{
	int i, p, block;

	if (m <= 0 || m > 31)
	{
		v[0] = 0;
		return;
	}
	if (m == 1)
	{
		v[1] = 0;
		p = n / 8;
		for (i = 0; i < p; i++)
			v[1] += ratsTable8_Hamming[data[i]];
		if (p = n % 8)
			v[1] += ratsTable8_Hamming[data[i] >> (8 - p)];
		v[0] = n - v[1];
		return;
	}

	p = (1 << m) - 1;
	for (i = 0; i <= p; i++)
		v[i] = 0;

	//Step-1: Get last 32-bit
	if (n % 8)
	{
		i = n / 8 - 4;
		block = (unsigned int)data[i] << 24 | (unsigned int)data[i + 1] << 16 | (unsigned int)data[i + 2] << 8 | (unsigned int)data[i + 3];
		block = (block << (n % 8)) | ((unsigned int)data[i + 1] >> (8 - (n % 8)));
	}
	else
	{
		i = n / 8 - 4;
		block = (unsigned int)data[i] << 24 | (unsigned int)data[i+1] << 16 | (unsigned int)data[i+2] << 8 | (unsigned int)data[i+3];
	}

	//Step-2: Get blocks count
	for (i = 0; i < n; i++)
	{
		block = (block << 1) + GET_EPSILON(data, i);
		v[block & p]++;
	}

	////Step-1: Get first block
	//if (m <= 8)
	//	block = (unsigned int)data[0] >> (8 - m);
	//else if (m <= 16)
	//	block = ((unsigned int)data[0] << 8 | (unsigned int)data[1]) >> (16 - m);
	//else if (m <= 24)
	//	block = ((unsigned int)data[0] << 16 | (unsigned int)data[1] << 8 | (unsigned int)data[2]) >> (24 - m);
	//else
	//	block = ((unsigned int)data[0] << 24 | (unsigned int)data[1] << 16 | (unsigned int)data[2] << 8 | (unsigned int)data[3]) >> (32 - m);
	//v[block]++;

	////Step-2: Get other blocks
	//for (i = m; i < n; i++)
	//{
	//	block = ((block << 1) + ((data[i / 8] >> (7 - i % 8)) & 1)) & p;
	//	v[block]++;
	//}
	//for (i = 0; i < m - 1; i++)
	//{
	//	block = ((block << 1) + ((data[i / 8] >> (7 - i % 8)) & 1)) & p;
	//	v[block]++;
	//}
}

int Serial(double alpha, unsigned char *data, int bits, int m, Serial_V *value)
{
	double	p_value1, p_value2, psim0, psim1, psim2, del1, del2;
	int i, p, n;
	unsigned int *Vm, Vm_size;

	//Step 1: Check input size
	if ((n = bits) < 100)
		return -1;//Error: Input size error.
	if ((m < 2) || (m >= log(n)/log(2) - 2))
		return -1;//Error: Input size error.

	Vm_size = 1 << m;
	if ((Vm = (unsigned int*)calloc(Vm_size, sizeof(unsigned int))) == NULL)
		return -2;//Error: Memery alloc error.

	//Step 2: Get statistic
	rats_get_serial_frequency(data, n, m, Vm);
	p = 1 << m;
	for (psim0 = 0.0, i = 0; i < p; i++)
		psim0 += (double)Vm[i] * (double)Vm[i];
	psim0 = (psim0 * (double)p / (double)n) - (double)n;

	rats_get_serial_frequency(data, n, m - 1, Vm);
	p = 1 << (m - 1);
	for (psim1 = 0.0, i = 0; i < p; i++)
		psim1 += (double)Vm[i] * (double)Vm[i];
	psim1 = (psim1 * (double)p / (double)n) - (double)n;

	if (m <= 2)
	{
		psim2 = 0.0;
	}
	else
	{
		rats_get_serial_frequency(data, n, m - 2, Vm);
		p = 1 << (m - 2);
		for (psim2 = 0.0, i = 0; i < p; i++)
			psim2 += (double)Vm[i] * (double)Vm[i];
		psim2 = (psim2 * (double)p / (double)n) - (double)n;
	}

	free(Vm);

	del1 = psim0 - psim1;
	del2 = psim0 - 2.0*psim1 + psim2;
	p_value1 = cephes_igamc(pow(2, m - 2), del1 / 2.0);
	p_value2 = cephes_igamc(pow(2, m - 3), del2 / 2.0);

	//TODO: return value & ¿ÉÓÅ»¯
	if (value)
	{
		value->p_value1 = p_value1;
		value->p_value2 = p_value2;
		value->del1 = del1;
		value->del2 = del2;
	}

	//Step 3: Decision
	if ((p_value1 < alpha) || (p_value2 < alpha))
		return 0;
	return 1;
}
