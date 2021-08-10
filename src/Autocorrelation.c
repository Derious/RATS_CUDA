
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "cephes.h"
#include "rats.h"



#define SQRT2		1.41421356237309504880
const unsigned char Autocorrelation_HammingWeight8bit[256] = {
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
//
//void Autocorrelation_Init8bitWeightTable()
//{
//	int i, j;
//
//	FILE			*fp;
//
//	if ((fp = fopen("tempfile/HammingWeight8bit.cpp", "w")) == NULL)
//	{
//		printf("cannot open file!/n");
//		exit(0);
//	}
//
//	fprintf(fp, "unsigned char HammingWeight8bit[256]={\n");
//	for (i = 0; i < 256; i++)
//	{
//		unsigned char index = i & 0xff;
//		unsigned char weight = 0;
//		for (j = 0; j < 8; j++)
//		{
//			if (index & 0x1)
//			{
//				weight++;
//			}
//			index >>= 1;
//		}
//
//		fprintf(fp, "0x%02X, ", weight);
//		if (i % 16 == 15)
//		{
//			fprintf(fp, "\n");
//		}
//	}
//
//	fprintf(fp, "};\n\n");
//	fclose(fp);
//}
//
//void Autocorrelation_CalA816(unsigned char *data, int bits, unsigned int d, unsigned int *A)
//{
//	int ucd = d / 8;
//	int len_data = bits / 8;
//	int len_acr = len_data - ucd;
//	unsigned char *t = (unsigned char *)malloc(len_acr * sizeof(unsigned char));
//
//	for (int i = 0; i < len_acr; i++)
//	{
//		t[i] = data[i] ^ data[i + ucd];
//		*A += Autocorrelation_HammingWeight8bit[t[i]];
//	}
//}
//
//void Autocorrelation_CalA12(unsigned char *data, int bits, unsigned int d, unsigned int *A)
//{
//	int i;
//	unsigned char getbit[8] =
//	{
//		//0b00000000,
//		0x80,
//		0x40,
//		0x20,
//		0x10,
//		0x08,
//		0x04,
//		0x02,
//		0x01
//	};
//	for (i = 0; i + d < bits; i++)
//	{
//		int index0 = i / 8;
//		int pos0 = i % 8;
//		int index1 = (i + d) / 8;
//		int pos1 = (i + d) % 8;
//		unsigned char xorresult = ((data[index0] & getbit[pos0]) >> (7 - pos0)) ^ ((data[index1] & getbit[pos1]) >> (7 - pos1));
//		if (xorresult)
//		{
//			*A++;
//		}
//	}
//}

//void Autocorrelation_CalA(unsigned char *data, int bits, unsigned int d, unsigned int &A)
//{
//
//}

int Autocorrelation(double alpha, unsigned char *data, int bits, int d, Autocorrelation_V *value)
{
	double p_value, v_obs;
	unsigned int A;
	int bytes, rbits, di, dj;
	int i, n;

	if (d < 0)
		return -1;		//ERROR: Paramter error.

	bytes = (bits - d) / 8;
	rbits = (bits - d) % 8;
	di = d / 8;
	dj = d % 8;

	A = 0;
	if (di == 0)
	{
		for (i = 0; i < bytes; i++)
			A += Autocorrelation_HammingWeight8bit[(data[i] ^ (data[i] << dj & 0xFF) ^ (data[i + 1] >> (8 - dj)))];
		if (rbits)
			A += Autocorrelation_HammingWeight8bit[(data[i] ^ (data[i] << dj & 0xFF) ^ (data[i + 1] >> (8 - dj))) >> (8 - rbits)];
	}
	else if (dj == 0)
	{
		for (i = 0; i < bytes; i++)
			A += Autocorrelation_HammingWeight8bit[(data[i] ^ data[i + di])];
		if (rbits)
			A += Autocorrelation_HammingWeight8bit[(data[i] ^ data[i + di]) >> (8 - rbits)];
	}
	else
	{
		for (i = 0; i < bytes; i++)
			A += Autocorrelation_HammingWeight8bit[(data[i] ^ (data[i + di] << dj & 0xFF) ^ (data[i + di + 1] >> (8 - dj)))];
		if (rbits)
			A += Autocorrelation_HammingWeight8bit[(data[i] ^ (data[i + di] << dj & 0xFF) ^ (data[i + di + 1] >> (8 - dj))) >> (8 - rbits)];
	}

	v_obs = (2 * (double)A - (double)(bits - d)) / sqrt((double)(bits - d));
	p_value = erfc(fabs(v_obs) / SQRT2);

	if (value)
	{
		value->p_value = p_value;
		value->v_obs = v_obs;
		value->A = A;
	}

	if (p_value < alpha)
		return 0;//未通过
	return 1;//通过
}

void Autocorrelation_demo()
{
	unsigned char a[] = { 0x5f, 0x79 };
	//unsigned char a[] = { 0x5f,0x79,0x63,0xd8,0x87,0x40,0x8a,0x9d,0xc2,0x4a,0xa8,0xe3,0x94,0x1f, 0x7f, 0x88, 0xca };
	//unsigned char a[] = "000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000";
	//unsigned char a[] = "erq823alsfjqo;wjfao;eiufalsfj;lsdfwr632/rwaw4rqur88 qwfo lofjqwo;-3=\"2wef]qwfqwr";

	Autocorrelation(0, a, (sizeof(a)) * 8, 1, NULL);
}

//int main()
//{
//	//Autocorrelation_Init8bitWeightTable();
//	Autocorrelation_demo();
//	system("pause");
//	return 0;
//}