#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

#include "cephes.h"
#include "rats.h"

static const unsigned char ratsTable8_Hamming[256] = {
	0, 1, 1, 2, 1, 2, 2, 3,	1, 2, 2, 3, 2, 3, 3, 4,		//0000-XXXX
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

int Runs(double alpha, unsigned char *data, int bits, Runs_V *value)
{
	int		S, i, nb, nr, V;
	double	pi, erfc_arg, p_value;

	nb = bits / 8;
	nr = bits % 8;

	S = 0;
	for (i = 0; i < nb; i++)
		S += ratsTable8_Hamming[data[i]];
	if (nr)
		S += ratsTable8_Hamming[data[i] >> (8 - nr)] + 8 - nr;

	pi = (double)S / (double)bits;

	if (fabs(pi - 0.5) >(2.0 / sqrt(bits)))
		p_value = 0.0;
	else 
	{
		V = 1;
		if (nr)
		{
			for (i = 0; i < nb; i++)
				V += ratsTable8_Hamming[data[i] ^ (data[i] << 1 & 0xFE) ^ (data[i + 1] >> 7 & 0x01)];
			V += ratsTable8_Hamming[(data[i] ^ (data[i] << 1 & 0xFE)) >> (9 - nr)];
		}
		else
		{
			for (i = 0; i < nb - 1; i++)
				V += ratsTable8_Hamming[data[i] ^ (data[i] << 1 & 0xFE) ^ (data[i + 1] >> 7 & 0x01)];
			V += ratsTable8_Hamming[(data[i] ^ (data[i] << 1 & 0xFE))];
		}

		erfc_arg = fabs((double)V - 2.0 * bits * pi * (1 - pi)) / (2.0 * pi * (1 - pi) * sqrt(2 * bits));
		p_value = erfc(erfc_arg);
	}

	//TODO: return value & ¿ÉÓÅ»¯
	if (value)
	{
		value->p_value = p_value;
		value->runs = V;
		value->pi = pi;
	}

	if (p_value < alpha)
		return 0;
	return 1;
}
