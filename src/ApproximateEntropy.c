#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#include "cephes.h"  
#include "rats.h"

#ifndef GET_EPSILON
#define GET_EPSILON(d, i)	(((d)[(i) / 8] >> (7 - (i) % 8)) & 1)
#endif // !GET_EPSILON

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
I T E M ---- A P P R O X I M A T E  E N T R O P Y   T E S T
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */



int ApproximateEntropy(double alpha, unsigned char *data, int bits, int m, ApproximateEntropy_V *value)
{
	int				i, j, k, r, blockSize, powLen, index;
	double			sum, numOfBlocks, ApEn[2], apen, v_obs, p_value;
	unsigned int	*P;

	//fprintf(stats[TEST_APEN], "\t\t\tAPPROXIMATE ENTROPY TEST\n");
	//fprintf(stats[TEST_APEN], "\t\t--------------------------------------------\n");
	//fprintf(stats[TEST_APEN], "\t\tCOMPUTATIONAL INFORMATION:\n");
	//fprintf(stats[TEST_APEN], "\t\t--------------------------------------------\n");
	//fprintf(stats[TEST_APEN], "\t\t(a) m (block length)    = %d\n", m);

	if (m > (int)(log(bits) / log(2) - 5)) 
		return -2;		//ERROR: 'm' exceeds recommended value.

	r = 0;

	for (blockSize = m; blockSize <= m + 1; blockSize++) {
		if (blockSize == 0) {
			ApEn[0] = 0.00;
			r++;
		}
		else {
			numOfBlocks = (double)bits;
			powLen = (int)pow(2, blockSize + 1) - 1;
			if ((P = (unsigned int*)calloc(powLen, sizeof(unsigned int))) == NULL)
				return -1;		//ERROR: Insufficient memory available.

			for (i = 1; i<powLen - 1; i++)
				P[i] = 0;
			for (i = 0; i < bits; i++) { /* COMPUTE FREQUENCY */
				k = 1;
				for (j = 0; j<blockSize; j++) {
					k <<= 1;
					if (GET_EPSILON(data, (i + j) % bits) == 1)
						k++;
				}
				P[k - 1]++;
			}
			/* DISPLAY FREQUENCY */
			sum = 0.0;
			index = (int)pow(2, blockSize) - 1;
			for (i = 0; i<(int)pow(2, blockSize); i++) {
				if (P[index] > 0)
					sum += P[index] * log(P[index] / numOfBlocks);
				index++;
			}
			sum /= numOfBlocks;
			ApEn[r] = sum;
			r++;
			free(P);
		}
	}
	apen = ApEn[0] - ApEn[1];

	v_obs = 2.0*bits*(log(2) - apen);
	p_value = cephes_igamc(pow(2, m - 1), v_obs / 2.0);

	//TODO: return value & ø…”≈ªØ

	if (p_value < alpha)
		return 0;
	return 1;
}