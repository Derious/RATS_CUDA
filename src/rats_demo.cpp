#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "rats.h"
#include "rats_app.h"

#define BITS	1000000


int rats_app_sts(unsigned char *data, int bits)
{
	double alpha = 0.01;
	clock_t time_s, time_e;
	//if (bits < 1000000)
	//	return 0;
	//bits = 1000000;

	//1-频数检测
	time_s = clock();
	if (!MonobitFrequency(alpha, data, bits, NULL))			return 0;
	time_e = clock();
	printf("RATS test flag is %d, cost time %.4f.\n", 1, (double)(time_e - time_s) / CLOCKS_PER_SEC);
	//2-块内频数检测
	time_s = clock();
	if (!BlockFrequency(alpha, data, bits, 128, NULL))		return 0;
	if (!BlockFrequency(alpha, data, bits, 100, NULL))		return 0;
	time_e = clock();
	printf("RATS test flag is %d, cost time %.4f.\n", 2, (double)(time_e - time_s) / CLOCKS_PER_SEC);
	//3-扑克检测
	time_s = clock();

	if (!Poker(alpha, data, bits, 4, NULL))					return 0;
	if (!Poker(alpha, data, bits, 8, NULL))					return 0;
	time_e = clock();
	printf("RATS test flag is %d, cost time %.4f.\n", 3, (double)(time_e - time_s) / CLOCKS_PER_SEC);
	//4-序列检测
	time_s = clock();
	if (!Serial(alpha, data, bits, 2, NULL))				return 0;
	if (!Serial(alpha, data, bits, 5, NULL))				return 0;
	if (!Serial(alpha, data, bits, 16, NULL))				return 0;
	time_e = clock();
	printf("RATS test flag is %d, cost time %.4f.\n", 4, (double)(time_e - time_s) / CLOCKS_PER_SEC);
	//5-游程总数检测
	time_s = clock();
	if (!Runs(alpha, data, bits, NULL))						return 0;
	time_e = clock();
	printf("RATS test flag is %d, cost time %.4f.\n", 5, (double)(time_e - time_s) / CLOCKS_PER_SEC);
	////6-游程分布检测
	time_s = clock();
	if (!RunsDistribution(alpha, data, bits, NULL))			return 0;
	time_e = clock();
	printf("RATS test flag is %d, cost time %.4f.\n", 6, (double)(time_e - time_s) / CLOCKS_PER_SEC);
	//7-块内最大1游程检测
	time_s = clock();
	if (!LongestRun(alpha, data, bits, NULL))				return 0;
	time_e = clock();
	printf("RATS test flag is %d, cost time %.4f.\n", 7, (double)(time_e - time_s) / CLOCKS_PER_SEC);
	//8-二元推导检测
	time_s = clock();
	if (!BinaryDerivative(alpha, data, bits, 3, NULL))		return 0;
	if (!BinaryDerivative(alpha, data, bits, 7, NULL))		return 0;
	time_e = clock();
	printf("RATS test flag is %d, cost time %.4f.\n", 8, (double)(time_e - time_s) / CLOCKS_PER_SEC);
	////9-自相关检测
	time_s = clock();
	if (!Autocorrelation(alpha, data, bits, 1, NULL))		return 0;
	if (!Autocorrelation(alpha, data, bits, 2, NULL))		return 0;
	if (!Autocorrelation(alpha, data, bits, 8, NULL))		return 0;
	if (!Autocorrelation(alpha, data, bits, 16, NULL))		return 0;
	time_e = clock();
	printf("RATS test flag is %d, cost time %.4f.\n", 9, (double)(time_e - time_s) / CLOCKS_PER_SEC);
	//10-矩阵秩检测
	time_s = clock();
	if (!MatrixRank(alpha, data, bits, 32, NULL))			return 0;
	time_e = clock();
	printf("RATS test flag is %d, cost time %.4f.\n", 10, (double)(time_e - time_s) / CLOCKS_PER_SEC);
	//11-累加和检测
	time_s = clock();
	if (!Cumulative(alpha, data, bits, NULL))				return 0;
	time_e = clock();
	printf("RATS test flag is %d, cost time %.4f.\n", 11, (double)(time_e - time_s) / CLOCKS_PER_SEC);
	//12-近似熵检测
	time_s = clock();
	if (!ApproximateEntropy(alpha, data, bits, 10, NULL))	return 0;
	time_e = clock();
	printf("RATS test flag is %d, cost time %.4f.\n", 12, (double)(time_e - time_s) / CLOCKS_PER_SEC);
	//13-线性复杂度检测
	time_s = clock();
	if (!LinearComplexity1(alpha, data, bits, 500, NULL))	return 0;
	time_e = clock();
	printf("RATS test flag is %d, cost time %.4f.\n", 13, (double)(time_e - time_s) / CLOCKS_PER_SEC);
	//14-Maurer通用统计检测
	time_s = clock();
	if (!MaurerUniversal(alpha, data, bits, NULL))			return 0;
	time_e = clock();
	printf("RATS test flag is %d, cost time %.4f.\n", 14, (double)(time_e - time_s) / CLOCKS_PER_SEC);
	//15-离散傅里叶检测
	time_s = clock();
	if (!DiscreteFourierTransform(alpha, data, bits, NULL))	return 0;
	time_e = clock();
	printf("RATS test flag is %d, cost time %.4f.\n", 15, (double)(time_e - time_s) / CLOCKS_PER_SEC);
	//16-离散傅里叶检测
	time_s = clock();
	if (!NonOverlappingTemplateMatchings(alpha, data, bits, 9))	return 0;
	time_e = clock();
	printf("RATS test flag is %d, cost time %.4f.\n", 16, (double)(time_e - time_s) / CLOCKS_PER_SEC);
	//17-离散傅里叶检测
	time_s = clock();
	if (!OverlappingTemplateMatchings(alpha, data, bits, 9))	return 0;
	time_e = clock();
	printf("RATS test flag is %d, cost time %.4f.\n", 17, (double)(time_e - time_s) / CLOCKS_PER_SEC);
	//18-离散傅里叶检测
	time_s = clock();
	if (!RandomExcursions(alpha, data, bits))	return 0;
	time_e = clock();
	printf("RATS test flag is %d, cost time %.4f.\n", 18, (double)(time_e - time_s) / CLOCKS_PER_SEC);

	//19-离散傅里叶检测
	time_s = clock();
	if (!RandomExcursionsVariant(alpha, data, bits))	return 0;
	time_e = clock();
	printf("RATS test flag is %d, cost time %.4f.\n", 19, (double)(time_e - time_s) / CLOCKS_PER_SEC);

	return 1;
}

int main()
{
	unsigned char data[BITS / 8];
	unsigned char temp[1024];
	int b, i, rlen;
	FILE *fp;
	int flag;
	clock_t time_s, time_e;
//printf("1111111111111\n");
	if (!(fp = fopen("../data/data.pi", "rb")))
	{
		printf("File open error!\n");
		goto END;
	}
	//printf("1111111111111\n");
	memset(data, 0, BITS / 8);
	b = 0;
	while (b < BITS)
	{
		rlen = fread(temp, 1, 1024, fp);
		for (i = 0; i < rlen && b < BITS; i++)
		{
			if (temp[i] == '0' || temp[i] == '1')
			{
				data[b / 8] = (data[b / 8] << 1) + (temp[i] - '0');
				b++;
			}
		}
	}
	#define ROUND	10
	//printf("1111111111111\n");
	time_s = clock();
	for (int i = 0; i < ROUND; i++)
	{
		flag = rats_app_sts(data, BITS);
		/* code */
	}
	

	time_e = clock();
	printf("RATS start test flag is %d, cost time %.4f.\n", flag, (double)(time_e - time_s) / CLOCKS_PER_SEC /ROUND );
	//printf("RATS start test flag is %d, cost time %.4f.\n", flag, (double)(time_e - time_s) / CLOCKS_PER_SEC  );
// #define ROUND	4
// 	time_s = clock();
// 	for (i = 0; i < ROUND; i++)
// 	{
// 		flag = rats_app_start(data, 1000000);
// 	}
// 	time_e = clock();
// 	printf("RATS start test flag is %d, cost time %.4f.\n", flag, (double)(time_e - time_s) / CLOCKS_PER_SEC / ROUND);
// #undef	ROUND

// #define ROUND	1000
// 	time_s = clock();
// 	for (i = 0; i < ROUND; i++)
// 	{
// 		flag = rats_app_timer(data, 20000);
// 	}
// 	time_e = clock();
// 	printf("RATS timer test flag is %d, cost time %.4f.\n", flag, (double)(time_e - time_s) / CLOCKS_PER_SEC / ROUND);
// #undef	ROUND

// #define ROUND	200000
// 	time_s = clock();
// 	for (i = 0; i < ROUND; i++)
// 	{
// 		flag = rats_app_using(data, 128);
// 	}
// 	time_e = clock();
// 	printf("RATS using test flag is %d, cost time %.4f.\n", flag, (double)(time_e - time_s) / CLOCKS_PER_SEC / ROUND);
// #undef	ROUND

	fclose(fp);
 END:
//	system("pause");
	return 0;	
}