#include <stdio.h>
#include "rats.h"

int rats_app_start(unsigned char *data, int bits)
{
	double alpha = 0.01;
	int n = 10 * 100000;
	unsigned char err[32] = { 0 };
	int i, j;

	// if (bits < n)
	// 	return 0;

	for (i = 0, j = 0; i < 1; i++)
	{
		//1-Ƶ�����
		if (!MonobitFrequency(alpha, data, bits, NULL))
			err[ 0]++;
		//2-����Ƶ�����
		if (!BlockFrequency(alpha, data, bits, 100, NULL))		
			err[ 1]++;
		//3-�˿˼��
		if (!Poker(alpha, data, bits, 4, NULL))				
			err[ 2]++;
		if (!Poker(alpha, data, bits, 8, NULL))				
			err[ 3]++;
		//4-���м��
		if (!Serial(alpha, data, bits, 2, NULL))				
			err[ 4]++;
		if (!Serial(alpha, data, bits, 5, NULL))				
			err[ 5]++;
		//5-�γ��������
		if (!Runs(alpha, data, bits, NULL))					
			err[ 6]++;
		//6-�γ̷ֲ����
		if (!RunsDistribution(alpha, data, bits, NULL))			
			err[ 7]++;
		//7-�������1�γ̼��
		if (!LongestRun(alpha, data, bits, NULL))					
			err[ 8]++;
		//8-��Ԫ�Ƶ����
		if (!BinaryDerivative(alpha, data, bits, 3, NULL))		
			err[ 9]++;
		if (!BinaryDerivative(alpha, data, bits, 7, NULL))			
			err[10]++;
		//9-����ؼ��
		if (!Autocorrelation(alpha, data, bits, 1, NULL))			
			err[11]++;
		if (!Autocorrelation(alpha, data, bits, 2, NULL))			
			err[12]++;
		if (!Autocorrelation(alpha, data, bits, 8, NULL))			
			err[13]++;
		if (!Autocorrelation(alpha, data, bits, 16, NULL))			
			err[14]++;
		//10-�����ȼ��
		if (!MatrixRank(alpha, data, bits, 32, NULL))				
			err[15]++;
		//11-�ۼӺͼ��
		if (!Cumulative(alpha, data, bits, NULL))					
			err[16]++;
		//12-�����ؼ��
		if (!ApproximateEntropy(alpha, data, bits, 5, NULL))		
			err[17]++;
		//13-���Ը��Ӷȼ��
		if (!LinearComplexity1(alpha, data, bits, 500, NULL))		
			err[18]++;
		//14-Maurerͨ��ͳ�Ƽ��
		if (!MaurerUniversal(alpha, data, bits, NULL))			
			err[19]++;
		//15-��ɢ����Ҷ���
		if (!DiscreteFourierTransform(alpha, data, bits, NULL))	
			err[20]++;

		data += bits / 8;
	}

	for (i = 0; i < 21; i++)
	{
		if (err[i] > 1)
			return -(i + 1);
	}

	return 1;
}

int rats_app_timer(unsigned char *data, int bits)
{
	double alpha = 0.01;
	int n = 10 * 20000;
	unsigned char err[32] = { 0 };
	int i;

	// if (bits < n)
	// 	return 0;

	for (i = 0; i < 10; i++)
	{
		//1-Ƶ�����
		if (!MonobitFrequency(alpha, data, bits, NULL))			
			err[0]++;
		//2-����Ƶ�����
		if (!BlockFrequency(alpha, data, bits, 100, NULL))			
			err[1]++;
		//3-�˿˼��
		if (!Poker(alpha, data, bits, 4, NULL))					
			err[2]++;
		if (!Poker(alpha, data, bits, 8, NULL))					
			err[3]++;
		//4-���м��
		if (!Serial(alpha, data, bits, 2, NULL))					
			err[4]++;
		if (!Serial(alpha, data, bits, 5, NULL))					
			err[5]++;
		//5-�γ��������
		if (!Runs(alpha, data, bits, NULL))						
			err[6]++;
		//6-�γ̷ֲ����
		if (!RunsDistribution(alpha, data, bits, NULL))			
			err[7]++;
		//7-�������1�γ̼��
		if (!LongestRun(alpha, data, bits, NULL))					
			err[8]++;
		//8-��Ԫ�Ƶ����
		if (!BinaryDerivative(alpha, data, bits, 3, NULL))			
			err[9]++;
		if (!BinaryDerivative(alpha, data, bits, 7, NULL))			
			err[10]++;
		//9-����ؼ��
		if (!Autocorrelation(alpha, data, bits, 1, NULL))			
			err[11]++;
		if (!Autocorrelation(alpha, data, bits, 2, NULL))			
			err[12]++;
		if (!Autocorrelation(alpha, data, bits, 8, NULL))			
			err[13]++;
		if (!Autocorrelation(alpha, data, bits, 16, NULL))			
			err[14]++;
		//10-�����ȼ��
		if (!MatrixRank(alpha, data, bits, 32, NULL))				
			err[15]++;
		//11-�ۼӺͼ��
		if (!Cumulative(alpha, data, bits, NULL))					
			err[16]++;
		//12-�����ؼ��
		if (!ApproximateEntropy(alpha, data, bits, 5, NULL))		
			err[17]++;
		////13-���Ը��Ӷȼ��
		//if (!LinearComplexity(alpha, data, 20000, 500, NULL))		err[18]++;
		////14-Maurerͨ��ͳ�Ƽ��
		//if (!MaurerUniversal(alpha, data, 20000, NULL))				err[19]++;
		////15-��ɢ����Ҷ���
		//if (!DiscreteFourierTransform(alpha, data, 20000, NULL))	err[20]++;

		data += 20000 / 8;
	}

	for (i = 0; i < 21; i++)
	{
		if (err[i] > 1)
			return -(i + 1);
	}

	return 1;
}

int rats_app_using(unsigned char *data, int bits)
{
	double alpha = 0.01;
	int ret = 0;

	// if (bits < 128)
	// 	return 0;
	// bits = 128;

	//1-Ƶ�����
	if (!MonobitFrequency(alpha, data, bits, NULL))			ret |= 0x00000001;

	//3-�˿˼��
	if (!Poker(alpha, data, bits, 2, NULL))					ret |= 0x00000004;

	if (ret)
		ret |= 0x80000000;
	else
		ret = 1;

	return ret;
}