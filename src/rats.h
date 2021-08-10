#ifndef _RATS_H_
#define _RATS_H_

#define MAX(x,y)             ((x) <  (y)  ? (y)  : (x))
#define MIN(x,y)             ((x) >  (y)  ? (y)  : (x))
#define isNonPositive(x)     ((x) <= 0.e0 ?   1  : 0)
#define isPositive(x)        ((x) >  0.e0 ?   1 : 0)
#define isNegative(x)        ((x) <  0.e0 ?   1 : 0)
#define isGreaterThanOne(x)  ((x) >  1.e0 ?   1 : 0)
#define isZero(x)            ((x) == 0.e0 ?   1 : 0)
#define isOne(x)             ((x) == 1.e0 ?   1 : 0)

#ifdef __cplusplus
extern "C"{	/* start of __cplusplus */
#endif

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
					ITEM-1: MonobitFrequency Test
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

typedef struct _st_value_MonobitFrequency
{
	double p_value;
	double s_obs;
	int sum;
} MonobitFrequency_V;

int MonobitFrequency(double alpha, unsigned char *data, int bits, MonobitFrequency_V *value);

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
					ITEM-2: BlockFrequency Test
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
typedef struct _st_value_BlocktFrequency
{
	double p_value;
	double v_obs;
} BlocktFrequency_V;

int BlockFrequency(double alpha, unsigned char *data, int bits, int M, BlocktFrequency_V *value);

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
					ITEM-3: Poker Test
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
typedef struct _st_value_Poker
{
	double p_value;
	double v_obs;
	unsigned int N;
} Poker_V;
int Poker(double alpha, unsigned char *data, int bits, int m, Poker_V *value);

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
					ITEM-4: Serial Test
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

typedef struct _st_value_Serial
{
	double p_value1;
	double p_value2;
	double del1;
	double del2;
} Serial_V;

int Serial(double alpha, unsigned char *data, int bits, int m, Serial_V *value);

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
					ITEM-5: Runs Test
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
typedef struct _st_value_Runs
{
	double p_value;
	double pi;
	int runs;
} Runs_V;

int Runs(double alpha, unsigned char *data, int bits, Runs_V *value);

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
				ITEM-6: R U N S  D I S T R I B U T I O N  T E S T
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

typedef struct _st_value_RunsDistribution
{
	double p_value;
	double v_obs;
	unsigned int k;							//�γ���ȡֵ: 1-k
	double e[32];							//��[1]��ʼ�洢
	unsigned int Runs0[1024];				//��[1]��ʼ�洢
	unsigned int Runs1[1024];				//��[1]��ʼ�洢
} RunsDistribution_V;

int RunsDistribution(double alpha, unsigned char *data, int bits, RunsDistribution_V *value);


/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
				ITEM-7: L O N G E S T  R U N  T E S T
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
typedef struct _st_value_LongstRun
{
	double p_value;
	double s_obs;
	int sum;
} LongstRun_V;

int LongestRun(double alpha, unsigned char *data, int bits, LongstRun_V *value);

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
				ITEM-8: B I N A R Y  D E R I V A T I V E   T E S T
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

typedef struct _st_value_BinaryDerivative
{
	double p_value;
	double v_obs;
	int sum;
} BinaryDerivative_V;

int BinaryDerivative(double alpha, unsigned char *data, int bits, int k, BinaryDerivative_V *value);

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
				ITEM-9: A U T O C O R R E L A T I O N  T E S T
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

typedef struct _st_value_Autocorrelation
{
	double p_value;
	double v_obs;
	unsigned int A;
} Autocorrelation_V;

int Autocorrelation(double alpha, unsigned char *data, int bits, int d, Autocorrelation_V *value);

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
				ITEM-10: M A T R I X R A N K  T E S T
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

typedef struct _st_value_MatrixRank
{
	double p_value;
	double v_obs;
	double sum;
	double fn;
	double expectedValue;
	double sigma;
} MatrixRank_V;

int MatrixRank(double alpha, unsigned char *data, int bits, int m, MatrixRank_V *value);

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
				ITEM-11: C U M U L A T I V E  T E S T
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

typedef struct _st_value_Cumulative
{
	unsigned int Z;
	double p_value;
} Cumulative_V;

int Cumulative(double alpha, unsigned char *data, int bits, Cumulative_V *value);

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
				ITEM-12: APPROXIMATE ENTROPY TEST
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
typedef struct _st_value_ApproximateEntropy
{
	double p_value;
	double v_obs;
	double sum;
	double fn;
	double expectedValue;
	double sigma;
} ApproximateEntropy_V;

int ApproximateEntropy(double alpha, unsigned char *data, int bits, int m, ApproximateEntropy_V *value);

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
				ITEM-13: LINEAR COMPLEXITY TEST
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
typedef struct _st_value_LinearComplexity
{
	double p_value;
	double v_obs;
	double sum;
	double fn;
	double expectedValue;
	double sigma;
} LinearComplexity_V;

double LinearComplexity(double alpha, unsigned char *data, int bits, int M, double u, int sign, LinearComplexity_V *value);

int LinearComplexity1(double alpha, unsigned char *data, int bits, int M, LinearComplexity_V *value);
int LinearComplexity0(double alpha, unsigned char *data, int bits, int M, LinearComplexity_V *value);

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
				ITEM-14: M A U R E R - U N I V E R S A L - T E S T
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

typedef struct _st_value_MaurerUniversal
{
	double p_value;
	double v_obs;
	double sum;
	double fn;
	double expectedValue;
	double sigma;
} MaurerUniversal_V;

int MaurerUniversal(double alpha, unsigned char *data, int bits, MaurerUniversal_V *value);

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
				ITEM-15: D I S C R E T E  F O U R I E R  T R A N S F O R M  T E S T
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

typedef struct _st_value_DFT
{
	double p_value;
	double v_obs;
	double sum;
	double fn;
	double expectedValue;
	double sigma;
} DFT_V;

int DiscreteFourierTransform(double alpha, unsigned char *data, int bits, DFT_V *value);

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
          N O N O V E R L A P P I N G  T E M P L A T E  T E S T
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
int NonOverlappingTemplateMatchings(double alpha, unsigned char *data, int bits, int m);


/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
               O V E R L A P P I N G  T E M P L A T E  T E S T
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
int OverlappingTemplateMatchings(double alpha, unsigned char *data, int bits, int m);

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
                     R A N D O M  E X C U R S I O N S  T E S T
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

int RandomExcursions(double alpha, unsigned char *data, int bits);

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
            R A N D O M  E X C U R S I O N S  V A R I A N T  T E S T
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

int RandomExcursionsVariant(double alpha, unsigned char *data, int bits);

#ifdef __cplusplus
}			/* end of __cplusplus */
#endif

#endif // !_RATS_H_
