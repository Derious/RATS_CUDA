#ifndef _RATS_APP_H_
#define _RATS_APP_H_

#ifdef __cplusplus
extern "C"{	/* start of __cplusplus */
#endif

//����������Լ캯��������1000000-bit�����ɹ�����Ϊ1������ʧ�ܣ�
int rats_app_start(unsigned char *data, int bits);

//����������Լ캯��������20000-bit�����ɹ�����Ϊ1������ʧ�ܣ�
int rats_app_timer(unsigned char *data, int bits);

//����������Լ캯��������128-bit�����ɹ�����Ϊ1������ʧ�ܣ�
int rats_app_using(unsigned char *data, int bits);

#ifdef __cplusplus
}			/* end of __cplusplus */
#endif

#endif // !_RATS_APP_H_
