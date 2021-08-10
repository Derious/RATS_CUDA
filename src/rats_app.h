#ifndef _RATS_APP_H_
#define _RATS_APP_H_

#ifdef __cplusplus
extern "C"{	/* start of __cplusplus */
#endif

//随机数开机自检函数（长度1000000-bit，检测成功返回为1，其他失败）
int rats_app_start(unsigned char *data, int bits);

//随机数定期自检函数（长度20000-bit，检测成功返回为1，其他失败）
int rats_app_timer(unsigned char *data, int bits);

//随机数工作自检函数（长度128-bit，检测成功返回为1，其他失败）
int rats_app_using(unsigned char *data, int bits);

#ifdef __cplusplus
}			/* end of __cplusplus */
#endif

#endif // !_RATS_APP_H_
