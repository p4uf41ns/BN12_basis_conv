#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include <time.h>
#include <sys/time.h>
#include <unistd.h>
#include <string.h>
#define X_length 114
#define X6_2_length 116

/*============================================================================*/
/* Field                                                                      */
/*============================================================================*/
typedef struct{
    mpz_t x0;
}Fp;

typedef struct{
    Fp x0,x1;
}Fp2;

typedef struct{
    Fp2 x0,x1,x2;
}Fp6;

typedef struct{
    Fp6 x0,x1;
}Fp12;
/*----------------------------------------------------------------------------*/
//Fp
void Fp_init(Fp *A);
void Fp_clear(Fp *A);
void Fp_printf(Fp *A,char *str);
void Fp_set(Fp *ANS,Fp *A);
void Fp_set_ui(Fp *ANS,unsigned long int UI);
void Fp_set_mpz(Fp *ANS,mpz_t A);
void Fp_set_neg(Fp *ANS,Fp *A);
void Fp_set_random(Fp *ANS,gmp_randstate_t state);
void Fp_mul(Fp *ANS,Fp *A,Fp *B);
void Fp_mul_ui(Fp *ANS,Fp *A,unsigned long int UI);
void Fp_mul_mpz(Fp *ANS,Fp *A,mpz_t B);
void Fp_add(Fp *ANS,Fp *A,Fp *B);
void Fp_add_ui(Fp *ANS,Fp *A,unsigned long int UI);
void Fp_add_mpz(Fp *ANS,Fp *A,mpz_t B);
void Fp_sub(Fp *ANS,Fp *A,Fp *B);
void Fp_sub_ui(Fp *ANS,Fp *A,unsigned long int UI);
void Fp_sub_mpz(Fp *ANS,Fp *A,mpz_t B);
void Fp_inv(Fp *ANS,Fp *A);
int  Fp_legendre(Fp *A);
int  Fp_isCNR(Fp *A);
void Fp_sqrt(Fp *ANS,Fp *A);
void Fp_pow(Fp *ANS,Fp *A,mpz_t scalar);
int  Fp_cmp(Fp *A,Fp *B);
int  Fp_cmp_ui(Fp *A,unsigned long int UI);
int  Fp_cmp_mpz(Fp *A,mpz_t B);
int  Fp_cmp_zero(Fp *A);
int  Fp_cmp_one(Fp *A);
/*----------------------------------------------------------------------------*/
//Fp2
void Fp2_init(Fp2 *A);
void Fp2_clear(Fp2 *A);
void Fp2_printf(Fp2 *A,char *str);
void Fp2_set(Fp2 *ANS,Fp2 *A);
void Fp2_set_ui(Fp2 *ANS,unsigned long int UI);
void Fp2_set_mpz(Fp2 *ANS,mpz_t A);
void Fp2_set_neg(Fp2 *ANS,Fp2 *A);
void Fp2_set_random(Fp2 *ANS,gmp_randstate_t state);
void Fp2_mul(Fp2 *ANS,Fp2 *A,Fp2 *B);
void Fp2_mul_ui(Fp2 *ANS,Fp2 *A,unsigned long int UI);
void Fp2_mul_mpz(Fp2 *ANS,Fp2 *A,mpz_t B);
void Fp2_mul_basis(Fp2 *ANS,Fp2 *A);
void Fp2_inv_basis(Fp2 *ANS,Fp2 *A);
void Fp2_sqr(Fp2 *ANS,Fp2 *A);
void Fp2_add(Fp2 *ANS,Fp2 *A,Fp2 *B);
void Fp2_add_ui(Fp2 *ANS,Fp2 *A,unsigned long int UI);
void Fp2_add_mpz(Fp2 *ANS,Fp2 *A,mpz_t B);
void Fp2_sub(Fp2 *ANS,Fp2 *A,Fp2 *B);
void Fp2_sub_ui(Fp2 *ANS,Fp2 *A,unsigned long int UI);
void Fp2_sub_mpz(Fp2 *ANS,Fp2 *A,mpz_t B);
void Fp2_inv(Fp2 *ANS,Fp2 *A);
int  Fp2_legendre(Fp2 *A);
int  Fp2_isCNR(Fp2 *A);
void Fp2_sqrt(Fp2 *ANS,Fp2 *A);
void Fp2_pow(Fp2 *ANS,Fp2 *A,mpz_t scalar);
int  Fp2_cmp(Fp2 *A,Fp2 *B);
int  Fp2_cmp_ui(Fp2 *A,unsigned long int UI);
int  Fp2_cmp_mpz(Fp2 *A,mpz_t B);
int  Fp2_cmp_zero(Fp2 *A);
int  Fp2_cmp_one(Fp2 *A);
/*----------------------------------------------------------------------------*/
//Fp6
void Fp6_init(Fp6 *A);
void Fp6_clear(Fp6 *A);
void Fp6_printf(Fp6 *A,char *str);
void Fp6_set(Fp6 *ANS,Fp6 *A);
void Fp6_set_ui(Fp6 *ANS,unsigned long int UI);
void Fp6_set_mpz(Fp6 *ANS,mpz_t A);
void Fp6_set_neg(Fp6 *ANS,Fp6 *A);
void Fp6_set_random(Fp6 *ANS,gmp_randstate_t state);
void Fp6_mul_ui(Fp6 *ANS,Fp6 *A,unsigned long int UI);
void Fp6_mul_mpz(Fp6 *ANS,Fp6 *A,mpz_t B);
void Fp6_mul_basis(Fp6 *ANS,Fp6 *A);
void Fp6_add(Fp6 *ANS,Fp6 *A,Fp6 *B);
void Fp6_add_ui(Fp6 *ANS,Fp6 *A,unsigned long int UI);
void Fp6_add_mpz(Fp6 *ANS,Fp6 *A,mpz_t B);
void Fp6_sub(Fp6 *ANS,Fp6 *A,Fp6 *B);
void Fp6_sub_ui(Fp6 *ANS,Fp6 *A,unsigned long int UI);
void Fp6_sub_mpz(Fp6 *ANS,Fp6 *A,mpz_t B);
int  Fp6_cmp(Fp6 *A,Fp6 *B);
int  Fp6_cmp_ui(Fp6 *A,unsigned long int UI);
int  Fp6_cmp_mpz(Fp6 *A,mpz_t B);
int  Fp6_cmp_zero(Fp6 *A);
//karatsuba (T'_1)
void Fp6_mul_karatsuba(Fp6 *ANS,Fp6 *A,Fp6 *B);
void Fp6_sqr_karatsuba(Fp6 *ANS,Fp6 *A);
void Fp6_inv_karatsuba(Fp6 *ANS,Fp6 *A);
int  Fp6_legendre_karatsuba(Fp6 *A);
int  Fp6_isCNR_karatsuba(Fp6 *A);
void Fp6_sqrt_karatsuba(Fp6 *ANS,Fp6 *A);
void Fp6_pow_karatsuba(Fp6 *ANS,Fp6 *A,mpz_t scalar);
int  Fp6_cmp_one_karatsuba(Fp6 *A);
//CVMA (T'_2)
void Fp6_mul_CVMA(Fp6 *ANS,Fp6 *A,Fp6 *B);
void Fp6_sqr_CVMA(Fp6 *ANS,Fp6 *A);
void Fp6_inv_CVMA(Fp6 *ANS,Fp6 *A);
int  Fp6_legendre_CVMA(Fp6 *A);
int  Fp6_isCNR_CVMA(Fp6 *A);
void Fp6_sqrt_CVMA(Fp6 *ANS,Fp6 *A);
void Fp6_pow_CVMA(Fp6 *ANS,Fp6 *A,mpz_t scalar);
int  Fp6_cmp_one_CVMA(Fp6 *A);
void Fp6_frobenius_map_p1_CVMA(Fp6 *ANS,Fp6 *A);
void Fp6_frobenius_map_p2_CVMA(Fp6 *ANS,Fp6 *A);
void Fp6_frobenius_map_p3_CVMA(Fp6 *ANS,Fp6 *A);
void Fp6_frobenius_map_p4_CVMA(Fp6 *ANS,Fp6 *A);
void Fp6_frobenius_map_p5_CVMA(Fp6 *ANS,Fp6 *A);
/*----------------------------------------------------------------------------*/
//Fp12
void Fp12_init(Fp12 *A);
void Fp12_clear(Fp12 *A);
void Fp12_printf(Fp12 *A,char *str);
void Fp12_set(Fp12 *ANS,Fp12 *A);
void Fp12_set_ui(Fp12 *ANS,unsigned long int UI);
void Fp12_set_mpz(Fp12 *ANS,mpz_t A);
void Fp12_set_neg(Fp12 *ANS,Fp12 *A);
void Fp12_set_random(Fp12 *ANS,gmp_randstate_t state);
void Fp12_mul_ui(Fp12 *ANS,Fp12 *A,unsigned long int UI);
void Fp12_mul_mpz(Fp12 *ANS,Fp12 *A,mpz_t B);
void Fp12_add(Fp12 *ANS,Fp12 *A,Fp12 *B);
void Fp12_add_ui(Fp12 *ANS,Fp12 *A,unsigned long int UI);
void Fp12_add_mpz(Fp12 *ANS,Fp12 *A,mpz_t B);
void Fp12_sub(Fp12 *ANS,Fp12 *A,Fp12 *B);
void Fp12_sub_ui(Fp12 *ANS,Fp12 *A,unsigned long int UI);
void Fp12_sub_mpz(Fp12 *ANS,Fp12 *A,mpz_t B);
int  Fp12_cmp(Fp12 *A,Fp12 *B);
int  Fp12_cmp_ui(Fp12 *A,unsigned long int UI);
int  Fp12_cmp_mpz(Fp12 *A,mpz_t B);
int  Fp12_cmp_zero(Fp12 *A);
//karatsuba (T'_1)
void Fp12_mul_karatsuba(Fp12 *ANS,Fp12 *A,Fp12 *B);
void Fp12_sqr_karatsuba(Fp12 *ANS,Fp12 *A);
void Fp12_inv_karatsuba(Fp12 *ANS,Fp12 *A);
int  Fp12_legendre_karatsuba(Fp12 *A);
void Fp12_sqrt_karatsuba(Fp12 *ANS,Fp12 *A);
void Fp12_pow_karatsuba(Fp12 *ANS,Fp12 *A,mpz_t scalar);
int  Fp12_cmp_one_karatsuba(Fp12 *A);
void Fp12_frobenius_map_p1_karatsuba(Fp12 *ANS,Fp12 *A);
void Fp12_frobenius_map_p2_karatsuba(Fp12 *ANS,Fp12 *A);
void Fp12_frobenius_map_p3_karatsuba(Fp12 *ANS,Fp12 *A);
void Fp12_frobenius_map_p4_karatsuba(Fp12 *ANS,Fp12 *A);
void Fp12_frobenius_map_p6_karatsuba(Fp12 *ANS,Fp12 *A);
void Fp12_frobenius_map_p8_karatsuba(Fp12 *ANS,Fp12 *A);
void Fp12_frobenius_map_p10_karatsuba(Fp12 *ANS,Fp12 *A);
//CVMA (T'_2)
void Fp12_mul_CVMA(Fp12 *ANS,Fp12 *A,Fp12 *B);
void Fp12_sqr_CVMA(Fp12 *ANS,Fp12 *A);
void Fp12_sqr_ctclotomic_CVMA(Fp12 *ANS,Fp12 *A);
void Fp12_inv_CVMA(Fp12 *ANS,Fp12 *A);
int  Fp12_legendre_CVMA(Fp12 *A);
void Fp12_sqrt_CVMA(Fp12 *ANS,Fp12 *A);
void Fp12_pow_CVMA(Fp12 *ANS,Fp12 *A,mpz_t scalar);
int  Fp12_cmp_one_CVMA(Fp12 *A);
void Fp12_frobenius_map_p1_CVMA(Fp12 *ANS,Fp12 *A);
void Fp12_frobenius_map_p2_CVMA(Fp12 *ANS,Fp12 *A);
void Fp12_frobenius_map_p3_CVMA(Fp12 *ANS,Fp12 *A);
void Fp12_frobenius_map_p4_CVMA(Fp12 *ANS,Fp12 *A);
void Fp12_frobenius_map_p6_CVMA(Fp12 *ANS,Fp12 *A);
void Fp12_frobenius_map_p8_CVMA(Fp12 *ANS,Fp12 *A);
void Fp12_frobenius_map_p10_CVMA(Fp12 *ANS,Fp12 *A);

/*============================================================================*/
/* Elliptic Curve                                                             */
/*============================================================================*/
typedef struct{
	Fp x,y;
	int infinity;
}EFp;

typedef struct{
    Fp2 x,y;
	int infinity;
}EFp2;

typedef struct{
    Fp6 x,y;
	int infinity;
}EFp6;

typedef struct{
    Fp12 x,y;
	int infinity;
}EFp12;
/*----------------------------------------------------------------------------*/
//EFp
void EFp_init(EFp *P);
void EFp_set(EFp *P,EFp *A);
void EFp_set_ui(EFp *ANS,unsigned long int UI);
void EFp_set_mpz(EFp *ANS,mpz_t A);
void EFp_set_neg(EFp *ANS,EFp *A);
void EFp_clear(EFp *P);
void EFp_printf(EFp *P,char *str);
void EFp_rational_point(EFp *P);
void EFp_ECD(EFp *ANS,EFp *P);
void EFp_ECA(EFp *ANS,EFp *P1,EFp *P2);
void EFp_SCM(EFp *ANS,EFp *P,mpz_t scalar);
void EFp_skew_frobenius_map_p2(EFp *ANS,EFp *A);
/*----------------------------------------------------------------------------*/
//EFp2
void EFp2_init(EFp2 *P);
void EFp2_set(EFp2 *ANS,EFp2 *A);
void EFp2_set_ui(EFp2 *ANS,unsigned long int UI);
void EFp2_set_mpz(EFp2 *ANS,mpz_t A);
void EFp2_set_neg(EFp2 *ANS,EFp2 *A);
void EFp2_clear(EFp2 *P);
void EFp2_printf(EFp2 *P,char *str);
void EFp2_rational_point(EFp2 *P);
void EFp2_ECD(EFp2 *ANS,EFp2 *P);
void EFp2_ECA(EFp2 *ANS,EFp2 *P1,EFp2 *P2);
void EFp2_SCM(EFp2 *ANS,EFp2 *P,mpz_t scalar);
void EFp2_skew_frobenius_map_p1_karatsuba(EFp2 *ANS,EFp2 *A);
void EFp2_skew_frobenius_map_p2_karatsuba(EFp2 *ANS,EFp2 *A);
void EFp2_skew_frobenius_map_p3_karatsuba(EFp2 *ANS,EFp2 *A);
void EFp2_skew_frobenius_map_p10_karatsuba(EFp2 *ANS,EFp2 *A);
/*----------------------------------------------------------------------------*/
//EFp6
void EFp6_init(EFp6 *P);
void EFp6_set(EFp6 *ANS,EFp6 *A);
void EFp6_set_ui(EFp6 *ANS,unsigned long int UI);
void EFp6_set_mpz(EFp6 *ANS,mpz_t A);
void EFp6_set_neg(EFp6 *ANS,EFp6 *A);
void EFp6_clear(EFp6 *P);
void EFp6_printf(EFp6 *P,char *str);
//karatsuba (T'_1)
void EFp6_rational_point_karatsuba(EFp6 *P);
void EFp6_ECD_karatsuba(EFp6 *ANS,EFp6 *P);
void EFp6_ECA_karatsuba(EFp6 *ANS,EFp6 *P1,EFp6 *P2);
void EFp6_SCM_karatsuba(EFp6 *ANS,EFp6 *P,mpz_t scalar);
//CVMA (T'_2)
void EFp6_rational_point_CVMA(EFp6 *P);
void EFp6_ECD_CVMA(EFp6 *ANS,EFp6 *P);
void EFp6_ECA_CVMA(EFp6 *ANS,EFp6 *P1,EFp6 *P2);
void EFp6_SCM_CVMA(EFp6 *ANS,EFp6 *P,mpz_t scalar);
/*----------------------------------------------------------------------------*/
//EFp12
void EFp12_init(EFp12 *P);
void EFp12_set(EFp12 *ANS,EFp12 *A);
void EFp12_set_ui(EFp12 *ANS,unsigned long int UI);
void EFp12_set_mpz(EFp12 *ANS,mpz_t A);
void EFp12_set_neg(EFp12 *ANS,EFp12 *A);
void EFp12_clear(EFp12 *P);
void EFp12_printf(EFp12 *P,char *str);
//karatsuba (T'_1)
void EFp12_rational_point_karatsuba(EFp12 *P);
void EFp12_generate_G1_karatsuba(EFp12 *P);
void EFp12_generate_G2_karatsuba(EFp12 *Q);
void EFp12_ECD_karatsuba(EFp12 *ANS,EFp12 *P);
void EFp12_ECA_karatsuba(EFp12 *ANS,EFp12 *P1,EFp12 *P2);
void EFp12_SCM_karatsuba(EFp12 *ANS,EFp12 *P,mpz_t scalar);
//CVMA (T'_2)
void EFp12_rational_point_CVMA(EFp12 *P);
void EFp12_generate_G1_CVMA(EFp12 *P);
void EFp12_generate_G2_CVMA(EFp12 *Q);
void EFp12_ECD_CVMA(EFp12 *ANS,EFp12 *P);
void EFp12_ECA_CVMA(EFp12 *ANS,EFp12 *P1,EFp12 *P2);
void EFp12_SCM_CVMA(EFp12 *ANS,EFp12 *P,mpz_t scalar);

/*============================================================================*/
/* Pairing                                                                    */
/*============================================================================*/
enum f_state{
	f_p1,f_p2,f_p3,f_p4,f_p5,f_p6,f_p7,f_p8,f_p9,f_p10,f_p11,f_p12
};
int X_binary[X_length+1];
int X6_2_binary[X6_2_length+1];
mpz_t X,prime,order,trace;
mpz_t EFp_total,EFp12_total;
mpz_t curve_b;
Fp2 frobenius_constant_karatsuba[12][6];
Fp2 frobenius_constant_CVMA[12];
Fp2 skew_frobenius_constant_karatsuba[12][2];
Fp2 Alpha_1,Alpha_1_inv;
mpz_t epsilon1,epsilon2;
mpz_t Two_inv;

typedef struct{
	Fp data[3][3];
}MATRIX;
MATRIX Matrix_karatsuba_to_CVMA;
MATRIX Matrix_CVMA_to_karatsuba;

Fp TMP1_FP,TMP2_FP,TMP3_FP,TMP4_FP,TMP5_FP;
Fp2 TMP1_FP2,TMP2_FP2,TMP3_FP2,TMP4_FP2,TMP5_FP2,TMP6_FP2,TMP7_FP2,TMP8_FP2,TMP9_FP2,TMP10_FP2,TMP11_FP2,TMP12_FP2,TMP13_FP2,TMP14_FP2;
Fp6 TMP1_FP6,TMP2_FP6,TMP3_FP6,TMP4_FP6;
Fp12 TMP1_FP12,TMP2_FP12,TMP3_FP12,TMP4_FP12,TMP5_FP12;
EFp TMP1_EFP,TMP2_EFP;
EFp2 TMP1_EFP2,TMP2_EFP2;
EFp6 TMP1_EFP6,TMP2_EFP6;
EFp12 TMP1_EFP12,TMP2_EFP12;

/*----------------------------------------------------------------------------*/
//Basis conversion
void karatsuba_to_CVMA(Fp12 *ANS,Fp12 *A);
void CVMA_to_karatsuba(Fp12 *ANS,Fp12 *A);
/*----------------------------------------------------------------------------*/
//twist
void EFp12_to_EFp2_karatsuba(EFp2 *ANS,EFp12 *A);
void EFp2_to_EFp12_karatsuba(EFp12 *ANS,EFp2 *A);
void EFp12_to_EFp_karatsuba(EFp *ANS,EFp12 *A);
void EFp_to_EFp12_karatsuba(EFp12 *ANS,EFp *A);
void EFp12_to_EFp_CVMA(EFp *ANS,EFp12 *A);
void EFp_to_EFp12_CVMA(EFp12 *ANS,EFp *A);
/*----------------------------------------------------------------------------*/
//Pseudo 8-sparse
void Pseudo_8_sparse_mapping(EFp *P,EFp2 *Q,Fp *L);
void Pseudo_8_sparse_mul_karatsuba(Fp12 *ANS,Fp12 *A,Fp12 *B);
void ff_ltt_karatsuba(Fp12 *f,EFp2 *T,EFp *P,Fp *L);
void f_ltq_karatsuba(Fp12 *f,EFp2 *T,EFp2 *Q,EFp *P,Fp *L);
/*----------------------------------------------------------------------------*/
//miller
void Miller_algo_for_plain_ate_karatsuba(Fp12 *ANS,EFp12 *Q,EFp12 *P);
void Miller_algo_for_opt_ate_karatsuba(Fp12 *ANS,EFp12 *Q,EFp12 *P);
void Miller_algo_for_x_ate_karatsuba(Fp12 *ANS,EFp12 *Q,EFp12 *P);
/*----------------------------------------------------------------------------*/
//final exp
void Fp12_pow_X_CVMA(Fp12 *ANS,Fp12 *A);
void Final_exp_easy_CVMA(Fp12 *ANS,Fp12 *f);
void Final_exp_hard_plain_CVMA(Fp12 *ANS,Fp12 *f);
void Final_exp_hard_optimal_CVMA(Fp12 *ANS,Fp12 *f);
/*----------------------------------------------------------------------------*/
//pairing
void Plain_ate_pairing(Fp12 *ANS,EFp12 *P,EFp12 *Q);
void Opt_ate_pairing(Fp12 *ANS,EFp12 *P,EFp12 *Q);
void X_ate_pairing(Fp12 *ANS,EFp12 *P,EFp12 *Q);
/*----------------------------------------------------------------------------*/
//JSF
void Joint_sparse_form(int **binary,mpz_t S[2],int *loop_length);
//G1 SCM
void EFp12_G1_SCM_plain_karatsuba(EFp12 *ANS,EFp12 *P,mpz_t scalar);
void EFp12_G1_SCM_2split_karatsuba(EFp12 *ANS,EFp12 *P,mpz_t scalar);
void EFp12_G1_SCM_2split_JSF_karatsuba(EFp12 *ANS,EFp12 *P,mpz_t scalar);
/*----------------------------------------------------------------------------*/
//G2 SCM
void EFp12_G2_SCM_plain_karatsuba(EFp12 *ANS,EFp12 *Q,mpz_t scalar);
void EFp12_G2_SCM_2split_karatsuba(EFp12 *ANS,EFp12 *Q,mpz_t scalar);
void EFp12_G2_SCM_2split_JSF_karatsuba(EFp12 *ANS,EFp12 *Q,mpz_t scalar);
void EFp12_G2_SCM_4split_karatsuba(EFp12 *ANS,EFp12 *Q,mpz_t scalar);
/*----------------------------------------------------------------------------*/
//G3 EXP
void Fp12_G3_EXP_plain_CVMA(Fp12 *ANS,Fp12 *A,mpz_t scalar);
void Fp12_G3_EXP_2split_CVMA(Fp12 *ANS,Fp12 *A,mpz_t scalar);
void Fp12_G3_EXP_2split_JSF_CVMA(Fp12 *ANS,Fp12 *A,mpz_t scalar);
void Fp12_G3_EXP_4split_CVMA(Fp12 *ANS,Fp12 *A,mpz_t scalar);
/*----------------------------------------------------------------------------*/
//Settings
void BN12_init();
void BN12_clear();
void BN12_print_parameters();
void init_parameters();
void generate_X();
int  generate_prime();
int  generate_order();
void generate_trace();
void weil();
void get_epsilon();
void get_Two_inv();
void set_basis();
void set_frobenius_constant_karatsuba();
void set_frobenius_constant_CVMA();
void set_curve_parameter();
void get_basis_conversion_matrix();
/*----------------------------------------------------------------------------*/
//Matrix operations
void Matrix_init(MATRIX *M);
void Matrix_set(MATRIX *ANS,MATRIX *M);
void Matrix_clear(MATRIX *M);
void Matrix_printf(MATRIX *M);
int Matrix_cmp(MATRIX *M,MATRIX *N);
void Matrix_inv(MATRIX *ANS,MATRIX *A);
void Matrix_mul(MATRIX *ANS,MATRIX *A,MATRIX *B);

/*============================================================================*/
/* Test                                                                       */
/*============================================================================*/
struct timeval tv_start,tv_end;
float BASIS_CONVERSION_CTOK,BASIS_CONVERSION_KTOC;
float MILLER_TATE,MILLER_PLAINATE,MILLER_OPTATE,MILLER_XATE;
float FINALEXP_EASY,FINALEXP_HARD_PLAIN,FINALEXP_HARD_OPT;
float G1SCM_PLAIN,G1SCM_2SPLIT,G1SCM_2SPLIT_JSF;
float G2SCM_PLAIN,G2SCM_2SPLIT,G2SCM_2SPLIT_JSF,G2SCM_4SPLIT;
float G3EXP_PLAIN,G3EXP_2SPLIT,G3EXP_2SPLIT_JSF,G3EXP_4SPLIT;
struct mpz_Cost{
    unsigned long int mpz_mul;
    unsigned long int mpz_mul_ui;
    unsigned long int mpz_sqr;
    unsigned long int mpz_add;
    unsigned long int mpz_add_ui;
    unsigned long int mpz_invert;
};
struct Fp_Cost{
    unsigned long int Fp_mul;
    unsigned long int Fp_mul_mpz;
    unsigned long int Fp_mul_ui;
    unsigned long int Fp_sqr;
    unsigned long int Fp_basis;
    unsigned long int Fp_add;
    unsigned long int Fp_add_mpz;
    unsigned long int Fp_add_ui;
    unsigned long int Fp_inv;
    unsigned long int Fp_neg;
};
struct mpz_Cost mpz_cost;
struct Fp_Cost Fp_cost;

/*----------------------------------------------------------------------------*/
//time
float timedifference_msec(struct timeval t0, struct timeval t1);
float timedifference_usec(struct timeval t0, struct timeval t1);
/*----------------------------------------------------------------------------*/
//cost
void Init_mpz_Cost(struct mpz_Cost *cost);
void Print_mpz_Cost(struct mpz_Cost *cost,char *str);
void Init_Fp_Cost(struct Fp_Cost *cost);
void Print_Fp_Cost(struct Fp_Cost *cost,char *str);
/*----------------------------------------------------------------------------*/
//test
void test_plain_ate_pairing();
void test_opt_ate_pairing();
void test_x_ate_pairing();
void test_G1_SCM();
void test_G2_SCM();
void test_G3_EXP();
void test_basis_conversion();
void computation_time();
void computation_cost();
