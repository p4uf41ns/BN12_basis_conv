#include "bn12.h"

/*============================================================================*/
/* main                                                                       */
/*============================================================================*/
int main(void){
    BN12_init();
    BN12_print_parameters();
    
    //test_plain_ate_pairing();         //Miller's Alg. (T'_1) and Final Exp. (T'_2)
    //test_opt_ate_pairing();           //Miller's Alg. (T'_1) and Final Exp. (T'_2)
    test_x_ate_pairing();             //Miller's Alg. (T'_1) and Final Exp. (T'_2)
    //test_G1_SCM();                    //Scalar multiplication in G1 (T'_1)
    //test_G2_SCM();                    //Scalar multiplication in G2 (T'_1)
    //test_G3_EXP();                    //Exponentiation in G3 (T'_2)
    //test_basis_conversion();          //Basis conversion (T'_1 to T'_2 / T'_2 to T'_1)
    //computation_time();
    //computation_cost();
    
    BN12_clear();
    return 0;
}

/*============================================================================*/
/* Field                                                                      */
/*============================================================================*/
//Fp
void Fp_init(Fp *A){
    mpz_init(A->x0);
}

void Fp_clear(Fp *A){
    mpz_clear(A->x0);
}

void Fp_printf(Fp *A,char *str){
    gmp_printf("%s%Zd",str,A->x0);
}

void Fp_set(Fp *ANS,Fp *A){
    mpz_set(ANS->x0,A->x0);
}

void Fp_set_ui(Fp *ANS,unsigned long int UI){
    mpz_set_ui(ANS->x0,UI);
}

void Fp_set_mpz(Fp *ANS,mpz_t A){
    mpz_set(ANS->x0,A);
}

//---------------------------------------------------------------//
//time

void Fp_set_neg(Fp *ANS,Fp *A){
    mpz_sub(ANS->x0,prime,A->x0);
    mpz_mod(ANS->x0,ANS->x0,prime);
}

void Fp_set_random(Fp *ANS,gmp_randstate_t state){
    mpz_urandomm(ANS->x0,state,prime);
}

void Fp_mul(Fp *ANS,Fp *A,Fp *B){
    mpz_mul(ANS->x0,A->x0,B->x0);
    mpz_mod(ANS->x0,ANS->x0,prime);
}

void Fp_mul_ui(Fp *ANS,Fp *A,unsigned long int UI){
    mpz_mul_ui(ANS->x0,A->x0,UI);
    mpz_mod(ANS->x0,ANS->x0,prime);
}

void Fp_mul_mpz(Fp *ANS,Fp *A,mpz_t B){
    mpz_mul(ANS->x0,A->x0,B);
    mpz_mod(ANS->x0,ANS->x0,prime);
}

void Fp_add(Fp *ANS,Fp *A,Fp *B){
    mpz_add(ANS->x0,A->x0,B->x0);
    mpz_mod(ANS->x0,ANS->x0,prime);
}

void Fp_add_ui(Fp *ANS,Fp *A,unsigned long int UI){
    mpz_add_ui(ANS->x0,A->x0,UI);
    mpz_mod(ANS->x0,ANS->x0,prime);
}

void Fp_add_mpz(Fp *ANS,Fp *A,mpz_t B){
    mpz_add(ANS->x0,A->x0,B);
    mpz_mod(ANS->x0,ANS->x0,prime);
}

void Fp_sub(Fp *ANS,Fp *A,Fp *B){
    mpz_sub(ANS->x0,A->x0,B->x0);
    mpz_mod(ANS->x0,ANS->x0,prime);
}

void Fp_sub_ui(Fp *ANS,Fp *A,unsigned long int UI){
    mpz_sub_ui(ANS->x0,A->x0,UI);
    mpz_mod(ANS->x0,ANS->x0,prime);
}

void Fp_sub_mpz(Fp *ANS,Fp *A,mpz_t B){
    mpz_sub(ANS->x0,A->x0,B);
    mpz_mod(ANS->x0,ANS->x0,prime);
}

void Fp_inv(Fp *ANS,Fp *A){
    mpz_invert(ANS->x0,A->x0,prime);
}

//---------------------------------------------------------------//
//cost
/*
void Fp_set_neg(Fp *ANS,Fp *A){
    mpz_cost.mpz_add++;
    Fp_cost.Fp_neg++;
    mpz_sub(ANS->x0,prime,A->x0);
    mpz_mod(ANS->x0,ANS->x0,prime);
}

void Fp_set_random(Fp *ANS,gmp_randstate_t state){
    mpz_urandomm(ANS->x0,state,prime);
}

void Fp_mul(Fp *ANS,Fp *A,Fp *B){
    if(mpz_cmp(A->x0,B->x0)==0){
        mpz_cost.mpz_sqr++;
        Fp_cost.Fp_sqr++;
    }else{
        mpz_cost.mpz_mul++;
        Fp_cost.Fp_mul++;
    }
    mpz_mul(ANS->x0,A->x0,B->x0);
    mpz_mod(ANS->x0,ANS->x0,prime);
}

void Fp_mul_ui(Fp *ANS,Fp *A,unsigned long int UI){
    Fp_cost.Fp_mul_ui++;
    mpz_cost.mpz_mul_ui++;
    mpz_mul_ui(ANS->x0,A->x0,UI);
    mpz_mod(ANS->x0,ANS->x0,prime);
}

void Fp_mul_mpz(Fp *ANS,Fp *A,mpz_t B){
    if(mpz_cmp(A->x0,B)==0){
        mpz_cost.mpz_sqr++;
        Fp_cost.Fp_sqr++;
    }else{
        mpz_cost.mpz_mul++;
        Fp_cost.Fp_mul_mpz++;
    }
    mpz_mul(ANS->x0,A->x0,B);
    mpz_mod(ANS->x0,ANS->x0,prime);
}

void Fp_add(Fp *ANS,Fp *A,Fp *B){
    Fp_cost.Fp_add++;
    mpz_cost.mpz_add++;
    mpz_add(ANS->x0,A->x0,B->x0);
    mpz_mod(ANS->x0,ANS->x0,prime);
}

void Fp_add_ui(Fp *ANS,Fp *A,unsigned long int UI){
    Fp_cost.Fp_add_ui++;
    mpz_cost.mpz_add_ui++;
    mpz_add_ui(ANS->x0,A->x0,UI);
    mpz_mod(ANS->x0,ANS->x0,prime);
}

void Fp_add_mpz(Fp *ANS,Fp *A,mpz_t B){
    Fp_cost.Fp_add_mpz++;
    mpz_cost.mpz_add++;
    mpz_add(ANS->x0,A->x0,B);
    mpz_mod(ANS->x0,ANS->x0,prime);
}

void Fp_sub(Fp *ANS,Fp *A,Fp *B){
    Fp_cost.Fp_add++;
    mpz_cost.mpz_add++;
    mpz_sub(ANS->x0,A->x0,B->x0);
    mpz_mod(ANS->x0,ANS->x0,prime);
}

void Fp_sub_ui(Fp *ANS,Fp *A,unsigned long int UI){
    Fp_cost.Fp_add_ui++;
    mpz_cost.mpz_add_ui++;
    mpz_sub_ui(ANS->x0,A->x0,UI);
    mpz_mod(ANS->x0,ANS->x0,prime);
}

void Fp_sub_mpz(Fp *ANS,Fp *A,mpz_t B){
    Fp_cost.Fp_add_mpz++;
    mpz_cost.mpz_add++;
    mpz_sub(ANS->x0,A->x0,B);
    mpz_mod(ANS->x0,ANS->x0,prime);
}

void Fp_inv(Fp *ANS,Fp *A){
    Fp_cost.Fp_inv++;
    mpz_cost.mpz_invert++;
    mpz_invert(ANS->x0,A->x0,prime);
}
*/
//---------------------------------------------------------------//

int  Fp_legendre(Fp *A){
    return mpz_legendre(A->x0,prime);
}

int  Fp_isCNR(Fp *A){
    Fp tmp;
    Fp_init(&tmp);
    mpz_t exp;
    mpz_init(exp);
    mpz_sub_ui(exp,prime,1);
    mpz_tdiv_q_ui(exp,exp,3);
    Fp_pow(&tmp,A,exp);
    if(Fp_cmp_one(&tmp)==0){
        mpz_clear(exp);
        Fp_clear(&tmp);
        return 1;
    }else{
        mpz_clear(exp);
        Fp_clear(&tmp);
        return -1;
    }
}

void Fp_sqrt(Fp *ANS,Fp *A){
    Fp x,y,t,k,n,tmp;
    Fp_init(&x);
    Fp_init(&y);
    Fp_init(&t);
    Fp_init(&k);
    Fp_init(&n);
    Fp_init(&tmp);
    unsigned long int e,m;
    mpz_t exp,q,z,result;
    mpz_init(exp);
    mpz_init(q);
    mpz_init(z);
    mpz_init(result);
    gmp_randstate_t state;
    gmp_randinit_default (state);
    gmp_randseed_ui(state,(unsigned long)time(NULL));
    Fp_set_random(&n,state);
    while(Fp_legendre(&n)!=-1){
        Fp_set_random(&n,state);
    }
    mpz_sub_ui(q,prime,1);
    mpz_mod_ui(result,q,2);
    e=0;
    while(mpz_cmp_ui(result,0)==0){
        mpz_tdiv_q_ui(q,q,2);
        mpz_mod_ui(result,q,2);
        e++;
    }
    Fp_pow(&y,&n,q);
    mpz_set_ui(z,e);
    mpz_sub_ui(exp,q,1);
    mpz_tdiv_q_ui(exp,exp,2);
    Fp_pow(&x,A,exp);
    Fp_mul(&tmp,&x,&x);
    Fp_mul(&k,&tmp,A);
    Fp_mul(&x,&x,A);
    while(mpz_cmp_ui(k.x0,1)!=0){
        m=1;
        mpz_ui_pow_ui(exp,2,m);
        Fp_pow(&tmp,&k,exp);
        while(mpz_cmp_ui(tmp.x0,1)!=0){
            m++;
            mpz_ui_pow_ui(exp,2,m);
            Fp_pow(&tmp,&k,exp);
        }
        mpz_sub_ui(exp,z,m);
        mpz_sub_ui(exp,exp,1);
        mpz_ui_pow_ui(result,2,mpz_get_ui(exp));
        Fp_pow(&t,&y,result);
        Fp_mul(&y,&t,&t);
        mpz_set_ui(z,m);
        Fp_mul(&x,&x,&t);
        Fp_mul(&k,&k,&y);
    }
    Fp_set(ANS,&x);
    mpz_clear(exp);
    mpz_clear(q);
    mpz_clear(z);
    mpz_clear(result);
    Fp_clear(&x);
    Fp_clear(&y);
    Fp_clear(&t);
    Fp_clear(&k);
    Fp_clear(&n);
    Fp_clear(&tmp);
}

void Fp_pow(Fp *ANS,Fp *A,mpz_t scalar){
    int i,length;
    length=(int)mpz_sizeinbase(scalar,2);
    char binary[length];
    mpz_get_str(binary,2,scalar);
    Fp tmp;
    Fp_init(&tmp);
    Fp_set(&tmp,A);
    for(i=1; binary[i]!='\0'; i++){
        Fp_mul(&tmp,&tmp,&tmp);
        if(binary[i]=='1'){
            Fp_mul(&tmp,A,&tmp);
        }
    }
    Fp_set(ANS,&tmp);
    Fp_clear(&tmp);
}


int  Fp_cmp(Fp *A,Fp *B){
    if(mpz_cmp(A->x0,B->x0)==0){
        return 0;   
    }
    return 1;
}

int  Fp_cmp_ui(Fp *A,unsigned long int UI){
    if(mpz_cmp_ui(A->x0,UI)==0){
        return 0;
    }
    return 1;
}

int  Fp_cmp_mpz(Fp *A,mpz_t B){
    if(mpz_cmp(A->x0,B)==0){
        return 0;
    }
    return 1;
}

int  Fp_cmp_zero(Fp *A){
    if(mpz_cmp_ui(A->x0,0)==0){
        return 0;
    }
    return 1;
}

int  Fp_cmp_one(Fp *A){
    if(mpz_cmp_ui(A->x0,1)==0){
        return 0;
    }
    return 1;
}

/*----------------------------------------------------------------------------*/
//Fp2
void Fp2_init(Fp2 *A){
    Fp_init(&A->x0);
    Fp_init(&A->x1);
}

void Fp2_clear(Fp2 *A){
    Fp_clear(&A->x0);
    Fp_clear(&A->x1);
}

void Fp2_printf(Fp2 *A,char *str){
    gmp_printf("%s(",str);
    Fp_printf(&A->x0,"");
    gmp_printf(",");
    Fp_printf(&A->x1,"");
    gmp_printf(")");
}

void Fp2_set(Fp2 *ANS,Fp2 *A){
    Fp_set(&ANS->x0,&A->x0);
    Fp_set(&ANS->x1,&A->x1);
}

void Fp2_set_ui(Fp2 *ANS,unsigned long int UI){
    Fp_set_ui(&ANS->x0,UI);
    Fp_set_ui(&ANS->x1,UI);
}    

void Fp2_set_mpz(Fp2 *ANS,mpz_t A){
    Fp_set_mpz(&ANS->x0,A);
    Fp_set_mpz(&ANS->x1,A);
}

void Fp2_set_neg(Fp2 *ANS,Fp2 *A){
    Fp_set_neg(&ANS->x0,&A->x0);
    Fp_set_neg(&ANS->x1,&A->x1);
}

void Fp2_set_random(Fp2 *ANS,gmp_randstate_t state){
    Fp_set_random(&ANS->x0,state);
    Fp_set_random(&ANS->x1,state);
}

void Fp2_mul(Fp2 *ANS,Fp2 *A,Fp2 *B){
	Fp_mul(&TMP1_FP,&A->x0,&B->x0);
	Fp_mul(&TMP2_FP,&A->x1,&B->x1);
	Fp_add(&TMP3_FP,&A->x0,&A->x1);
	Fp_add(&TMP4_FP,&B->x0,&B->x1);
	Fp_sub(&ANS->x0,&TMP1_FP,&TMP2_FP);
	Fp_mul(&ANS->x1,&TMP3_FP,&TMP4_FP);
	Fp_sub(&ANS->x1,&ANS->x1,&TMP1_FP);
	Fp_sub(&ANS->x1,&ANS->x1,&TMP2_FP);
}

void Fp2_mul_ui(Fp2 *ANS,Fp2 *A,unsigned long int UI){
    Fp_mul_ui(&ANS->x0,&A->x0,UI);
    Fp_mul_ui(&ANS->x1,&A->x1,UI);
}

void Fp2_mul_mpz(Fp2 *ANS,Fp2 *A,mpz_t B){
    Fp_mul_mpz(&ANS->x0,&A->x0,B);
    Fp_mul_mpz(&ANS->x1,&A->x1,B);
}

void Fp2_mul_basis(Fp2 *ANS,Fp2 *A){
    Fp_set(&TMP1_FP,&A->x0);
    Fp_sub(&ANS->x0,&TMP1_FP,&A->x1);
    Fp_add(&ANS->x1,&TMP1_FP,&A->x1);
}

void Fp2_inv_basis(Fp2 *ANS,Fp2 *A){
    Fp_set(&TMP1_FP,&A->x0);
    Fp_set(&TMP2_FP,&A->x1);
    Fp_add(&ANS->x0,&TMP1_FP,&TMP2_FP);
    Fp_mul_mpz(&ANS->x0,&ANS->x0,Alpha_1_inv.x0.x0);
    Fp_sub(&ANS->x1,&TMP2_FP,&TMP1_FP);
    Fp_mul_mpz(&ANS->x1,&ANS->x1,Alpha_1_inv.x0.x0);
}

void Fp2_sqr(Fp2 *ANS,Fp2 *A){
    Fp_add(&TMP1_FP,&A->x0,&A->x1);
	Fp_sub(&TMP2_FP,&A->x0,&A->x1);
	Fp_mul(&ANS->x1,&A->x0,&A->x1);
	Fp_add(&ANS->x1,&ANS->x1,&ANS->x1);
	Fp_mul(&ANS->x0,&TMP1_FP,&TMP2_FP);
}

void Fp2_add(Fp2 *ANS,Fp2 *A,Fp2 *B){
    Fp_add(&ANS->x0,&A->x0,&B->x0);
    Fp_add(&ANS->x1,&A->x1,&B->x1);
}

void Fp2_add_ui(Fp2 *ANS,Fp2 *A,unsigned long int UI){
    Fp_add_ui(&ANS->x0,&A->x0,UI);
    Fp_add_ui(&ANS->x1,&A->x1,UI);
}

void Fp2_add_mpz(Fp2 *ANS,Fp2 *A,mpz_t B){
    Fp_add_mpz(&ANS->x0,&A->x0,B);
    Fp_add_mpz(&ANS->x1,&A->x1,B);
}

void Fp2_sub(Fp2 *ANS,Fp2 *A,Fp2 *B){
    Fp_sub(&ANS->x0,&A->x0,&B->x0);
    Fp_sub(&ANS->x1,&A->x1,&B->x1);
}

void Fp2_sub_ui(Fp2 *ANS,Fp2 *A,unsigned long int UI){
    Fp_sub_ui(&ANS->x0,&A->x0,UI);
    Fp_sub_ui(&ANS->x1,&A->x1,UI);
}

void Fp2_sub_mpz(Fp2 *ANS,Fp2 *A,mpz_t B){
    Fp_sub_mpz(&ANS->x0,&A->x0,B);
    Fp_sub_mpz(&ANS->x1,&A->x1,B);
}

void Fp2_inv(Fp2 *ANS,Fp2 *A){
    Fp_set(&TMP1_FP,&A->x0);
    Fp_set_neg(&TMP2_FP,&A->x1);
    Fp_mul(&TMP3_FP,&TMP1_FP,&A->x0);
    Fp_mul(&TMP4_FP,&TMP2_FP,&A->x1);
    Fp_sub(&TMP3_FP,&TMP3_FP,&TMP4_FP);
    Fp_inv(&TMP3_FP,&TMP3_FP);
    Fp_mul(&ANS->x0,&TMP1_FP,&TMP3_FP);
    Fp_mul(&ANS->x1,&TMP2_FP,&TMP3_FP);
}

int  Fp2_legendre(Fp2 *A){
    Fp2 tmp;
    Fp2_init(&tmp);
    
    mpz_t exp;
    mpz_init(exp);
    mpz_pow_ui(exp,prime,2);
    mpz_sub_ui(exp,exp,1);
    mpz_tdiv_q_ui(exp,exp,2);
    Fp2_pow(&tmp,A,exp);
    
    if(Fp2_cmp_one(&tmp)==0){
        mpz_clear(exp);
        Fp2_clear(&tmp);
        return 1;
    }else{
        mpz_clear(exp);
        Fp2_clear(&tmp);
        return -1;
    }
}

int  Fp2_isCNR(Fp2 *A){
    Fp2 tmp;
    Fp2_init(&tmp);
    mpz_t exp;
    mpz_init(exp);
    
    mpz_pow_ui(exp,prime,2);
    mpz_sub_ui(exp,exp,1);
    mpz_tdiv_q_ui(exp,exp,3);
    Fp2_pow(&tmp,A,exp);
    
    if(Fp2_cmp_one(&tmp)==0){
        mpz_clear(exp);
        Fp2_clear(&tmp);
        return 1;
    }else{
        mpz_clear(exp);
        Fp2_clear(&tmp);
        return -1;
    }

}

void Fp2_sqrt(Fp2 *ANS,Fp2 *A){
    Fp2 x,y,t,k,n,tmp;
    Fp2_init(&x);
    Fp2_init(&y);
    Fp2_init(&t);
    Fp2_init(&k);
    Fp2_init(&n);
    Fp2_init(&tmp);
    unsigned long int e,m;
    mpz_t exp,q,z,result;
    mpz_init(exp);
    mpz_init(q);
    mpz_init(z);
    mpz_init(result);
    gmp_randstate_t state;
	gmp_randinit_default(state);
	gmp_randseed_ui(state,(unsigned long)time(NULL));
    
    Fp2_set_random(&n,state);
    while(Fp2_legendre(&n)!=-1){
        Fp2_set_random(&n,state);
    }
    mpz_pow_ui(q,prime,2);
    mpz_sub_ui(q,q,1);
    mpz_mod_ui(result,q,2);
    e=0;
    while(mpz_cmp_ui(result,0)==0){
        mpz_tdiv_q_ui(q,q,2);
        mpz_mod_ui(result,q,2);
        e++;
    }
    Fp2_pow(&y,&n,q);
    mpz_set_ui(z,e);
    mpz_sub_ui(exp,q,1);
    mpz_tdiv_q_ui(exp,exp,2);
    Fp2_pow(&x,A,exp);
    Fp2_mul(&tmp,&x,&x);
    Fp2_mul(&k,&tmp,A);
    Fp2_mul(&x,&x,A);
    while(Fp2_cmp_one(&k)!=0){
        m=1;
        mpz_ui_pow_ui(exp,2,m);
        Fp2_pow(&tmp,&k,exp);
        while(Fp2_cmp_one(&tmp)!=0){
            m++;
            mpz_ui_pow_ui(exp,2,m);
            Fp2_pow(&tmp,&k,exp);
        }
        mpz_sub_ui(exp,z,m);
        mpz_sub_ui(exp,exp,1);
        mpz_ui_pow_ui(result,2,mpz_get_ui(exp));
        Fp2_pow(&t,&y,result);
        Fp2_mul(&y,&t,&t);
        mpz_set_ui(z,m);
        Fp2_mul(&x,&x,&t);
        Fp2_mul(&k,&k,&y);
    }
    Fp2_set(ANS,&x);
    
    mpz_clear(exp);
    mpz_clear(q);
    mpz_clear(z);
    mpz_clear(result);
    Fp2_clear(&x);
    Fp2_clear(&y);
    Fp2_clear(&t);
    Fp2_clear(&k);
    Fp2_clear(&n);
    Fp2_clear(&tmp);
}

void Fp2_pow(Fp2 *ANS,Fp2 *A,mpz_t scalar){
    int i,length;
    length=(int)mpz_sizeinbase(scalar,2);
    char binary[length];
    mpz_get_str(binary,2,scalar);
    Fp2 tmp;
    Fp2_init(&tmp);
    Fp2_set(&tmp,A);
    
    for(i=1; binary[i]!='\0'; i++){
        Fp2_sqr(&tmp,&tmp);
        if(binary[i]=='1'){
            Fp2_mul(&tmp,A,&tmp);
        }
    }
    
    Fp2_set(ANS,&tmp);
    Fp2_clear(&tmp);
}


int  Fp2_cmp(Fp2 *A,Fp2 *B){
    if(Fp_cmp(&A->x0,&B->x0)==0 && Fp_cmp(&A->x1,&B->x1)==0){
        return 0;   
    }
    return 1;
}

int  Fp2_cmp_ui(Fp2 *A,unsigned long int UI){
    if(Fp_cmp_ui(&A->x0,UI)==0 && Fp_cmp_ui(&A->x1,UI)==0){
        return 0;
    }
    return 1;
}

int  Fp2_cmp_mpz(Fp2 *A,mpz_t B){
    if(Fp_cmp_mpz(&A->x0,B)==0 && Fp_cmp_mpz(&A->x1,B)==0){
        return 0;
    }
    return 1;
}

int  Fp2_cmp_zero(Fp2 *A){
    if(Fp_cmp_zero(&A->x0)==0 && Fp_cmp_zero(&A->x1)==0){
        return 0;
    }
    return 1;
}

int  Fp2_cmp_one(Fp2 *A){
    if(Fp_cmp_one(&A->x0)==0 && Fp_cmp_zero(&A->x1)==0){
        return 0;
    }
    return 1;
}

/*----------------------------------------------------------------------------*/
//Fp6
void Fp6_init(Fp6 *A){
    Fp2_init(&A->x0);
    Fp2_init(&A->x1);
    Fp2_init(&A->x2);
}

void Fp6_clear(Fp6 *A){
    Fp2_clear(&A->x0);
    Fp2_clear(&A->x1);
    Fp2_clear(&A->x2);
}

void Fp6_printf(Fp6 *A,char *str){
    gmp_printf("%s(",str);
    Fp2_printf(&A->x0,"");
    gmp_printf(",");
    Fp2_printf(&A->x1,"");
    gmp_printf(",");
    Fp2_printf(&A->x2,"");
    gmp_printf(")");
}

void Fp6_set(Fp6 *ANS,Fp6 *A){
    Fp2_set(&ANS->x0,&A->x0);
    Fp2_set(&ANS->x1,&A->x1);
    Fp2_set(&ANS->x2,&A->x2);
}

void Fp6_set_ui(Fp6 *ANS,unsigned long int UI){
    Fp2_set_ui(&ANS->x0,UI);
    Fp2_set_ui(&ANS->x1,UI);
    Fp2_set_ui(&ANS->x2,UI);
}

void Fp6_set_mpz(Fp6 *ANS,mpz_t A){
    Fp2_set_mpz(&ANS->x0,A);
    Fp2_set_mpz(&ANS->x1,A);
    Fp2_set_mpz(&ANS->x2,A);
}

void Fp6_set_neg(Fp6 *ANS,Fp6 *A){
    Fp2_set_neg(&ANS->x0,&A->x0);
    Fp2_set_neg(&ANS->x1,&A->x1);
    Fp2_set_neg(&ANS->x2,&A->x2);
}

void Fp6_set_random(Fp6 *ANS,gmp_randstate_t state){
    Fp2_set_random(&ANS->x0,state);
    Fp2_set_random(&ANS->x1,state);
    Fp2_set_random(&ANS->x2,state);
}

void Fp6_mul_ui(Fp6 *ANS,Fp6 *A,unsigned long int UI){
    Fp2_mul_ui(&ANS->x0,&A->x0,UI);
    Fp2_mul_ui(&ANS->x1,&A->x1,UI);
    Fp2_mul_ui(&ANS->x2,&A->x2,UI);
}

void Fp6_mul_mpz(Fp6 *ANS,Fp6 *A,mpz_t B){
    Fp2_mul_mpz(&ANS->x0,&A->x0,B);
    Fp2_mul_mpz(&ANS->x1,&A->x1,B);
    Fp2_mul_mpz(&ANS->x2,&A->x2,B);
}

void Fp6_mul_basis(Fp6 *ANS,Fp6 *A){
    Fp2_mul_basis(&ANS->x0,&A->x0);
    Fp2_mul_basis(&ANS->x1,&A->x1);
    Fp2_mul_basis(&ANS->x2,&A->x2);
}

void Fp6_add(Fp6 *ANS,Fp6 *A,Fp6 *B){
    Fp2_add(&ANS->x0,&A->x0,&B->x0);
    Fp2_add(&ANS->x1,&A->x1,&B->x1);
    Fp2_add(&ANS->x2,&A->x2,&B->x2);
}

void Fp6_add_ui(Fp6 *ANS,Fp6 *A,unsigned long int UI){
    Fp2_add_ui(&ANS->x0,&A->x0,UI);
    Fp2_add_ui(&ANS->x1,&A->x1,UI);
    Fp2_add_ui(&ANS->x2,&A->x2,UI);
}

void Fp6_add_mpz(Fp6 *ANS,Fp6 *A,mpz_t B){
    Fp2_add_mpz(&ANS->x0,&A->x0,B);
    Fp2_add_mpz(&ANS->x1,&A->x1,B);
    Fp2_add_mpz(&ANS->x2,&A->x2,B);
}

void Fp6_sub(Fp6 *ANS,Fp6 *A,Fp6 *B){
    Fp2_sub(&ANS->x0,&A->x0,&B->x0);
    Fp2_sub(&ANS->x1,&A->x1,&B->x1);
    Fp2_sub(&ANS->x2,&A->x2,&B->x2);
}

void Fp6_sub_ui(Fp6 *ANS,Fp6 *A,unsigned long int UI){
    Fp2_sub_ui(&ANS->x0,&A->x0,UI);
    Fp2_sub_ui(&ANS->x1,&A->x1,UI);
    Fp2_sub_ui(&ANS->x2,&A->x2,UI);
}

void Fp6_sub_mpz(Fp6 *ANS,Fp6 *A,mpz_t B){
    Fp2_sub_mpz(&ANS->x0,&A->x0,B);
    Fp2_sub_mpz(&ANS->x1,&A->x1,B);
    Fp2_sub_mpz(&ANS->x2,&A->x2,B);
}

int  Fp6_cmp(Fp6 *A,Fp6 *B){
    if(Fp2_cmp(&A->x0,&B->x0)==0 && Fp2_cmp(&A->x1,&B->x1)==0 && Fp2_cmp(&A->x2,&B->x2)==0){
        return 0;   
    }
    return 1;
}

int  Fp6_cmp_ui(Fp6 *A,unsigned long int UI){
    if(Fp2_cmp_ui(&A->x0,UI)==0 && Fp2_cmp_ui(&A->x1,UI)==0 && Fp2_cmp_ui(&A->x2,UI)==0){
        return 0;
    }
    return 1;
}

int  Fp6_cmp_mpz(Fp6 *A,mpz_t B){
    if(Fp2_cmp_mpz(&A->x0,B)==0 && Fp2_cmp_mpz(&A->x1,B)==0 && Fp2_cmp_mpz(&A->x2,B)==0){
        return 0;
    }
    return 1;
}

int  Fp6_cmp_zero(Fp6 *A){
    if(Fp2_cmp_zero(&A->x0)==0 && Fp2_cmp_zero(&A->x1)==0 && Fp2_cmp_zero(&A->x2)==0){
        return 0;
    }
    return 1;
}

/*----------------------------------------------------------------------------*/
//Fp6 karatsuba
void Fp6_mul_karatsuba(Fp6 *ANS,Fp6 *A,Fp6 *B){
    Fp2_mul(&TMP2_FP2,&A->x0,&B->x0);
    Fp2_mul(&TMP3_FP2,&A->x1,&B->x1);
    Fp2_mul(&TMP4_FP2,&A->x2,&B->x2);
    Fp2_add(&TMP5_FP2,&A->x1,&A->x2);
    Fp2_add(&TMP1_FP2,&B->x1,&B->x2);
    Fp2_mul(&TMP5_FP2,&TMP5_FP2,&TMP1_FP2);
    Fp2_add(&TMP6_FP2,&A->x0,&A->x1);
    Fp2_add(&TMP1_FP2,&B->x0,&B->x1);
    Fp2_mul(&TMP6_FP2,&TMP6_FP2,&TMP1_FP2);
    Fp2_add(&TMP7_FP2,&A->x0,&A->x2);
    Fp2_add(&TMP1_FP2,&B->x0,&B->x2);
    Fp2_mul(&TMP7_FP2,&TMP7_FP2,&TMP1_FP2);
    Fp2_sub(&ANS->x0,&TMP5_FP2,&TMP3_FP2);
    Fp2_sub(&ANS->x0,&ANS->x0,&TMP4_FP2);
    Fp2_add(&ANS->x0,&ANS->x0,&ANS->x0);
    Fp2_add(&ANS->x0,&ANS->x0,&TMP2_FP2);
    Fp2_add(&ANS->x1,&TMP4_FP2,&TMP4_FP2);
    Fp2_add(&ANS->x1,&ANS->x1,&TMP6_FP2);
    Fp2_sub(&ANS->x1,&ANS->x1,&TMP2_FP2);
    Fp2_sub(&ANS->x1,&ANS->x1,&TMP3_FP2);
    Fp2_add(&ANS->x2,&TMP3_FP2,&TMP7_FP2);
    Fp2_sub(&ANS->x2,&ANS->x2,&TMP2_FP2);
    Fp2_sub(&ANS->x2,&ANS->x2,&TMP4_FP2);
}

void Fp6_sqr_karatsuba(Fp6 *ANS,Fp6 *A){
    Fp2_sqr(&TMP2_FP2,&A->x0);
    Fp2_add(&TMP4_FP2,&A->x0,&A->x2);
    Fp2_add(&TMP3_FP2,&TMP4_FP2,&A->x1);
    Fp2_sqr(&TMP3_FP2,&TMP3_FP2);
    Fp2_sub(&TMP4_FP2,&TMP4_FP2,&A->x1);
    Fp2_sqr(&TMP4_FP2,&TMP4_FP2);
    Fp2_mul(&TMP5_FP2,&A->x1,&A->x2);
    Fp2_add(&TMP5_FP2,&TMP5_FP2,&TMP5_FP2);
    Fp2_sqr(&TMP6_FP2,&A->x2);
    Fp2_add(&TMP7_FP2,&TMP3_FP2,&TMP4_FP2);
    Fp2_mul_mpz(&TMP7_FP2,&TMP7_FP2,Two_inv);
    Fp2_add(&ANS->x0,&TMP5_FP2,&TMP5_FP2);
    Fp2_add(&ANS->x0,&ANS->x0,&TMP2_FP2);
    Fp2_add(&ANS->x1,&TMP6_FP2,&TMP6_FP2);
    Fp2_add(&ANS->x1,&ANS->x1,&TMP3_FP2);
    Fp2_sub(&ANS->x1,&ANS->x1,&TMP5_FP2);
    Fp2_sub(&ANS->x1,&ANS->x1,&TMP7_FP2);
    Fp2_sub(&ANS->x2,&TMP7_FP2,&TMP6_FP2);
    Fp2_sub(&ANS->x2,&ANS->x2,&TMP2_FP2);
}

void Fp6_inv_karatsuba(Fp6 *ANS,Fp6 *A){
    Fp2_sqr(&TMP1_FP2,&A->x0);
    Fp2_mul(&TMP4_FP2,&A->x1,&A->x2);
    Fp2_add(&TMP4_FP2,&TMP4_FP2,&TMP4_FP2);
    Fp2_sub(&TMP1_FP2,&TMP1_FP2,&TMP4_FP2);
    Fp2_sqr(&TMP2_FP2,&A->x2);
    Fp2_add(&TMP2_FP2,&TMP2_FP2,&TMP2_FP2);
    Fp2_mul(&TMP4_FP2,&A->x0,&A->x1);
    Fp2_sub(&TMP2_FP2,&TMP2_FP2,&TMP4_FP2);
    Fp2_sqr(&TMP3_FP2,&A->x1);
    Fp2_mul(&TMP4_FP2,&A->x0,&A->x2);
    Fp2_sub(&TMP3_FP2,&TMP3_FP2,&TMP4_FP2);
    Fp2_mul(&TMP5_FP2,&A->x1,&TMP3_FP2);
    Fp2_mul(&TMP4_FP2,&A->x2,&TMP2_FP2);
    Fp2_add(&TMP5_FP2,&TMP5_FP2,&TMP4_FP2);
    Fp2_add(&TMP5_FP2,&TMP5_FP2,&TMP5_FP2);
    Fp2_mul(&TMP4_FP2,&A->x0,&TMP1_FP2);
    Fp2_add(&TMP5_FP2,&TMP5_FP2,&TMP4_FP2);
    Fp2_inv(&TMP5_FP2,&TMP5_FP2);
    Fp2_mul(&ANS->x0,&TMP5_FP2,&TMP1_FP2);
    Fp2_mul(&ANS->x1,&TMP5_FP2,&TMP2_FP2);
    Fp2_mul(&ANS->x2,&TMP5_FP2,&TMP3_FP2);
}

int  Fp6_legendre_karatsuba(Fp6 *A){
    mpz_t exp;
    mpz_init(exp);
    Fp6 tmp;
    Fp6_init(&tmp);
    
    mpz_pow_ui(exp,prime,6);
    mpz_sub_ui(exp,exp,1);
    mpz_tdiv_q_ui(exp,exp,2);
    Fp6_pow_karatsuba(&tmp,A,exp);
    
    if(Fp6_cmp_one_karatsuba(&tmp)==0){
        mpz_clear(exp);
        Fp6_clear(&tmp);
        return 1;
    }else{
        mpz_clear(exp);
        Fp6_clear(&tmp);
        return -1;
    }
}

int  Fp6_isCNR_karatsuba(Fp6 *A){
    Fp6 tmp;
    Fp6_init(&tmp);
    mpz_t exp;
    mpz_init(exp);
    
    mpz_pow_ui(exp,prime,6);
    mpz_sub_ui(exp,exp,1);
    mpz_tdiv_q_ui(exp,exp,3);
    Fp6_pow_karatsuba(&tmp,A,exp);
    
    if(Fp6_cmp_one_karatsuba(&tmp)==0){
        mpz_clear(exp);
        Fp6_clear(&tmp);
        return 1;
    }else{
        mpz_clear(exp);
        Fp6_clear(&tmp);
        return -1;
    }
}

void Fp6_sqrt_karatsuba(Fp6 *ANS,Fp6 *A){
    Fp6 tmp1,tmp2;
    Fp6_init(&tmp1);
    Fp6_init(&tmp2);
    mpz_t exp,buf;
    mpz_init(exp);
    mpz_init(buf);
    
    Fp2_set(&tmp1.x0,&A->x0);
    Fp2_mul_mpz(&tmp1.x1,&A->x1,frobenius_constant_karatsuba[f_p4][1].x0.x0);
    Fp2_mul_mpz(&tmp1.x2,&A->x2,frobenius_constant_karatsuba[f_p4][2].x0.x0);
    Fp2_set(&tmp2.x0,&A->x0);
    Fp2_mul_mpz(&tmp2.x1,&A->x1,frobenius_constant_karatsuba[f_p2][1].x0.x0);
    Fp2_mul_mpz(&tmp2.x2,&A->x2,frobenius_constant_karatsuba[f_p2][2].x0.x0);
    Fp6_mul_karatsuba(&tmp1,&tmp1,&tmp2);
    Fp6_mul_karatsuba(&tmp1,&tmp1,A);
    Fp6_set_ui(&tmp2,0);
    Fp2_sqrt(&tmp2.x0,&tmp1.x0);
    Fp2_inv(&tmp2.x0,&tmp2.x0);
    Fp2_set(&tmp2.x0,&tmp2.x0);
    mpz_pow_ui(exp,prime,8);
    mpz_pow_ui(buf,prime,4);
    mpz_add(exp,exp,buf);
    mpz_add_ui(exp,exp,2);
    mpz_tdiv_q_ui(exp,exp,2);
    Fp6_pow_karatsuba(&tmp1,A,exp);
    Fp6_mul_karatsuba(&tmp1,&tmp1,&tmp2);
    Fp6_set(ANS,&tmp1);
    
    mpz_clear(exp);
    mpz_clear(buf);
    Fp6_clear(&tmp1);
    Fp6_clear(&tmp2);  
}

void Fp6_pow_karatsuba(Fp6 *ANS,Fp6 *A,mpz_t scalar){
    int i,length;
    length=(int)mpz_sizeinbase(scalar,2);
    char binary[length];
    mpz_get_str(binary,2,scalar);
    Fp6 tmp;
    Fp6_init(&tmp);
    Fp6_set(&tmp,A);
    
    for(i=1; binary[i]!='\0'; i++){
        Fp6_sqr_karatsuba(&tmp,&tmp);
        if(binary[i]=='1'){
            Fp6_mul_karatsuba(&tmp,A,&tmp);
        }
    }
    
    Fp6_set(ANS,&tmp);
    Fp6_clear(&tmp);
}

int  Fp6_cmp_one_karatsuba(Fp6 *A){
    if(Fp2_cmp_one(&A->x0)==0 && Fp2_cmp_zero(&A->x1)==0 && Fp2_cmp_zero(&A->x2)==0){
        return 0;
    }
    return 1;
}

/*----------------------------------------------------------------------------*/
//Fp6 CVMA
void Fp6_mul_CVMA(Fp6 *ANS,Fp6 *A,Fp6 *B){
    Fp2_mul(&TMP1_FP2,&A->x0,&B->x0);
    Fp2_mul(&TMP2_FP2,&A->x1,&B->x1);
    Fp2_mul(&TMP3_FP2,&A->x2,&B->x2);
    Fp2_sub(&TMP7_FP2,&A->x0,&A->x1);
    Fp2_sub(&TMP4_FP2,&B->x1,&B->x0);
    Fp2_mul(&TMP4_FP2,&TMP4_FP2,&TMP7_FP2);
    Fp2_sub(&TMP7_FP2,&A->x1,&A->x2);
    Fp2_sub(&TMP5_FP2,&B->x2,&B->x1);
    Fp2_mul(&TMP5_FP2,&TMP5_FP2,&TMP7_FP2);
    Fp2_sub(&TMP7_FP2,&A->x0,&A->x2);
    Fp2_sub(&TMP6_FP2,&B->x2,&B->x0);
    Fp2_mul(&TMP6_FP2,&TMP6_FP2,&TMP7_FP2);
    Fp2_add(&ANS->x0,&TMP4_FP2,&TMP5_FP2);
    Fp2_sub(&ANS->x0,&ANS->x0,&TMP1_FP2);
    Fp2_add(&ANS->x1,&TMP5_FP2,&TMP6_FP2);
    Fp2_sub(&ANS->x1,&ANS->x1,&TMP2_FP2);
    Fp2_add(&ANS->x2,&TMP4_FP2,&TMP6_FP2);
    Fp2_sub(&ANS->x2,&ANS->x2,&TMP3_FP2);
}

void Fp6_sqr_CVMA(Fp6 *ANS,Fp6 *A){
    Fp2_sub(&TMP1_FP2,&A->x1,&A->x0);
    Fp2_add(&TMP2_FP2,&TMP1_FP2,&A->x2);
    Fp2_sqr(&TMP2_FP2,&TMP2_FP2);
    Fp2_sub(&TMP3_FP2,&A->x2,&A->x0);
    Fp2_sub(&TMP5_FP2,&A->x1,&A->x2);
    Fp2_mul(&TMP3_FP2,&TMP3_FP2,&TMP5_FP2);
    Fp2_sqr(&TMP4_FP2,&TMP1_FP2);
    Fp2_sqr(&TMP5_FP2,&A->x1);
    Fp2_add(&TMP4_FP2,&TMP4_FP2,&TMP5_FP2);
    Fp2_add(&TMP5_FP2,&A->x0,&A->x2);
    Fp2_mul(&TMP5_FP2,&TMP5_FP2,&TMP1_FP2);
    Fp2_add(&TMP5_FP2,&TMP5_FP2,&TMP3_FP2);
    Fp2_sub(&ANS->x0,&TMP5_FP2,&TMP4_FP2);
    Fp2_add(&ANS->x1,&TMP3_FP2,&TMP3_FP2);
    Fp2_sub(&ANS->x1,&ANS->x1,&TMP4_FP2);
    Fp2_sub(&ANS->x2,&TMP5_FP2,&TMP2_FP2);
}

void Fp6_inv_CVMA(Fp6 *ANS,Fp6 *A){
    Fp2_sub(&TMP4_FP2,&A->x2,&A->x0);
    Fp2_sub(&TMP5_FP2,&A->x0,&A->x1);
    Fp2_sub(&TMP6_FP2,&A->x2,&A->x1);
    Fp2_sqr(&TMP7_FP2,&TMP4_FP2);
    Fp2_sqr(&TMP8_FP2,&TMP5_FP2);
    Fp2_sqr(&TMP9_FP2,&TMP6_FP2);
    Fp2_mul(&TMP10_FP2,&A->x1,&A->x2);
    Fp2_mul(&TMP11_FP2,&A->x0,&A->x2);
    Fp2_mul(&TMP12_FP2,&A->x0,&A->x1);
    Fp2_sub(&TMP1_FP2,&TMP7_FP2,&TMP10_FP2);
    Fp2_sub(&TMP2_FP2,&TMP8_FP2,&TMP11_FP2);
    Fp2_sub(&TMP3_FP2,&TMP9_FP2,&TMP12_FP2);
    Fp2_add(&TMP14_FP2,&TMP5_FP2,&A->x0);
    Fp2_mul(&TMP14_FP2,&TMP14_FP2,&TMP1_FP2);
    Fp2_add(&TMP13_FP2,&TMP5_FP2,&TMP6_FP2);
    Fp2_mul(&TMP13_FP2,&TMP13_FP2,&TMP2_FP2);
    Fp2_sub(&TMP14_FP2,&TMP13_FP2,&TMP14_FP2);
    Fp2_mul(&TMP13_FP2,&TMP6_FP2,&TMP3_FP2);
    Fp2_sub(&TMP14_FP2,&TMP14_FP2,&TMP13_FP2);
    Fp2_inv(&TMP14_FP2,&TMP14_FP2);
    Fp2_set_neg(&TMP14_FP2,&TMP14_FP2);
    Fp2_mul(&ANS->x0,&TMP1_FP2,&TMP14_FP2);
    Fp2_mul(&ANS->x1,&TMP2_FP2,&TMP14_FP2);
    Fp2_mul(&ANS->x2,&TMP3_FP2,&TMP14_FP2);
}

int  Fp6_legendre_CVMA(Fp6 *A){
    mpz_t exp;
    mpz_init(exp);
    Fp6 tmp;
    Fp6_init(&tmp);
    
    mpz_pow_ui(exp,prime,6);
    mpz_sub_ui(exp,exp,1);
    mpz_tdiv_q_ui(exp,exp,2);
    Fp6_pow_CVMA(&tmp,A,exp);
    
    if(Fp6_cmp_one_CVMA(&tmp)==0){
        mpz_clear(exp);
        Fp6_clear(&tmp);
        return 1;
    }else{
        mpz_clear(exp);
        Fp6_clear(&tmp);
        return -1;
    }
}

int  Fp6_isCNR_CVMA(Fp6 *A){
    Fp6 tmp;
    Fp6_init(&tmp);
    mpz_t exp;
    mpz_init(exp);
    
    mpz_pow_ui(exp,prime,6);
    mpz_sub_ui(exp,exp,1);
    mpz_tdiv_q_ui(exp,exp,3);
    Fp6_pow_CVMA(&tmp,A,exp);
    
    if(Fp6_cmp_one_CVMA(&tmp)==0){
        mpz_clear(exp);
        Fp6_clear(&tmp);
        return 1;
    }else{
        mpz_clear(exp);
        Fp6_clear(&tmp);
        return -1;
    }
}

void Fp6_sqrt_CVMA(Fp6 *ANS,Fp6 *A){
    Fp6 tmp1,tmp2;
    Fp6_init(&tmp1);
    Fp6_init(&tmp2);
    mpz_t exp,buf;
    mpz_init(exp);
    mpz_init(buf);
    
    Fp6_frobenius_map_p4_CVMA(&tmp1,A);
    Fp6_frobenius_map_p2_CVMA(&tmp2,A);
    Fp6_mul_CVMA(&tmp1,&tmp1,&tmp2);
    Fp6_mul_CVMA(&tmp1,&tmp1,A);
    Fp2_sqrt(&tmp2.x0,&tmp1.x0);
    Fp2_inv(&tmp2.x0,&tmp2.x0);
    Fp2_set_neg(&tmp2.x0,&tmp2.x0);
    Fp2_set(&tmp2.x1,&tmp2.x0);
    Fp2_set(&tmp2.x2,&tmp2.x0);
    mpz_pow_ui(exp,prime,4);
    mpz_pow_ui(buf,prime,2);
    mpz_add(exp,exp,buf);
    mpz_add_ui(exp,exp,2);
    mpz_tdiv_q_ui(exp,exp,2);
    Fp6_pow_CVMA(&tmp1,A,exp);
    Fp6_mul_CVMA(&tmp1,&tmp1,&tmp2);
    Fp6_set(ANS,&tmp1);
    
    mpz_clear(exp);
    mpz_clear(buf);
    Fp6_clear(&tmp1);
    Fp6_clear(&tmp2); 
}

void Fp6_pow_CVMA(Fp6 *ANS,Fp6 *A,mpz_t scalar){
    int i,length;
    length=(int)mpz_sizeinbase(scalar,2);
    char binary[length];
    mpz_get_str(binary,2,scalar);
    Fp6 tmp;
    Fp6_init(&tmp);
    Fp6_set(&tmp,A);
    
    for(i=1; binary[i]!='\0'; i++){
        Fp6_sqr_CVMA(&tmp,&tmp);
        if(binary[i]=='1'){
            Fp6_mul_CVMA(&tmp,A,&tmp);
        }
    }
    
    Fp6_set(ANS,&tmp);
    Fp6_clear(&tmp);
}

int  Fp6_cmp_one_CVMA(Fp6 *A){
    Fp6 tmp;
    Fp6_init(&tmp);
    Fp6_set_neg(&tmp,A);
    
    if(Fp2_cmp_one(&tmp.x0)==0 && Fp2_cmp_one(&tmp.x1)==0 && Fp2_cmp_one(&tmp.x2)==0){
        Fp6_clear(&tmp);
        return 0;
    }
    
    Fp6_clear(&tmp);
    return 1;
}

void Fp6_frobenius_map_p1_CVMA(Fp6 *ANS,Fp6 *A){
    Fp6_set(&TMP1_FP6,A);
    Fp_set(&ANS->x0.x0,&TMP1_FP6.x2.x0);
    Fp_set_neg(&ANS->x0.x1,&TMP1_FP6.x2.x1);
    Fp_set(&ANS->x1.x0,&TMP1_FP6.x0.x0);
    Fp_set_neg(&ANS->x1.x1,&TMP1_FP6.x0.x1);
    Fp_set(&ANS->x2.x0,&TMP1_FP6.x1.x0);
    Fp_set_neg(&ANS->x2.x1,&TMP1_FP6.x1.x1);
}

void Fp6_frobenius_map_p2_CVMA(Fp6 *ANS,Fp6 *A){
    Fp6_set(&TMP1_FP6,A);
    Fp2_set(&ANS->x0,&TMP1_FP6.x1);
    Fp2_set(&ANS->x1,&TMP1_FP6.x2);
    Fp2_set(&ANS->x2,&TMP1_FP6.x0);
}

void Fp6_frobenius_map_p3_CVMA(Fp6 *ANS,Fp6 *A){
    Fp_set(&ANS->x0.x0,&A->x0.x0);
    Fp_set_neg(&ANS->x0.x1,&A->x0.x1);
    Fp_set(&ANS->x1.x0,&A->x1.x0);
    Fp_set_neg(&ANS->x1.x1,&A->x1.x1);
    Fp_set(&ANS->x2.x0,&A->x2.x0);
    Fp_set_neg(&ANS->x2.x1,&A->x2.x1);
}

void Fp6_frobenius_map_p4_CVMA(Fp6 *ANS,Fp6 *A){
    Fp6_set(&TMP1_FP6,A);
    Fp2_set(&ANS->x0,&TMP1_FP6.x2);
    Fp2_set(&ANS->x1,&TMP1_FP6.x0);
    Fp2_set(&ANS->x2,&TMP1_FP6.x1);
}

void Fp6_frobenius_map_p5_CVMA(Fp6 *ANS,Fp6 *A){
    Fp6_set(&TMP1_FP6,A);
    Fp_set(&ANS->x0.x0,&TMP1_FP6.x1.x0);
    Fp_set_neg(&ANS->x0.x1,&TMP1_FP6.x1.x1);
    Fp_set(&ANS->x1.x0,&TMP1_FP6.x2.x0);
    Fp_set_neg(&ANS->x1.x1,&TMP1_FP6.x2.x1);
    Fp_set(&ANS->x2.x0,&TMP1_FP6.x0.x0);
    Fp_set_neg(&ANS->x2.x1,&TMP1_FP6.x0.x1);
}

/*----------------------------------------------------------------------------*/
//Fp12
void Fp12_init(Fp12 *A){
    Fp6_init(&A->x0);
    Fp6_init(&A->x1);
}

void Fp12_clear(Fp12 *A){
    Fp6_clear(&A->x0);
    Fp6_clear(&A->x1);
}

void Fp12_printf(Fp12 *A,char *str){
    gmp_printf("%s(",str);
    Fp6_printf(&A->x0,"");
    gmp_printf(",");
    Fp6_printf(&A->x1,"");
    gmp_printf(")");
}

void Fp12_set(Fp12 *ANS,Fp12 *A){
    Fp6_set(&ANS->x0,&A->x0);
    Fp6_set(&ANS->x1,&A->x1);
}

void Fp12_set_ui(Fp12 *ANS,unsigned long int UI){
    Fp6_set_ui(&ANS->x0,UI);
    Fp6_set_ui(&ANS->x1,UI);
}

void Fp12_set_mpz(Fp12 *ANS,mpz_t A){
    Fp6_set_mpz(&ANS->x0,A);
    Fp6_set_mpz(&ANS->x1,A);    
}

void Fp12_set_neg(Fp12 *ANS,Fp12 *A){
    Fp6_set_neg(&ANS->x0,&A->x0);
    Fp6_set_neg(&ANS->x1,&A->x1);
}

void Fp12_set_random(Fp12 *ANS,gmp_randstate_t state){
    Fp6_set_random(&ANS->x0,state);
    Fp6_set_random(&ANS->x1,state);
}

void Fp12_mul_ui(Fp12 *ANS,Fp12 *A,unsigned long int UI){
    Fp6_mul_ui(&ANS->x0,&A->x0,UI);
    Fp6_mul_ui(&ANS->x1,&A->x1,UI);
}

void Fp12_mul_mpz(Fp12 *ANS,Fp12 *A,mpz_t B){
    Fp6_mul_mpz(&ANS->x0,&A->x0,B);
    Fp6_mul_mpz(&ANS->x1,&A->x1,B);
}

void Fp12_add(Fp12 *ANS,Fp12 *A,Fp12 *B){
    Fp6_add(&ANS->x0,&A->x0,&B->x0);
    Fp6_add(&ANS->x1,&A->x1,&B->x1);
}

void Fp12_add_ui(Fp12 *ANS,Fp12 *A,unsigned long int UI){
    Fp6_add_ui(&ANS->x0,&A->x0,UI);
    Fp6_add_ui(&ANS->x1,&A->x1,UI);
}

void Fp12_add_mpz(Fp12 *ANS,Fp12 *A,mpz_t B){
    Fp6_add_mpz(&ANS->x0,&ANS->x0,B);
    Fp6_add_mpz(&ANS->x1,&ANS->x1,B);
}

void Fp12_sub(Fp12 *ANS,Fp12 *A,Fp12 *B){
    Fp6_sub(&ANS->x0,&A->x0,&B->x0);
    Fp6_sub(&ANS->x1,&A->x1,&B->x1);
}

void Fp12_sub_ui(Fp12 *ANS,Fp12 *A,unsigned long int UI){
    Fp6_sub_ui(&ANS->x0,&ANS->x0,UI);
    Fp6_sub_ui(&ANS->x1,&ANS->x1,UI);
}

void Fp12_sub_mpz(Fp12 *ANS,Fp12 *A,mpz_t B){
    Fp6_sub_mpz(&ANS->x0,&ANS->x0,B);
    Fp6_sub_mpz(&ANS->x1,&ANS->x1,B);
}

int  Fp12_cmp(Fp12 *A,Fp12 *B){
    if(Fp6_cmp(&A->x0,&B->x0)==0 && Fp6_cmp(&A->x1,&B->x1)==0){
        return 0;   
    }
    return 1;
}

int  Fp12_cmp_ui(Fp12 *A,unsigned long int UI){
    if(Fp6_cmp_ui(&A->x0,UI)==0 && Fp6_cmp_ui(&A->x1,UI)==0){
        return 0;
    }
    return 1;
}

int  Fp12_cmp_mpz(Fp12 *A,mpz_t B){
    if(Fp6_cmp_mpz(&A->x0,B)==0 && Fp6_cmp_mpz(&A->x1,B)==0){
        return 0;
    }
    return 1;
}

int  Fp12_cmp_zero(Fp12 *A){
    if(Fp6_cmp_zero(&A->x0)==0 && Fp6_cmp_zero(&A->x1)==0){
        return 0;
    }
    return 1;
}

/*----------------------------------------------------------------------------*/
//Fp12 karatsuba
void Fp12_mul_karatsuba(Fp12 *ANS,Fp12 *A,Fp12 *B){
	Fp6_mul_karatsuba(&TMP2_FP6,&A->x1,&B->x1);
	Fp6_add(&TMP1_FP6,&A->x0,&A->x1);
	Fp6_add(&ANS->x1,&B->x0,&B->x1);
	Fp6_mul_karatsuba(&ANS->x1,&TMP1_FP6,&ANS->x1);
	Fp6_mul_karatsuba(&TMP1_FP6,&A->x0,&B->x0);
	Fp6_mul_basis(&ANS->x0,&TMP2_FP6);
	Fp6_add(&ANS->x0,&ANS->x0,&TMP1_FP6);
	Fp6_sub(&ANS->x1,&ANS->x1,&TMP1_FP6);
	Fp6_sub(&ANS->x1,&ANS->x1,&TMP2_FP6);
}

void Fp12_sqr_karatsuba(Fp12 *ANS,Fp12 *A){
	Fp6_add(&TMP1_FP6,&A->x0,&A->x1);
	Fp6_mul_basis(&TMP2_FP6,&A->x1);
	Fp6_add(&TMP2_FP6,&TMP2_FP6,&A->x0);
	Fp6_mul_karatsuba(&TMP3_FP6,&A->x0,&A->x1);
	Fp6_mul_karatsuba(&ANS->x0,&TMP1_FP6,&TMP2_FP6);
	Fp6_sub(&ANS->x0,&ANS->x0,&TMP3_FP6);
	Fp6_mul_basis(&TMP1_FP6,&TMP3_FP6);
	Fp6_sub(&ANS->x0,&ANS->x0,&TMP1_FP6);
	Fp6_add(&ANS->x1,&TMP3_FP6,&TMP3_FP6);
}

void Fp12_inv_karatsuba(Fp12 *ANS,Fp12 *A){
    Fp6_set(&TMP1_FP6,&A->x0);
    Fp6_set_neg(&TMP2_FP6,&A->x1);
    Fp6_mul_karatsuba(&TMP3_FP6,&TMP1_FP6,&A->x0);
    Fp6_mul_karatsuba(&TMP4_FP6,&TMP2_FP6,&A->x1);
    Fp6_mul_basis(&TMP4_FP6,&TMP4_FP6);
    Fp6_add(&TMP3_FP6,&TMP3_FP6,&TMP4_FP6);
    Fp6_inv_karatsuba(&TMP3_FP6,&TMP3_FP6);
    Fp6_mul_karatsuba(&ANS->x0,&TMP1_FP6,&TMP3_FP6);
    Fp6_mul_karatsuba(&ANS->x1,&TMP2_FP6,&TMP3_FP6);
}

int  Fp12_legendre_karatsuba(Fp12 *A){
    mpz_t exp;
    mpz_init(exp);
    Fp12 tmp;
    Fp12_init(&tmp);
    
    mpz_pow_ui(exp,prime,12);
    mpz_sub_ui(exp,exp,1);
    mpz_tdiv_q_ui(exp,exp,2);
    Fp12_pow_karatsuba(&tmp,A,exp);
    
    if(Fp12_cmp_one_karatsuba(&tmp)==0){
        mpz_clear(exp);
        Fp12_clear(&tmp);
        return 1;
    }else{
        mpz_clear(exp);
        Fp12_clear(&tmp);
        return -1;
    }
}

void Fp12_sqrt_karatsuba(Fp12 *ANS,Fp12 *A){
    Fp12 x,y,t,k,n,tmp;
    Fp12_init(&x);
    Fp12_init(&y);
    Fp12_init(&t);
    Fp12_init(&k);
    Fp12_init(&n);
    Fp12_init(&tmp);
    unsigned long int e,m;
    mpz_t exp,q,z,result;
    mpz_init(exp);
    mpz_init(q);
    mpz_init(z);
    mpz_init(result);
    gmp_randstate_t state;
	gmp_randinit_default (state);
	gmp_randseed_ui(state,(unsigned long)time(NULL));
    
    Fp12_set_random(&n,state);
    while(Fp12_legendre_karatsuba(&n)!=-1){
        Fp12_set_random(&n,state);
    }
    mpz_pow_ui(q,prime,12);
    mpz_sub_ui(q,q,1);
    mpz_mod_ui(result,q,2);
    e=0;
    while(mpz_cmp_ui(result,0)==0){
        mpz_tdiv_q_ui(q,q,2);
        mpz_mod_ui(result,q,2);
        e++;
    }
    Fp12_pow_karatsuba(&y,&n,q);
    mpz_set_ui(z,e);    
    mpz_sub_ui(exp,q,1);
    mpz_tdiv_q_ui(exp,exp,2);
    Fp12_pow_karatsuba(&x,A,exp);
    Fp12_mul_karatsuba(&tmp,&x,&x);
    Fp12_mul_karatsuba(&k,&tmp,A);
    Fp12_mul_karatsuba(&x,&x,A);
    while(Fp12_cmp_one_karatsuba(&k)!=0){
        m=1;
        mpz_ui_pow_ui(exp,2,m);
        Fp12_pow_karatsuba(&tmp,&k,exp);
        while(Fp12_cmp_one_karatsuba(&tmp)!=0){
            m++;
            mpz_ui_pow_ui(exp,2,m);
            Fp12_pow_karatsuba(&tmp,&k,exp);
        }
        mpz_sub_ui(exp,z,m);
        mpz_sub_ui(exp,exp,1);
        mpz_ui_pow_ui(result,2,mpz_get_ui(exp));
        Fp12_pow_karatsuba(&t,&y,result);
        Fp12_mul_karatsuba(&y,&t,&t);
        mpz_set_ui(z,m);
        Fp12_mul_karatsuba(&x,&x,&t); 
        Fp12_mul_karatsuba(&k,&k,&y);
    }
    Fp12_set(ANS,&x);
    
    mpz_clear(q);
    mpz_clear(z);
    mpz_clear(exp);
    mpz_clear(result);
    Fp12_clear(&x);
    Fp12_clear(&y);
    Fp12_clear(&t);
    Fp12_clear(&k);
    Fp12_clear(&n);
    Fp12_clear(&tmp);
}

void Fp12_pow_karatsuba(Fp12 *ANS,Fp12 *A,mpz_t scalar){
    int i,length;
    length=(int)mpz_sizeinbase(scalar,2);
    char binary[length];
    mpz_get_str(binary,2,scalar);
    Fp12 tmp;
    Fp12_init(&tmp);
    Fp12_set(&tmp,A);
    
    for(i=1; binary[i]!='\0'; i++){
        Fp12_sqr_karatsuba(&tmp,&tmp);
        if(binary[i]=='1'){
            Fp12_mul_karatsuba(&tmp,A,&tmp);
        }
    }
    
    Fp12_set(ANS,&tmp);
    Fp12_clear(&tmp);
}

int  Fp12_cmp_one_karatsuba(Fp12 *A){
    if(Fp6_cmp_one_karatsuba(&A->x0)==0 && Fp6_cmp_zero(&A->x1)==0){
        return 0;
    }
    return 1;
}

void Fp12_frobenius_map_p1_karatsuba(Fp12 *ANS,Fp12 *A){
	Fp_set(&ANS->x0.x0.x0,&A->x0.x0.x0);
	Fp_set_neg(&ANS->x0.x0.x1,&A->x0.x0.x1);
	Fp_set(&ANS->x0.x1.x0,&A->x0.x1.x0);
	Fp_set_neg(&ANS->x0.x1.x1,&A->x0.x1.x1);
	Fp2_mul_mpz(&ANS->x0.x1,&ANS->x0.x1,frobenius_constant_karatsuba[f_p1][1].x0.x0);
	Fp_set(&ANS->x0.x2.x0,&A->x0.x2.x0);
	Fp_set_neg(&ANS->x0.x2.x1,&A->x0.x2.x1);
	Fp2_mul_mpz(&ANS->x0.x2,&ANS->x0.x2,frobenius_constant_karatsuba[f_p1][2].x0.x0);
	Fp_set(&ANS->x1.x0.x0,&A->x1.x0.x0);
	Fp_set_neg(&ANS->x1.x0.x1,&A->x1.x0.x1);
	Fp2_mul_basis(&ANS->x1.x0,&ANS->x1.x0);
	Fp2_mul_mpz(&ANS->x1.x0,&ANS->x1.x0,frobenius_constant_karatsuba[f_p1][3].x0.x0);
	Fp_set(&ANS->x1.x1.x0,&A->x1.x1.x0);
	Fp_set_neg(&ANS->x1.x1.x1,&A->x1.x1.x1);
	Fp2_mul_basis(&ANS->x1.x1,&ANS->x1.x1);
	Fp2_mul_mpz(&ANS->x1.x1,&ANS->x1.x1,frobenius_constant_karatsuba[f_p1][4].x0.x0);
	Fp_set(&ANS->x1.x2.x0,&A->x1.x2.x0);
	Fp_set_neg(&ANS->x1.x2.x1,&A->x1.x2.x1);
	Fp2_mul_basis(&ANS->x1.x2,&ANS->x1.x2);
	Fp2_mul_mpz(&ANS->x1.x2,&ANS->x1.x2,frobenius_constant_karatsuba[f_p1][5].x0.x0);
}

void Fp12_frobenius_map_p2_karatsuba(Fp12 *ANS,Fp12 *A){
	Fp2_set(&ANS->x0.x0,&A->x0.x0);
	Fp2_mul_mpz(&ANS->x0.x1,&A->x0.x1,frobenius_constant_karatsuba[f_p2][1].x0.x0);
	Fp2_mul_mpz(&ANS->x0.x2,&A->x0.x2,frobenius_constant_karatsuba[f_p2][2].x0.x0);
	Fp2_set_neg(&ANS->x1.x0,&A->x1.x0);
	Fp2_mul_mpz(&ANS->x1.x1,&A->x1.x1,frobenius_constant_karatsuba[f_p2][4].x0.x0);
	Fp2_mul_mpz(&ANS->x1.x2,&A->x1.x2,frobenius_constant_karatsuba[f_p2][5].x0.x0);
}

void Fp12_frobenius_map_p3_karatsuba(Fp12 *ANS,Fp12 *A){
	Fp_set(&ANS->x0.x0.x0,&A->x0.x0.x0);
	Fp_set_neg(&ANS->x0.x0.x1,&A->x0.x0.x1);
	Fp_set(&ANS->x0.x1.x0,&A->x0.x1.x0);
	Fp_set_neg(&ANS->x0.x1.x1,&A->x0.x1.x1);
	Fp_set(&ANS->x0.x2.x0,&A->x0.x2.x0);
	Fp_set_neg(&ANS->x0.x2.x1,&A->x0.x2.x1);
	Fp_set(&ANS->x1.x0.x0,&A->x1.x0.x0);
	Fp_set_neg(&ANS->x1.x0.x1,&A->x1.x0.x1);
	Fp2_mul_basis(&ANS->x1.x0,&ANS->x1.x0);
	Fp2_mul_mpz(&ANS->x1.x0,&ANS->x1.x0,frobenius_constant_karatsuba[f_p3][3].x0.x0);
	Fp_set(&ANS->x1.x1.x0,&A->x1.x1.x0);
	Fp_set_neg(&ANS->x1.x1.x1,&A->x1.x1.x1);
	Fp2_mul_basis(&ANS->x1.x1,&ANS->x1.x1);
	Fp2_mul_mpz(&ANS->x1.x1,&ANS->x1.x1,frobenius_constant_karatsuba[f_p3][4].x0.x0);
	Fp_set(&ANS->x1.x2.x0,&A->x1.x2.x0);
	Fp_set_neg(&ANS->x1.x2.x1,&A->x1.x2.x1);
	Fp2_mul_basis(&ANS->x1.x2,&ANS->x1.x2);
	Fp2_mul_mpz(&ANS->x1.x2,&ANS->x1.x2,frobenius_constant_karatsuba[f_p3][5].x0.x0);
}

void Fp12_frobenius_map_p4_karatsuba(Fp12 *ANS,Fp12 *A){
	Fp2_set(&ANS->x0.x0,&A->x0.x0);
	Fp2_mul_mpz(&ANS->x0.x1,&A->x0.x1,frobenius_constant_karatsuba[f_p4][1].x0.x0);
	Fp2_mul_mpz(&ANS->x0.x2,&A->x0.x2,frobenius_constant_karatsuba[f_p4][2].x0.x0);
	Fp2_set(&ANS->x1.x0,&A->x1.x0);
	Fp2_mul_mpz(&ANS->x1.x1,&A->x1.x1,frobenius_constant_karatsuba[f_p4][4].x0.x0);
	Fp2_mul_mpz(&ANS->x1.x2,&A->x1.x2,frobenius_constant_karatsuba[f_p4][5].x0.x0);
}

void Fp12_frobenius_map_p6_karatsuba(Fp12 *ANS,Fp12 *A){
	Fp6_set(&ANS->x0,&A->x0);
	Fp6_set_neg(&ANS->x1,&A->x1);
}

void Fp12_frobenius_map_p8_karatsuba(Fp12 *ANS,Fp12 *A){
	Fp2_set(&ANS->x0.x0,&A->x0.x0);
	Fp2_mul_mpz(&ANS->x0.x1,&A->x0.x1,frobenius_constant_karatsuba[f_p8][1].x0.x0);
	Fp2_mul_mpz(&ANS->x0.x2,&A->x0.x2,frobenius_constant_karatsuba[f_p8][2].x0.x0);
	Fp2_set(&ANS->x1.x0,&A->x1.x0);
	Fp2_mul_mpz(&ANS->x1.x1,&A->x1.x1,frobenius_constant_karatsuba[f_p8][4].x0.x0);
	Fp2_mul_mpz(&ANS->x1.x2,&A->x1.x2,frobenius_constant_karatsuba[f_p8][5].x0.x0);
}

void Fp12_frobenius_map_p10_karatsuba(Fp12 *ANS,Fp12 *A){
	Fp2_set(&ANS->x0.x0,&A->x0.x0);
	Fp2_mul_mpz(&ANS->x0.x1,&A->x0.x1,frobenius_constant_karatsuba[f_p10][1].x0.x0);
	Fp2_mul_mpz(&ANS->x0.x2,&A->x0.x2,frobenius_constant_karatsuba[f_p10][2].x0.x0);
	Fp2_set_neg(&ANS->x1.x0,&A->x1.x0);
	Fp2_mul_mpz(&ANS->x1.x1,&A->x1.x1,frobenius_constant_karatsuba[f_p10][4].x0.x0);
	Fp2_mul_mpz(&ANS->x1.x2,&A->x1.x2,frobenius_constant_karatsuba[f_p10][5].x0.x0);
}

/*----------------------------------------------------------------------------*/
//Fp12 CVMA
void Fp12_mul_CVMA(Fp12 *ANS,Fp12 *A,Fp12 *B){
    Fp6_mul_CVMA(&TMP2_FP6,&A->x1,&B->x1);
    Fp6_add(&TMP1_FP6,&A->x0,&A->x1);
    Fp6_add(&ANS->x1,&B->x0,&B->x1);
    Fp6_mul_CVMA(&ANS->x1,&TMP1_FP6,&ANS->x1);
    Fp6_mul_CVMA(&TMP1_FP6,&A->x0,&B->x0);
    Fp6_mul_basis(&ANS->x0,&TMP2_FP6);
    Fp6_add(&ANS->x0,&ANS->x0,&TMP1_FP6);
    Fp6_sub(&ANS->x1,&ANS->x1,&TMP1_FP6);
    Fp6_sub(&ANS->x1,&ANS->x1,&TMP2_FP6);
}

void Fp12_sqr_CVMA(Fp12 *ANS,Fp12 *A){
    Fp6_add(&TMP1_FP6,&A->x0,&A->x1);
    Fp6_mul_basis(&TMP2_FP6,&A->x1);
    Fp6_add(&TMP2_FP6,&TMP2_FP6,&A->x0);
    Fp6_mul_CVMA(&TMP3_FP6,&A->x0,&A->x1);
    Fp6_mul_CVMA(&ANS->x0,&TMP1_FP6,&TMP2_FP6);
    Fp6_sub(&ANS->x0,&ANS->x0,&TMP3_FP6);
    Fp6_mul_basis(&TMP1_FP6,&TMP3_FP6);
    Fp6_sub(&ANS->x0,&ANS->x0,&TMP1_FP6);
    Fp6_add(&ANS->x1,&TMP3_FP6,&TMP3_FP6);
}

void Fp12_sqr_cyclotomic_CVMA(Fp12 *ANS,Fp12 *A){
    Fp6_add(&TMP1_FP6,&A->x0,&A->x1);
    Fp6_sqr_CVMA(&TMP1_FP6,&TMP1_FP6);
    Fp6_sqr_CVMA(&TMP2_FP6,&A->x1);
    Fp6_mul_basis(&ANS->x1,&TMP2_FP6);
    Fp6_add(&ANS->x0,&ANS->x1,&ANS->x1);
    Fp_sub_ui(&ANS->x0.x0.x0,&ANS->x0.x0.x0,1);
    Fp_sub_ui(&ANS->x0.x1.x0,&ANS->x0.x1.x0,1);
    Fp_sub_ui(&ANS->x0.x2.x0,&ANS->x0.x2.x0,1);
    Fp6_sub(&ANS->x1,&TMP1_FP6,&ANS->x1);
    Fp6_sub(&ANS->x1,&ANS->x1,&TMP2_FP6);
    Fp_add_ui(&ANS->x1.x0.x0,&ANS->x1.x0.x0,1);
    Fp_add_ui(&ANS->x1.x1.x0,&ANS->x1.x1.x0,1);
    Fp_add_ui(&ANS->x1.x2.x0,&ANS->x1.x2.x0,1);
}

void Fp12_inv_CVMA(Fp12 *ANS,Fp12 *A){
    Fp6_set(&TMP1_FP6,&A->x0);
    Fp6_set_neg(&TMP2_FP6,&A->x1);
    Fp6_mul_CVMA(&TMP3_FP6,&TMP1_FP6,&A->x0);
    Fp6_mul_CVMA(&TMP4_FP6,&TMP2_FP6,&A->x1);
    Fp6_mul_basis(&TMP4_FP6,&TMP4_FP6);
    Fp6_add(&TMP3_FP6,&TMP3_FP6,&TMP4_FP6);
    Fp6_inv_CVMA(&TMP3_FP6,&TMP3_FP6);
    Fp6_mul_CVMA(&ANS->x0,&TMP1_FP6,&TMP3_FP6);
    Fp6_mul_CVMA(&ANS->x1,&TMP2_FP6,&TMP3_FP6);
}

int  Fp12_legendre_CVMA(Fp12 *A){
    mpz_t exp;
    mpz_init(exp);
    Fp12 tmp;
    Fp12_init(&tmp);
    
    mpz_pow_ui(exp,prime,12);
    mpz_sub_ui(exp,exp,1);
    mpz_tdiv_q_ui(exp,exp,2);
    Fp12_pow_CVMA(&tmp,A,exp);
    
    if(Fp12_cmp_one_CVMA(&tmp)==0){
        mpz_clear(exp);
        Fp12_clear(&tmp);
        return 1;
    }else{
        mpz_clear(exp);
        Fp12_clear(&tmp);
        return -1;
    }
}

void Fp12_sqrt_CVMA(Fp12 *ANS,Fp12 *A){
    Fp12 x,y,t,k,n,tmp;
    Fp12_init(&x);
    Fp12_init(&y);
    Fp12_init(&t);
    Fp12_init(&k);
    Fp12_init(&n);
    Fp12_init(&tmp);
    unsigned long int e,m;
    mpz_t exp,q,z,result;
    mpz_init(exp);
    mpz_init(q);
    mpz_init(z);
    mpz_init(result);
    gmp_randstate_t state;
    gmp_randinit_default (state);
    gmp_randseed_ui(state,(unsigned long)time(NULL));
    
    Fp12_set_random(&n,state);
    while(Fp12_legendre_CVMA(&n)!=-1){
        Fp12_set_random(&n,state);
    }
    mpz_pow_ui(q,prime,12);
    mpz_sub_ui(q,q,1);
    mpz_mod_ui(result,q,2);
    e=0;
    while(mpz_cmp_ui(result,0)==0){
        mpz_tdiv_q_ui(q,q,2);
        mpz_mod_ui(result,q,2);
        e++;
    }
    Fp12_pow_CVMA(&y,&n,q);
    mpz_set_ui(z,e);    
    mpz_sub_ui(exp,q,1);
    mpz_tdiv_q_ui(exp,exp,2);
    Fp12_pow_CVMA(&x,A,exp);
    Fp12_mul_CVMA(&tmp,&x,&x);
    Fp12_mul_CVMA(&k,&tmp,A);
    Fp12_mul_CVMA(&x,&x,A);
    while(Fp12_cmp_one_CVMA(&k)!=0){
        m=1;
        mpz_ui_pow_ui(exp,2,m);
        Fp12_pow_CVMA(&tmp,&k,exp);
        while(Fp12_cmp_one_CVMA(&tmp)!=0){
            m++;
            mpz_ui_pow_ui(exp,2,m);
            Fp12_pow_CVMA(&tmp,&k,exp);
        }
        mpz_sub_ui(exp,z,m);
        mpz_sub_ui(exp,exp,1);
        mpz_ui_pow_ui(result,2,mpz_get_ui(exp));
        Fp12_pow_CVMA(&t,&y,result);
        Fp12_mul_CVMA(&y,&t,&t);
        mpz_set_ui(z,m);
        Fp12_mul_CVMA(&x,&x,&t); 
        Fp12_mul_CVMA(&k,&k,&y);
    }
    Fp12_set(ANS,&x);
    
    mpz_clear(q);
    mpz_clear(z);
    mpz_clear(exp);
    mpz_clear(result);
    Fp12_clear(&x);
    Fp12_clear(&y);
    Fp12_clear(&t);
    Fp12_clear(&k);
    Fp12_clear(&n);
    Fp12_clear(&tmp);
}

void Fp12_pow_CVMA(Fp12 *ANS,Fp12 *A,mpz_t scalar){
    int i,length;
    length=(int)mpz_sizeinbase(scalar,2);
    char binary[length];
    mpz_get_str(binary,2,scalar);
    Fp12 tmp;
    Fp12_init(&tmp);
    Fp12_set(&tmp,A);
    
    for(i=1; binary[i]!='\0'; i++){
        Fp12_sqr_CVMA(&tmp,&tmp);
        if(binary[i]=='1'){
            Fp12_mul_CVMA(&tmp,A,&tmp);
        }
    }
    
    Fp12_set(ANS,&tmp);
    Fp12_clear(&tmp);
}

int  Fp12_cmp_one_CVMA(Fp12 *A){
    if(Fp6_cmp_one_CVMA(&A->x0)==0 && Fp6_cmp_zero(&A->x1)==0){
        return 0;
    }
    return 1;
}

void Fp12_frobenius_map_p1_CVMA(Fp12 *ANS,Fp12 *A){
    Fp6_frobenius_map_p1_CVMA(&ANS->x0,&A->x0);
    Fp6_frobenius_map_p1_CVMA(&ANS->x1,&A->x1);
    Fp2_mul_basis(&ANS->x1.x0,&ANS->x1.x0);
    Fp2_mul_mpz(&ANS->x1.x0,&ANS->x1.x0,frobenius_constant_CVMA[f_p1].x0.x0);
    Fp2_mul_basis(&ANS->x1.x1,&ANS->x1.x1);
    Fp2_mul_mpz(&ANS->x1.x1,&ANS->x1.x1,frobenius_constant_CVMA[f_p1].x0.x0);
    Fp2_mul_basis(&ANS->x1.x2,&ANS->x1.x2);
    Fp2_mul_mpz(&ANS->x1.x2,&ANS->x1.x2,frobenius_constant_CVMA[f_p1].x0.x0);
}

void Fp12_frobenius_map_p2_CVMA(Fp12 *ANS,Fp12 *A){
    Fp6_frobenius_map_p2_CVMA(&ANS->x0,&A->x0);
    Fp6_frobenius_map_p2_CVMA(&ANS->x1,&A->x1);
    Fp2_mul_mpz(&ANS->x1.x0,&ANS->x1.x0,frobenius_constant_CVMA[f_p2].x0.x0);
    Fp2_mul_mpz(&ANS->x1.x1,&ANS->x1.x1,frobenius_constant_CVMA[f_p2].x0.x0);
    Fp2_mul_mpz(&ANS->x1.x2,&ANS->x1.x2,frobenius_constant_CVMA[f_p2].x0.x0);
}

void Fp12_frobenius_map_p3_CVMA(Fp12 *ANS,Fp12 *A){
    Fp6_frobenius_map_p3_CVMA(&ANS->x0,&A->x0);
    Fp6_frobenius_map_p3_CVMA(&ANS->x1,&A->x1);
    Fp2_mul_basis(&ANS->x1.x0,&ANS->x1.x0);
    Fp2_mul_mpz(&ANS->x1.x0,&ANS->x1.x0,frobenius_constant_CVMA[f_p3].x0.x0);
    Fp2_mul_basis(&ANS->x1.x1,&ANS->x1.x1);
    Fp2_mul_mpz(&ANS->x1.x1,&ANS->x1.x1,frobenius_constant_CVMA[f_p3].x0.x0);
    Fp2_mul_basis(&ANS->x1.x2,&ANS->x1.x2);
    Fp2_mul_mpz(&ANS->x1.x2,&ANS->x1.x2,frobenius_constant_CVMA[f_p3].x0.x0);
}

void Fp12_frobenius_map_p4_CVMA(Fp12 *ANS,Fp12 *A){
    Fp6_frobenius_map_p4_CVMA(&ANS->x0,&A->x0);
    Fp6_frobenius_map_p4_CVMA(&ANS->x1,&A->x1);
}

void Fp12_frobenius_map_p6_CVMA(Fp12 *ANS,Fp12 *A){
    Fp6_set(&ANS->x0,&A->x0);
    Fp6_set_neg(&ANS->x1,&A->x1);
}

void Fp12_frobenius_map_p8_CVMA(Fp12 *ANS,Fp12 *A){
    Fp6_frobenius_map_p2_CVMA(&ANS->x0,&A->x0);
    Fp6_frobenius_map_p2_CVMA(&ANS->x1,&A->x1);
}

void Fp12_frobenius_map_p10_CVMA(Fp12 *ANS,Fp12 *A){
    Fp6_frobenius_map_p4_CVMA(&ANS->x0,&A->x0);
    Fp6_frobenius_map_p4_CVMA(&ANS->x1,&A->x1);
    Fp6_set_neg(&ANS->x1,&ANS->x1);
}

/*============================================================================*/
/* Elliptic curve                                                             */
/*============================================================================*/
//EFp
void EFp_init(EFp *P){
    Fp_init(&P->x);
    Fp_init(&P->y);
    P->infinity=0;
}

void EFp_set(EFp *ANS,EFp *A){
    Fp_set(&ANS->x,&A->x);
    Fp_set(&ANS->y,&A->y);
    ANS->infinity=A->infinity;
}

void EFp_set_ui(EFp *ANS,unsigned long int UI){
    Fp_set_ui(&ANS->x,UI);
    Fp_set_ui(&ANS->y,UI);
    ANS->infinity=0;
}

void EFp_set_mpz(EFp *ANS,mpz_t A){
    Fp_set_mpz(&ANS->x,A);
    Fp_set_mpz(&ANS->y,A);
    ANS->infinity=0;
}

void EFp_set_neg(EFp *ANS,EFp *A){
    Fp_set(&ANS->x,&A->x);
    Fp_set_neg(&ANS->y,&A->y);
    ANS->infinity=A->infinity;
}

void EFp_clear(EFp *P){
    Fp_clear(&P->x);
    Fp_clear(&P->y);
}

void EFp_printf(EFp *P,char *str){
    printf("%s",str);
    if(P->infinity==0){
        printf("(");
        Fp_printf(&P->x,"");
        printf(",");
        Fp_printf(&P->y,"");
        printf(")");
    }else{
        printf("0");
    }
}

void EFp_rational_point(EFp *P){
    Fp tmp1,tmp2,tmp_x;
    Fp_init(&tmp1);
    Fp_init(&tmp2);
    Fp_init(&tmp_x);
    gmp_randstate_t state;
	gmp_randinit_default(state);
	gmp_randseed_ui(state,(unsigned long)time(NULL));
	
    while(1){
        Fp_set_random(&P->x,state);
        Fp_mul(&tmp1,&P->x,&P->x);
        Fp_mul(&tmp2,&tmp1,&P->x);
        Fp_add_ui(&tmp_x,&tmp2,2);
        if(Fp_legendre(&tmp_x)==1){
            Fp_sqrt(&P->y,&tmp_x);
            break;
        }
    }
    
    Fp_clear(&tmp1);
    Fp_clear(&tmp2);
    Fp_clear(&tmp_x);
}

void EFp_ECD(EFp *ANS,EFp *P){
    if(Fp_cmp_zero(&P->y)==0){
        ANS->infinity=1;
        return;
    }
    
    EFp_set(&TMP1_EFP,P);
    Fp_add(&TMP1_FP,&TMP1_EFP.y,&TMP1_EFP.y);
    Fp_inv(&TMP1_FP,&TMP1_FP);
    Fp_mul(&TMP2_FP,&TMP1_EFP.x,&TMP1_EFP.x);
    Fp_add(&TMP3_FP,&TMP2_FP,&TMP2_FP);
    Fp_add(&TMP2_FP,&TMP2_FP,&TMP3_FP);
    Fp_mul(&TMP3_FP,&TMP1_FP,&TMP2_FP);
    Fp_mul(&TMP1_FP,&TMP3_FP,&TMP3_FP);
    Fp_add(&TMP2_FP,&TMP1_EFP.x,&TMP1_EFP.x);
    Fp_sub(&ANS->x,&TMP1_FP,&TMP2_FP);
    Fp_sub(&TMP1_FP,&TMP1_EFP.x,&ANS->x);
    Fp_mul(&TMP2_FP,&TMP3_FP,&TMP1_FP);
    Fp_sub(&ANS->y,&TMP2_FP,&TMP1_EFP.y);
}

void EFp_ECA(EFp *ANS,EFp *P1,EFp *P2){
    if(P1->infinity==1){
        EFp_set(ANS,P2);
        return;
    }else if(P2->infinity==1){
        EFp_set(ANS,P1);
        return;
    }else if(Fp_cmp(&P1->x,&P2->x)==0){
        if(Fp_cmp(&P1->y,&P2->y)!=0){
            ANS->infinity=1;
            return;
        }else{
            EFp_ECD(ANS,P1);
            return;
        }
    }
    
    EFp_set(&TMP1_EFP,P1);
    EFp_set(&TMP2_EFP,P2);
    Fp_sub(&TMP1_FP,&TMP2_EFP.x,&TMP1_EFP.x);
    Fp_inv(&TMP1_FP,&TMP1_FP);
    Fp_sub(&TMP2_FP,&TMP2_EFP.y,&TMP1_EFP.y);
    Fp_mul(&TMP3_FP,&TMP1_FP,&TMP2_FP);
    Fp_mul(&TMP1_FP,&TMP3_FP,&TMP3_FP);
    Fp_sub(&TMP2_FP,&TMP1_FP,&TMP1_EFP.x);
    Fp_sub(&ANS->x,&TMP2_FP,&TMP2_EFP.x);
    Fp_sub(&TMP1_FP,&TMP1_EFP.x,&ANS->x);
    Fp_mul(&TMP2_FP,&TMP3_FP,&TMP1_FP);
    Fp_sub(&ANS->y,&TMP2_FP,&TMP1_EFP.y);
}

void EFp_SCM(EFp *ANS,EFp *P,mpz_t scalar){
    if(mpz_cmp_ui(scalar,0)==0){
        ANS->infinity=1;
        return;
    }else if(mpz_cmp_ui(scalar,1)==0){
        EFp_set(ANS,P);
        return;
    }
    
    EFp Tmp_P,Next_P;
    EFp_init(&Tmp_P);
    EFp_set(&Tmp_P,P);
    EFp_init(&Next_P);
    int i,length;
    length=(int)mpz_sizeinbase(scalar,2);
    char binary[length];
    mpz_get_str(binary,2,scalar);
    
    EFp_set(&Next_P,&Tmp_P);
    for(i=1; binary[i]!='\0'; i++){
        EFp_ECD(&Next_P,&Next_P);
        if(binary[i]=='1'){
            EFp_ECA(&Next_P,&Next_P,&Tmp_P);
        }
    }
    
    EFp_set(ANS,&Next_P);
    
    EFp_clear(&Next_P);
    EFp_clear(&Tmp_P);
}

void EFp_skew_frobenius_map_p2(EFp *ANS,EFp *A){
    Fp_mul_mpz(&ANS->x,&A->x,epsilon1);
    Fp_set_neg(&ANS->y,&A->y);
}

/*----------------------------------------------------------------------------*/
//EFp2
void EFp2_init(EFp2 *P){
    Fp2_init(&P->x);
    Fp2_init(&P->y);
    P->infinity=0;
}

void EFp2_set(EFp2 *ANS,EFp2 *A){
    Fp2_set(&ANS->x,&A->x);
    Fp2_set(&ANS->y,&A->y);
    ANS->infinity=A->infinity;
}

void EFp2_set_ui(EFp2 *ANS,unsigned long int UI){
    Fp2_set_ui(&ANS->x,UI);
    Fp2_set_ui(&ANS->y,UI);
    ANS->infinity=0;
}

void EFp2_set_mpz(EFp2 *ANS,mpz_t A){
    Fp2_set_mpz(&ANS->x,A);
    Fp2_set_mpz(&ANS->y,A);
    ANS->infinity=0;
}

void EFp2_set_neg(EFp2 *ANS,EFp2 *A){
    Fp2_set(&ANS->x,&A->x);
    Fp2_set_neg(&ANS->y,&A->y);
    ANS->infinity=A->infinity;
}

void EFp2_clear(EFp2 *P){
    Fp2_clear(&P->x);
    Fp2_clear(&P->y);
}

void EFp2_printf(EFp2 *P,char *str){
    printf("%s",str);
    if(P->infinity==0){
        printf("(");
        Fp2_printf(&P->x,"");
        printf(",");
        Fp2_printf(&P->y,"");
        printf(")");
    }else{
        printf("0");
    }
}

void EFp2_rational_point(EFp2 *P){
    Fp2 tmp1,tmp2;
    Fp2_init(&tmp1);
    Fp2_init(&tmp2);
    gmp_randstate_t state;
	gmp_randinit_default(state);
	gmp_randseed_ui(state,(unsigned long)time(NULL));
    
    while(1){
        Fp2_set_random(&P->x,state);
        Fp2_sqr(&tmp1,&P->x);
        Fp2_mul(&tmp2,&tmp1,&P->x);
        mpz_add_ui(tmp2.x0.x0,tmp2.x0.x0,2);
        if(Fp2_legendre(&tmp2)==1){
            Fp2_sqrt(&P->y,&tmp2);
            break;
        }
    }
    
    Fp2_clear(&tmp1);
    Fp2_clear(&tmp2);
}

void EFp2_ECD(EFp2 *ANS,EFp2 *P){
    if(Fp2_cmp_zero(&P->y)==0){
        ANS->infinity=1;
        return;
    }
    
    EFp2_set(&TMP1_EFP2,P);
    Fp2_add(&TMP1_FP2,&TMP1_EFP2.y,&TMP1_EFP2.y);
    Fp2_inv(&TMP1_FP2,&TMP1_FP2);
    Fp2_sqr(&TMP2_FP2,&TMP1_EFP2.x);
    Fp2_add(&TMP3_FP2,&TMP2_FP2,&TMP2_FP2);
    Fp2_add(&TMP2_FP2,&TMP2_FP2,&TMP3_FP2);
    Fp2_mul(&TMP3_FP2,&TMP1_FP2,&TMP2_FP2);
    Fp2_sqr(&TMP1_FP2,&TMP3_FP2);
    Fp2_add(&TMP2_FP2,&TMP1_EFP2.x,&TMP1_EFP2.x);
    Fp2_sub(&ANS->x,&TMP1_FP2,&TMP2_FP2);
    Fp2_sub(&TMP1_FP2,&TMP1_EFP2.x,&ANS->x);
    Fp2_mul(&TMP2_FP2,&TMP3_FP2,&TMP1_FP2);
    Fp2_sub(&ANS->y,&TMP2_FP2,&TMP1_EFP2.y);
}

void EFp2_ECA(EFp2 *ANS,EFp2 *P1,EFp2 *P2){
    if(P1->infinity==1){
        EFp2_set(ANS,P2);
        return;
    }else if(P2->infinity==1){
        EFp2_set(ANS,P1);
        return;
    }else if(Fp2_cmp(&P1->x,&P2->x)==0){
        if(Fp2_cmp(&P1->y,&P2->y)!=0){
            ANS->infinity=1;
            return;
        }else{
            EFp2_ECD(ANS,P1);
            return;
        }
    }
    
    EFp2_set(&TMP1_EFP2,P1);
    EFp2_set(&TMP2_EFP2,P2);
    Fp2_sub(&TMP1_FP2,&TMP2_EFP2.x,&TMP1_EFP2.x);
    Fp2_inv(&TMP1_FP2,&TMP1_FP2);
    Fp2_sub(&TMP2_FP2,&TMP2_EFP2.y,&TMP1_EFP2.y);
    Fp2_mul(&TMP3_FP2,&TMP1_FP2,&TMP2_FP2);
    Fp2_sqr(&TMP1_FP2,&TMP3_FP2);
    Fp2_sub(&TMP2_FP2,&TMP1_FP2,&TMP1_EFP2.x);
    Fp2_sub(&ANS->x,&TMP2_FP2,&TMP2_EFP2.x);
    Fp2_sub(&TMP1_FP2,&TMP1_EFP2.x,&ANS->x);
    Fp2_mul(&TMP2_FP2,&TMP3_FP2,&TMP1_FP2);
    Fp2_sub(&ANS->y,&TMP2_FP2,&TMP1_EFP2.y);
}

void EFp2_SCM(EFp2 *ANS,EFp2 *P,mpz_t scalar){
    if(mpz_cmp_ui(scalar,0)==0){
        ANS->infinity=1;
        return;
    }else if(mpz_cmp_ui(scalar,1)==0){
        EFp2_set(ANS,P);
        return;
    }
    
    EFp2 Tmp_P,Next_P;
    EFp2_init(&Tmp_P);
    EFp2_set(&Tmp_P,P);
    EFp2_init(&Next_P);
    int i,length;
    length=(int)mpz_sizeinbase(scalar,2);
    char binary[length];
    mpz_get_str(binary,2,scalar);
    
    EFp2_set(&Next_P,&Tmp_P);
    for(i=1; binary[i]!='\0'; i++){
        EFp2_ECD(&Next_P,&Next_P);
        if(binary[i]=='1'){
            EFp2_ECA(&Next_P,&Next_P,&Tmp_P);
        }
    }
    EFp2_set(ANS,&Next_P);
    
    EFp2_clear(&Next_P);
    EFp2_clear(&Tmp_P);
}

void EFp2_skew_frobenius_map_p1_karatsuba(EFp2 *ANS,EFp2 *A){
	Fp_set(&TMP1_FP,&A->x.x1);
	Fp_set(&ANS->x.x1,&A->x.x0);
	Fp_set(&ANS->x.x0,&TMP1_FP);
	Fp2_mul_mpz(&ANS->x,&ANS->x,skew_frobenius_constant_karatsuba[f_p1][0].x1.x0);
	Fp_set(&ANS->y.x0,&A->y.x0);
	Fp_set_neg(&ANS->y.x1,&A->y.x1);
	Fp2_mul(&ANS->y,&ANS->y,&skew_frobenius_constant_karatsuba[f_p1][1]);
}

void EFp2_skew_frobenius_map_p2_karatsuba(EFp2 *ANS,EFp2 *A){
	Fp2_mul_mpz(&ANS->x,&A->x,skew_frobenius_constant_karatsuba[f_p2][0].x0.x0);
	Fp2_mul_mpz(&ANS->y,&A->y,skew_frobenius_constant_karatsuba[f_p2][1].x0.x0);
}

void EFp2_skew_frobenius_map_p3_karatsuba(EFp2 *ANS,EFp2 *A){
    Fp_set(&TMP1_FP,&A->x.x1);
	Fp_set(&ANS->x.x1,&A->x.x0);
	Fp_set(&ANS->x.x0,&TMP1_FP);
	Fp2_mul_mpz(&ANS->x,&ANS->x,skew_frobenius_constant_karatsuba[f_p3][0].x1.x0);
	Fp_set(&ANS->y.x0,&A->y.x0);
	Fp_set_neg(&ANS->y.x1,&A->y.x1);
	Fp2_mul(&ANS->y,&ANS->y,&skew_frobenius_constant_karatsuba[f_p3][1]);
}

void EFp2_skew_frobenius_map_p10_karatsuba(EFp2 *ANS,EFp2 *A){
	Fp2_mul_mpz(&ANS->x,&A->x,skew_frobenius_constant_karatsuba[f_p10][0].x0.x0);
	Fp2_mul_mpz(&ANS->y,&A->y,skew_frobenius_constant_karatsuba[f_p10][1].x0.x0);
}

/*----------------------------------------------------------------------------*/
//EFp6
void EFp6_init(EFp6 *P){
    Fp6_init(&P->x);
    Fp6_init(&P->y);
    P->infinity=0;
}

void EFp6_set(EFp6 *ANS,EFp6 *A){
    Fp6_set(&ANS->x,&A->x);
    Fp6_set(&ANS->y,&A->y);
    ANS->infinity=A->infinity;
}

void EFp6_set_ui(EFp6 *ANS,unsigned long int UI){
    Fp6_set_ui(&ANS->x,UI);
    Fp6_set_ui(&ANS->y,UI);
    ANS->infinity=0;
}

void EFp6_set_mpz(EFp6 *ANS,mpz_t A){
    Fp6_set_mpz(&ANS->x,A);
    Fp6_set_mpz(&ANS->y,A);
    ANS->infinity=0;
}

void EFp6_set_neg(EFp6 *ANS,EFp6 *A){
    Fp6_set(&ANS->x,&A->x);
    Fp6_set_neg(&ANS->y,&A->y);
    ANS->infinity=A->infinity;
}

void EFp6_clear(EFp6 *P){
    Fp6_clear(&P->x);
    Fp6_clear(&P->y);
}

void EFp6_printf(EFp6 *P,char *str){
    printf("%s",str);
    if(P->infinity==0){
        printf("(");
        Fp6_printf(&P->x,"");
        printf(",");
        Fp6_printf(&P->y,"");
        printf(")");
    }else{
        printf("0");
    }
}

/*----------------------------------------------------------------------------*/
//EFp6 karatsuba
void EFp6_rational_point_karatsuba(EFp6 *P){
    Fp6 tmp1,tmp2;
    Fp6_init(&tmp1);
    Fp6_init(&tmp2);
    gmp_randstate_t state;
	gmp_randinit_default(state);
	gmp_randseed_ui(state,(unsigned long)time(NULL));
    
    while(1){
        Fp6_set_random(&P->x,state);
        Fp6_sqr_karatsuba(&tmp1,&P->x);
        Fp6_mul_karatsuba(&tmp2,&tmp1,&P->x);
        mpz_add_ui(tmp2.x0.x0.x0,tmp2.x0.x0.x0,2);
        if(Fp6_legendre_karatsuba(&tmp2)==1){
            Fp6_sqrt_karatsuba(&P->y,&tmp2);
            break;
        }
    }
    
    Fp6_clear(&tmp1);
    Fp6_clear(&tmp2);
}

void EFp6_ECD_karatsuba(EFp6 *ANS,EFp6 *P){
    if(Fp6_cmp_zero(&P->y)==0){
        ANS->infinity=1;
        return;
    }
    
    EFp6_set(&TMP1_EFP6,P);
    Fp6_add(&TMP1_FP6,&TMP1_EFP6.y,&TMP1_EFP6.y);
    Fp6_inv_karatsuba(&TMP1_FP6,&TMP1_FP6);
    Fp6_sqr_karatsuba(&TMP2_FP6,&TMP1_EFP6.x);
    Fp6_add(&TMP3_FP6,&TMP2_FP6,&TMP2_FP6);
    Fp6_add(&TMP2_FP6,&TMP2_FP6,&TMP3_FP6);
    Fp6_mul_karatsuba(&TMP3_FP6,&TMP1_FP6,&TMP2_FP6);
    Fp6_sqr_karatsuba(&TMP1_FP6,&TMP3_FP6);
    Fp6_add(&TMP2_FP6,&TMP1_EFP6.x,&TMP1_EFP6.x);
    Fp6_sub(&ANS->x,&TMP1_FP6,&TMP2_FP6);
    Fp6_sub(&TMP1_FP6,&TMP1_EFP6.x,&ANS->x);
    Fp6_mul_karatsuba(&TMP2_FP6,&TMP3_FP6,&TMP1_FP6);
    Fp6_sub(&ANS->y,&TMP2_FP6,&TMP1_EFP6.y);
}

void EFp6_ECA_karatsuba(EFp6 *ANS,EFp6 *P1,EFp6 *P2){
    if(P1->infinity==1){
        EFp6_set(ANS,P2);
        return;
    }else if(P2->infinity==1){
        EFp6_set(ANS,P1);
        return;
    }else if(Fp6_cmp(&P1->x,&P2->x)==0){
        if(Fp6_cmp(&P1->y,&P2->y)!=0){
            ANS->infinity=1;
            return;
        }else{
            EFp6_ECD_karatsuba(ANS,P1);
            return;
        }
    }
    
    EFp6_set(&TMP1_EFP6,P1);
    EFp6_set(&TMP2_EFP6,P2);
    Fp6_sub(&TMP1_FP6,&TMP2_EFP6.x,&TMP1_EFP6.x);
    Fp6_inv_karatsuba(&TMP1_FP6,&TMP1_FP6);
    Fp6_sub(&TMP2_FP6,&TMP2_EFP6.y,&TMP1_EFP6.y);
    Fp6_mul_karatsuba(&TMP3_FP6,&TMP1_FP6,&TMP2_FP6);
    Fp6_sqr_karatsuba(&TMP1_FP6,&TMP3_FP6);
    Fp6_sub(&TMP2_FP6,&TMP1_FP6,&TMP1_EFP6.x);
    Fp6_sub(&ANS->x,&TMP2_FP6,&TMP2_EFP6.x);
    Fp6_sub(&TMP1_FP6,&TMP1_EFP6.x,&ANS->x);
    Fp6_mul_karatsuba(&TMP2_FP6,&TMP3_FP6,&TMP1_FP6);
    Fp6_sub(&ANS->y,&TMP2_FP6,&TMP1_EFP6.y);
}

void EFp6_SCM_karatsuba(EFp6 *ANS,EFp6 *P,mpz_t scalar){
    if(mpz_cmp_ui(scalar,0)==0){
        ANS->infinity=1;
        return;
    }else if(mpz_cmp_ui(scalar,1)==0){
        EFp6_set(ANS,P);
        return;
    }
    
    EFp6 Tmp_P,Next_P;
    EFp6_init(&Tmp_P);
    EFp6_set(&Tmp_P,P);
    EFp6_init(&Next_P);
    int i,length;
    length=(int)mpz_sizeinbase(scalar,2);
    char binary[length];
    mpz_get_str(binary,2,scalar);
    
    EFp6_set(&Next_P,&Tmp_P);
    for(i=1; binary[i]!='\0'; i++){
        EFp6_ECD_karatsuba(&Next_P,&Next_P);
        if(binary[i]=='1'){
            EFp6_ECA_karatsuba(&Next_P,&Next_P,&Tmp_P);
        }
    }
    EFp6_set(ANS,&Next_P);
    
    EFp6_clear(&Next_P);
    EFp6_clear(&Tmp_P);
}

/*----------------------------------------------------------------------------*/
//EFp6 CVMA
void EFp6_rational_point_CVMA(EFp6 *P){
    Fp6 tmp1,tmp2;
    Fp6_init(&tmp1);
    Fp6_init(&tmp2);
    gmp_randstate_t state;
    gmp_randinit_default(state);
    gmp_randseed_ui(state,(unsigned long)time(NULL));
    
    while(1){
        Fp6_set_random(&P->x,state);
        Fp6_sqr_CVMA(&tmp1,&P->x);
        Fp6_mul_CVMA(&tmp2,&tmp1,&P->x);
        Fp_sub_ui(&tmp2.x0.x0,&tmp2.x0.x0,2);
        Fp_sub_ui(&tmp2.x1.x0,&tmp2.x1.x0,2);
        Fp_sub_ui(&tmp2.x2.x0,&tmp2.x2.x0,2);
        if(Fp6_legendre_CVMA(&tmp2)==1){
            Fp6_sqrt_CVMA(&P->y,&tmp2);
            break;
        }
    }
    
    Fp6_clear(&tmp1);
    Fp6_clear(&tmp2);
}

void EFp6_ECD_CVMA(EFp6 *ANS,EFp6 *P){
    if(Fp6_cmp_zero(&P->y)==0){
        ANS->infinity=1;
        return;
    }
    
    EFp6_set(&TMP1_EFP6,P);
    Fp6_add(&TMP1_FP6,&TMP1_EFP6.y,&TMP1_EFP6.y);
    Fp6_inv_CVMA(&TMP1_FP6,&TMP1_FP6);
    Fp6_sqr_CVMA(&TMP2_FP6,&TMP1_EFP6.x);
    Fp6_add(&TMP3_FP6,&TMP2_FP6,&TMP2_FP6);
    Fp6_add(&TMP2_FP6,&TMP2_FP6,&TMP3_FP6);
    Fp6_mul_CVMA(&TMP3_FP6,&TMP1_FP6,&TMP2_FP6);
    Fp6_sqr_CVMA(&TMP1_FP6,&TMP3_FP6);
    Fp6_add(&TMP2_FP6,&TMP1_EFP6.x,&TMP1_EFP6.x);
    Fp6_sub(&ANS->x,&TMP1_FP6,&TMP2_FP6);
    Fp6_sub(&TMP1_FP6,&TMP1_EFP6.x,&ANS->x);
    Fp6_mul_CVMA(&TMP2_FP6,&TMP3_FP6,&TMP1_FP6);
    Fp6_sub(&ANS->y,&TMP2_FP6,&TMP1_EFP6.y);
}

void EFp6_ECA_CVMA(EFp6 *ANS,EFp6 *P1,EFp6 *P2){
    if(P1->infinity==1){
        EFp6_set(ANS,P2);
        return;
    }else if(P2->infinity==1){
        EFp6_set(ANS,P1);
        return;
    }else if(Fp6_cmp(&P1->x,&P2->x)==0){
        if(Fp6_cmp(&P1->y,&P2->y)!=0){
            ANS->infinity=1;
            return;
        }else{
            EFp6_ECD_CVMA(ANS,P1);
            return;
        }
    }
    
    EFp6_set(&TMP1_EFP6,P1);
    EFp6_set(&TMP2_EFP6,P2);
    Fp6_sub(&TMP1_FP6,&TMP2_EFP6.x,&TMP1_EFP6.x);
    Fp6_inv_CVMA(&TMP1_FP6,&TMP1_FP6);
    Fp6_sub(&TMP2_FP6,&TMP2_EFP6.y,&TMP1_EFP6.y);
    Fp6_mul_CVMA(&TMP3_FP6,&TMP1_FP6,&TMP2_FP6);
    Fp6_sqr_CVMA(&TMP1_FP6,&TMP3_FP6);
    Fp6_sub(&TMP2_FP6,&TMP1_FP6,&TMP1_EFP6.x);
    Fp6_sub(&ANS->x,&TMP2_FP6,&TMP2_EFP6.x);
    Fp6_sub(&TMP1_FP6,&TMP1_EFP6.x,&ANS->x);
    Fp6_mul_CVMA(&TMP2_FP6,&TMP3_FP6,&TMP1_FP6);
    Fp6_sub(&ANS->y,&TMP2_FP6,&TMP1_EFP6.y);
}

void EFp6_SCM_CVMA(EFp6 *ANS,EFp6 *P,mpz_t scalar){
    if(mpz_cmp_ui(scalar,0)==0){
        ANS->infinity=1;
        return;
    }else if(mpz_cmp_ui(scalar,1)==0){
        EFp6_set(ANS,P);
        return;
    }
    
    EFp6 Tmp_P,Next_P;
    EFp6_init(&Tmp_P);
    EFp6_set(&Tmp_P,P);
    EFp6_init(&Next_P);
    int i,length;
    length=(int)mpz_sizeinbase(scalar,2);
    char binary[length];
    mpz_get_str(binary,2,scalar);
    
    EFp6_set(&Next_P,&Tmp_P);
    for(i=1; binary[i]!='\0'; i++){
        EFp6_ECD_CVMA(&Next_P,&Next_P);
        if(binary[i]=='1'){
            EFp6_ECA_CVMA(&Next_P,&Next_P,&Tmp_P);
        }
    }
    EFp6_set(ANS,&Next_P);
    
    EFp6_clear(&Next_P);
    EFp6_clear(&Tmp_P);
}

/*----------------------------------------------------------------------------*/
//EFp12
void EFp12_init(EFp12 *P){
    Fp12_init(&P->x);
    Fp12_init(&P->y);
    P->infinity=0;
}

void EFp12_set(EFp12 *ANS,EFp12 *A){
    Fp12_set(&ANS->x,&A->x);
    Fp12_set(&ANS->y,&A->y);
    ANS->infinity=A->infinity;
}

void EFp12_set_ui(EFp12 *ANS,unsigned long int UI){
    Fp12_set_ui(&ANS->x,UI);
    Fp12_set_ui(&ANS->y,UI);
    ANS->infinity=0;
}

void EFp12_set_mpz(EFp12 *ANS,mpz_t A){
    Fp12_set_mpz(&ANS->x,A);
    Fp12_set_mpz(&ANS->y,A);
    ANS->infinity=0;
}

void EFp12_set_neg(EFp12 *ANS,EFp12 *A){
    Fp12_set(&ANS->x,&A->x);
    Fp12_set_neg(&ANS->y,&A->y);
    ANS->infinity=A->infinity;
}

void EFp12_clear(EFp12 *P){
    Fp12_clear(&P->x);
    Fp12_clear(&P->y);
}

void EFp12_printf(EFp12 *P,char *str){
    printf("%s",str);
    if(P->infinity==0){
        printf("(");
        Fp12_printf(&P->x,"");
        printf(",");
        Fp12_printf(&P->y,"");
        printf(")");
    }else{
        printf("0");
    }
}

/*----------------------------------------------------------------------------*/
//EFp12 karatsuba
void EFp12_rational_point_karatsuba(EFp12 *P){
    Fp12 tmp1,tmp2;
    Fp12_init(&tmp1);
    Fp12_init(&tmp2);
    gmp_randstate_t state;
	gmp_randinit_default(state);
	gmp_randseed_ui(state,(unsigned long)time(NULL));
    
    while(1){
        Fp12_set_random(&P->x,state);
        Fp12_sqr_karatsuba(&tmp1,&P->x);
        Fp12_mul_karatsuba(&tmp2,&tmp1,&P->x);
        mpz_add_ui(tmp2.x0.x0.x0.x0,tmp2.x0.x0.x0.x0,2);
        if(Fp12_legendre_karatsuba(&tmp2)==1){
            Fp12_sqrt_karatsuba(&P->y,&tmp2);
            break;
        }
    }
    
    Fp12_clear(&tmp1);
    Fp12_clear(&tmp2);
}

void EFp12_generate_G1_karatsuba(EFp12 *P){
    EFp tmp_P;
    EFp_init(&tmp_P);
    
    EFp_rational_point(&tmp_P);
    EFp12_set_ui(P,0);
    Fp_set(&P->x.x0.x0.x0,&tmp_P.x);
    Fp_set(&P->y.x0.x0.x0,&tmp_P.y);
    P->infinity=tmp_P.infinity;
    
    EFp_clear(&tmp_P);
}

void EFp12_generate_G2_karatsuba(EFp12 *Q){
    EFp12 random_P,P,frobenius_P;
    EFp12_init(&random_P);
    EFp12_init(&P);
    EFp12_init(&frobenius_P);
    mpz_t exp;
    mpz_init(exp);
    
    EFp12_rational_point_karatsuba(&random_P);
    mpz_pow_ui(exp,order,2);
    mpz_tdiv_q(exp,EFp12_total,exp);
    EFp12_SCM_karatsuba(&P,&random_P,exp);
    Fp12_frobenius_map_p1_karatsuba(&frobenius_P.x,&P.x);
    Fp12_frobenius_map_p1_karatsuba(&frobenius_P.y,&P.y);
    EFp12_set_neg(&P,&P);
    EFp12_ECA_karatsuba(Q,&P,&frobenius_P);
    
    mpz_clear(exp);
    EFp12_clear(&random_P);
    EFp12_clear(&P);
    EFp12_clear(&frobenius_P);
}

void EFp12_ECD_karatsuba(EFp12 *ANS,EFp12 *P){
    if(Fp12_cmp_zero(&P->y)==0){
        ANS->infinity=1;
        return;
    }
    
    EFp12_set(&TMP1_EFP12,P);
    Fp12_add(&TMP1_FP12,&TMP1_EFP12.y,&TMP1_EFP12.y);
    Fp12_inv_karatsuba(&TMP1_FP12,&TMP1_FP12);
    Fp12_sqr_karatsuba(&TMP2_FP12,&TMP1_EFP12.x);
    Fp12_add(&TMP3_FP12,&TMP2_FP12,&TMP2_FP12);
    Fp12_add(&TMP2_FP12,&TMP2_FP12,&TMP3_FP12);
    Fp12_mul_karatsuba(&TMP3_FP12,&TMP1_FP12,&TMP2_FP12);
    Fp12_sqr_karatsuba(&TMP1_FP12,&TMP3_FP12);
    Fp12_add(&TMP2_FP12,&TMP1_EFP12.x,&TMP1_EFP12.x);
    Fp12_sub(&ANS->x,&TMP1_FP12,&TMP2_FP12);
    Fp12_sub(&TMP1_FP12,&TMP1_EFP12.x,&ANS->x);
    Fp12_mul_karatsuba(&TMP2_FP12,&TMP3_FP12,&TMP1_FP12);
    Fp12_sub(&ANS->y,&TMP2_FP12,&TMP1_EFP12.y);
}

void EFp12_ECA_karatsuba(EFp12 *ANS,EFp12 *P1,EFp12 *P2){
    if(P1->infinity==1){
        EFp12_set(ANS,P2);
        return;
    }else if(P2->infinity==1){
        EFp12_set(ANS,P1);
        return;
    }else if(Fp12_cmp(&P1->x,&P2->x)==0){
        if(Fp12_cmp(&P1->y,&P2->y)!=0){
            ANS->infinity=1;
            return;
        }else{
            EFp12_ECD_karatsuba(ANS,P1);
            return;
        }
    }
    
    EFp12_set(&TMP1_EFP12,P1);
    EFp12_set(&TMP2_EFP12,P2);
    Fp12_sub(&TMP1_FP12,&TMP2_EFP12.x,&TMP1_EFP12.x);
    Fp12_inv_karatsuba(&TMP1_FP12,&TMP1_FP12);
    Fp12_sub(&TMP2_FP12,&TMP2_EFP12.y,&TMP1_EFP12.y);
    Fp12_mul_karatsuba(&TMP3_FP12,&TMP1_FP12,&TMP2_FP12);
    Fp12_sqr_karatsuba(&TMP1_FP12,&TMP3_FP12);
    Fp12_sub(&TMP2_FP12,&TMP1_FP12,&TMP1_EFP12.x);
    Fp12_sub(&ANS->x,&TMP2_FP12,&TMP2_EFP12.x);
    Fp12_sub(&TMP1_FP12,&TMP1_EFP12.x,&ANS->x);
    Fp12_mul_karatsuba(&TMP2_FP12,&TMP3_FP12,&TMP1_FP12);
    Fp12_sub(&ANS->y,&TMP2_FP12,&TMP1_EFP12.y);
}

void EFp12_SCM_karatsuba(EFp12 *ANS,EFp12 *P,mpz_t scalar){
    if(mpz_cmp_ui(scalar,0)==0){
        ANS->infinity=1;
        return;
    }else if(mpz_cmp_ui(scalar,1)==0){
        EFp12_set(ANS,P);
        return;
    }
    
    EFp12 Tmp_P,Next_P;
    EFp12_init(&Tmp_P);
    EFp12_set(&Tmp_P,P);
    EFp12_init(&Next_P);
    int i,length;
    length=(int)mpz_sizeinbase(scalar,2);
    char binary[length];
    mpz_get_str(binary,2,scalar);
    
    EFp12_set(&Next_P,&Tmp_P);
    for(i=1; binary[i]!='\0'; i++){
        EFp12_ECD_karatsuba(&Next_P,&Next_P);
        if(binary[i]=='1'){
            EFp12_ECA_karatsuba(&Next_P,&Next_P,&Tmp_P);
        }
    }
    EFp12_set(ANS,&Next_P);
    
    EFp12_clear(&Next_P);
    EFp12_clear(&Tmp_P);
}

/*----------------------------------------------------------------------------*/
//EFp12 CVMA
void EFp12_rational_point_CVMA(EFp12 *P){
    Fp12 tmp1,tmp2;
    Fp12_init(&tmp1);
    Fp12_init(&tmp2);
    gmp_randstate_t state;
    gmp_randinit_default(state);
    gmp_randseed_ui(state,(unsigned long)time(NULL));
    
    while(1){
        Fp12_set_random(&P->x,state);
        Fp12_sqr_CVMA(&tmp1,&P->x);
        Fp12_mul_CVMA(&tmp2,&tmp1,&P->x);
        Fp_sub_ui(&tmp2.x0.x0.x0,&tmp2.x0.x0.x0,2);
        Fp_sub_ui(&tmp2.x0.x1.x0,&tmp2.x0.x1.x0,2);
        Fp_sub_ui(&tmp2.x0.x2.x0,&tmp2.x0.x2.x0,2);
        if(Fp12_legendre_CVMA(&tmp2)==1){
            Fp12_sqrt_CVMA(&P->y,&tmp2);
            break;
        }
    }
    
    Fp12_clear(&tmp1);
    Fp12_clear(&tmp2);
}

void EFp12_generate_G1_CVMA(EFp12 *P){
    EFp tmp_P;
    EFp_init(&tmp_P);
    
    EFp_rational_point(&tmp_P);
    EFp12_set_ui(P,0);
    Fp_set_neg(&P->x.x0.x0.x0,&tmp_P.x);
    Fp_set_neg(&P->x.x0.x1.x0,&tmp_P.x);
    Fp_set_neg(&P->x.x0.x2.x0,&tmp_P.x);
    Fp_set_neg(&P->y.x0.x0.x0,&tmp_P.y);
    Fp_set_neg(&P->y.x0.x1.x0,&tmp_P.y);
    Fp_set_neg(&P->y.x0.x2.x0,&tmp_P.y);
    P->infinity=tmp_P.infinity;
    
    EFp_clear(&tmp_P);
}

void EFp12_generate_G2_CVMA(EFp12 *Q){
    EFp12 random_P,P,frobenius_P;
    EFp12_init(&random_P);
    EFp12_init(&P);
    EFp12_init(&frobenius_P);
    mpz_t exp;
    mpz_init(exp);

    EFp12_rational_point_CVMA(&random_P);
    mpz_pow_ui(exp,order,2);
    mpz_tdiv_q(exp,EFp12_total,exp);
    EFp12_SCM_CVMA(&P,&random_P,exp);
    Fp12_frobenius_map_p1_CVMA(&frobenius_P.x,&P.x);
    Fp12_frobenius_map_p1_CVMA(&frobenius_P.y,&P.y);
    EFp12_set_neg(&P,&P);
    EFp12_ECA_CVMA(Q,&P,&frobenius_P);
    
    mpz_clear(exp);
    EFp12_clear(&random_P);
    EFp12_clear(&P);
    EFp12_clear(&frobenius_P);
}

void EFp12_ECD_CVMA(EFp12 *ANS,EFp12 *P){
    if(Fp12_cmp_zero(&P->y)==0){
        ANS->infinity=1;
        return;
    }
    
    EFp12_set(&TMP1_EFP12,P);
    Fp12_add(&TMP1_FP12,&TMP1_EFP12.y,&TMP1_EFP12.y);
    Fp12_inv_CVMA(&TMP1_FP12,&TMP1_FP12);
    Fp12_sqr_CVMA(&TMP2_FP12,&TMP1_EFP12.x);
    Fp12_add(&TMP3_FP12,&TMP2_FP12,&TMP2_FP12);
    Fp12_add(&TMP2_FP12,&TMP2_FP12,&TMP3_FP12);
    Fp12_mul_CVMA(&TMP3_FP12,&TMP1_FP12,&TMP2_FP12);
    Fp12_sqr_CVMA(&TMP1_FP12,&TMP3_FP12);
    Fp12_add(&TMP2_FP12,&TMP1_EFP12.x,&TMP1_EFP12.x);
    Fp12_sub(&ANS->x,&TMP1_FP12,&TMP2_FP12);
    Fp12_sub(&TMP1_FP12,&TMP1_EFP12.x,&ANS->x);
    Fp12_mul_CVMA(&TMP2_FP12,&TMP3_FP12,&TMP1_FP12);
    Fp12_sub(&ANS->y,&TMP2_FP12,&TMP1_EFP12.y);
}

void EFp12_ECA_CVMA(EFp12 *ANS,EFp12 *P1,EFp12 *P2){
    if(P1->infinity==1){
        EFp12_set(ANS,P2);
        return;
    }else if(P2->infinity==1){
        EFp12_set(ANS,P1);
        return;
    }else if(Fp12_cmp(&P1->x,&P2->x)==0){
        if(Fp12_cmp(&P1->y,&P2->y)!=0){
            ANS->infinity=1;
            return;
        }else{
            EFp12_ECD_CVMA(ANS,P1);
            return;
        }
    }
    
    EFp12_set(&TMP1_EFP12,P1);
    EFp12_set(&TMP2_EFP12,P2);
    Fp12_sub(&TMP1_FP12,&TMP2_EFP12.x,&TMP1_EFP12.x);
    Fp12_inv_CVMA(&TMP1_FP12,&TMP1_FP12);
    Fp12_sub(&TMP2_FP12,&TMP2_EFP12.y,&TMP1_EFP12.y);
    Fp12_mul_CVMA(&TMP3_FP12,&TMP1_FP12,&TMP2_FP12);
    Fp12_sqr_CVMA(&TMP1_FP12,&TMP3_FP12);
    Fp12_sub(&TMP2_FP12,&TMP1_FP12,&TMP1_EFP12.x);
    Fp12_sub(&ANS->x,&TMP2_FP12,&TMP2_EFP12.x);
    Fp12_sub(&TMP1_FP12,&TMP1_EFP12.x,&ANS->x);
    Fp12_mul_CVMA(&TMP2_FP12,&TMP3_FP12,&TMP1_FP12);
    Fp12_sub(&ANS->y,&TMP2_FP12,&TMP1_EFP12.y);
}

void EFp12_SCM_CVMA(EFp12 *ANS,EFp12 *P,mpz_t scalar){
    if(mpz_cmp_ui(scalar,0)==0){
        ANS->infinity=1;
        return;
    }else if(mpz_cmp_ui(scalar,1)==0){
        EFp12_set(ANS,P);
        return;
    }
    
    EFp12 Tmp_P,Next_P;
    EFp12_init(&Tmp_P);
    EFp12_set(&Tmp_P,P);
    EFp12_init(&Next_P);
    int i,length;
    length=(int)mpz_sizeinbase(scalar,2);
    char binary[length];
    mpz_get_str(binary,2,scalar);
    
    EFp12_set(&Next_P,&Tmp_P);
    for(i=1; binary[i]!='\0'; i++){
        EFp12_ECD_CVMA(&Next_P,&Next_P);
        if(binary[i]=='1'){
            EFp12_ECA_CVMA(&Next_P,&Next_P,&Tmp_P);
        }
    }
    EFp12_set(ANS,&Next_P);
    
    EFp12_clear(&Next_P);
    EFp12_clear(&Tmp_P);
}

/*============================================================================*/
/* Pairing                                                                    */
/*============================================================================*/
//Basis conversion
void karatsuba_to_CVMA(Fp12 *ANS,Fp12 *A){
    gettimeofday(&tv_start,NULL);
    
    Fp2 tmp0,tmp1;
    Fp2_init(&tmp0);
    Fp2_init(&tmp1);
    Fp12 ans;
    Fp12_init(&ans);
    
    Fp2_mul_mpz(&tmp0,&A->x0.x0,Matrix_karatsuba_to_CVMA.data[0][0].x0);
    Fp2_mul_mpz(&tmp1,&A->x0.x1,Matrix_karatsuba_to_CVMA.data[1][0].x0);
    Fp2_mul_mpz(&ans.x0.x0,&A->x0.x2,Matrix_karatsuba_to_CVMA.data[2][0].x0);
    Fp2_add(&ans.x0.x0,&ans.x0.x0,&tmp0);
    Fp2_add(&ans.x0.x0,&ans.x0.x0,&tmp1);
    Fp2_mul_mpz(&tmp0,&A->x0.x0,Matrix_karatsuba_to_CVMA.data[0][1].x0);
    Fp2_mul_mpz(&tmp1,&A->x0.x1,Matrix_karatsuba_to_CVMA.data[1][1].x0);
    Fp2_mul_mpz(&ans.x0.x1,&A->x0.x2,Matrix_karatsuba_to_CVMA.data[2][1].x0);
    Fp2_add(&ans.x0.x1,&ans.x0.x1,&tmp0);
    Fp2_add(&ans.x0.x1,&ans.x0.x1,&tmp1);
    Fp2_mul_mpz(&tmp0,&A->x0.x0,Matrix_karatsuba_to_CVMA.data[0][2].x0);
    Fp2_mul_mpz(&tmp1,&A->x0.x1,Matrix_karatsuba_to_CVMA.data[1][2].x0);
    Fp2_mul_mpz(&ans.x0.x2,&A->x0.x2,Matrix_karatsuba_to_CVMA.data[2][2].x0);
    Fp2_add(&ans.x0.x2,&ans.x0.x2,&tmp0);
    Fp2_add(&ans.x0.x2,&ans.x0.x2,&tmp1);
    Fp2_mul_mpz(&tmp0,&A->x1.x0,Matrix_karatsuba_to_CVMA.data[0][0].x0);
    Fp2_mul_mpz(&tmp1,&A->x1.x1,Matrix_karatsuba_to_CVMA.data[1][0].x0);
    Fp2_mul_mpz(&ans.x1.x0,&A->x1.x2,Matrix_karatsuba_to_CVMA.data[2][0].x0);
    Fp2_add(&ans.x1.x0,&ans.x1.x0,&tmp0);
    Fp2_add(&ans.x1.x0,&ans.x1.x0,&tmp1);
    Fp2_mul_mpz(&tmp0,&A->x1.x0,Matrix_karatsuba_to_CVMA.data[0][1].x0);
    Fp2_mul_mpz(&tmp1,&A->x1.x1,Matrix_karatsuba_to_CVMA.data[1][1].x0);
    Fp2_mul_mpz(&ans.x1.x1,&A->x1.x2,Matrix_karatsuba_to_CVMA.data[2][1].x0);
    Fp2_add(&ans.x1.x1,&ans.x1.x1,&tmp0);
    Fp2_add(&ans.x1.x1,&ans.x1.x1,&tmp1);
    Fp2_mul_mpz(&tmp0,&A->x1.x0,Matrix_karatsuba_to_CVMA.data[0][2].x0);
    Fp2_mul_mpz(&tmp1,&A->x1.x1,Matrix_karatsuba_to_CVMA.data[1][2].x0);
    Fp2_mul_mpz(&ans.x1.x2,&A->x1.x2,Matrix_karatsuba_to_CVMA.data[2][2].x0);
    Fp2_add(&ans.x1.x2,&ans.x1.x2,&tmp0);
    Fp2_add(&ans.x1.x2,&ans.x1.x2,&tmp1);
    Fp12_set(ANS,&ans);
    
    Fp2_clear(&tmp0);
    Fp2_clear(&tmp1);
    Fp12_clear(&ans);
    
    gettimeofday(&tv_end,NULL);
    BASIS_CONVERSION_KTOC=timedifference_msec(tv_start,tv_end);
}

void CVMA_to_karatsuba(Fp12 *ANS,Fp12 *A){
    gettimeofday(&tv_start,NULL);
    
    Fp2 tmp0,tmp1;
    Fp2_init(&tmp0);
    Fp2_init(&tmp1);
    Fp12 ans;
    Fp12_init(&ans);
    
    Fp2_mul_mpz(&tmp0,&A->x0.x0,Matrix_CVMA_to_karatsuba.data[0][0].x0);
    Fp2_mul_mpz(&tmp1,&A->x0.x1,Matrix_CVMA_to_karatsuba.data[1][0].x0);
    Fp2_mul_mpz(&ans.x0.x0,&A->x0.x2,Matrix_CVMA_to_karatsuba.data[2][0].x0);
    Fp2_add(&ans.x0.x0,&ans.x0.x0,&tmp0);
    Fp2_add(&ans.x0.x0,&ans.x0.x0,&tmp1);
    Fp2_mul_mpz(&tmp0,&A->x0.x0,Matrix_CVMA_to_karatsuba.data[0][1].x0);
    Fp2_mul_mpz(&tmp1,&A->x0.x1,Matrix_CVMA_to_karatsuba.data[1][1].x0);
    Fp2_mul_mpz(&ans.x0.x1,&A->x0.x2,Matrix_CVMA_to_karatsuba.data[2][1].x0);
    Fp2_add(&ans.x0.x1,&ans.x0.x1,&tmp0);
    Fp2_add(&ans.x0.x1,&ans.x0.x1,&tmp1);
    Fp2_mul_mpz(&tmp0,&A->x0.x0,Matrix_CVMA_to_karatsuba.data[0][2].x0);
    Fp2_mul_mpz(&tmp1,&A->x0.x1,Matrix_CVMA_to_karatsuba.data[1][2].x0);
    Fp2_mul_mpz(&ans.x0.x2,&A->x0.x2,Matrix_CVMA_to_karatsuba.data[2][2].x0);
    Fp2_add(&ans.x0.x2,&ans.x0.x2,&tmp0);
    Fp2_add(&ans.x0.x2,&ans.x0.x2,&tmp1);
    Fp2_mul_mpz(&tmp0,&A->x1.x0,Matrix_CVMA_to_karatsuba.data[0][0].x0);
    Fp2_mul_mpz(&tmp1,&A->x1.x1,Matrix_CVMA_to_karatsuba.data[1][0].x0);
    Fp2_mul_mpz(&ans.x1.x0,&A->x1.x2,Matrix_CVMA_to_karatsuba.data[2][0].x0);
    Fp2_add(&ans.x1.x0,&ans.x1.x0,&tmp0);
    Fp2_add(&ans.x1.x0,&ans.x1.x0,&tmp1);
    Fp2_mul_mpz(&tmp0,&A->x1.x0,Matrix_CVMA_to_karatsuba.data[0][1].x0);
    Fp2_mul_mpz(&tmp1,&A->x1.x1,Matrix_CVMA_to_karatsuba.data[1][1].x0);
    Fp2_mul_mpz(&ans.x1.x1,&A->x1.x2,Matrix_CVMA_to_karatsuba.data[2][1].x0);
    Fp2_add(&ans.x1.x1,&ans.x1.x1,&tmp0);
    Fp2_add(&ans.x1.x1,&ans.x1.x1,&tmp1);
    Fp2_mul_mpz(&tmp0,&A->x1.x0,Matrix_CVMA_to_karatsuba.data[0][2].x0);
    Fp2_mul_mpz(&tmp1,&A->x1.x1,Matrix_CVMA_to_karatsuba.data[1][2].x0);
    Fp2_mul_mpz(&ans.x1.x2,&A->x1.x2,Matrix_CVMA_to_karatsuba.data[2][2].x0);
    Fp2_add(&ans.x1.x2,&ans.x1.x2,&tmp0);
    Fp2_add(&ans.x1.x2,&ans.x1.x2,&tmp1);
    Fp12_set(ANS,&ans);
    
    Fp2_clear(&tmp0);
    Fp2_clear(&tmp1);
    Fp12_clear(&ans);
    
    gettimeofday(&tv_end,NULL);
    BASIS_CONVERSION_CTOK=timedifference_msec(tv_start,tv_end);
}
/*----------------------------------------------------------------------------*/
//twist
void EFp12_to_EFp2_karatsuba(EFp2 *ANS,EFp12 *A){
    Fp2_inv_basis(&ANS->x,&A->x.x0.x2);
    Fp2_inv_basis(&ANS->y,&A->y.x1.x0);
    Fp2_mul_mpz(&ANS->y,&ANS->y,Two_inv);
    ANS->infinity=A->infinity;
}

void EFp2_to_EFp12_karatsuba(EFp12 *ANS,EFp2 *A){
    Fp12_set_ui(&ANS->x,0);
    Fp12_set_ui(&ANS->y,0);
    Fp2_mul_basis(&ANS->x.x0.x2,&A->x);
    Fp2_mul_basis(&ANS->y.x1.x0,&A->y);
    Fp2_add(&ANS->y.x1.x0,&ANS->y.x1.x0,&ANS->y.x1.x0);
    ANS->infinity=A->infinity;
}

void EFp12_to_EFp_karatsuba(EFp *ANS,EFp12 *A){
    Fp_set(&ANS->x,&A->x.x0.x0.x0);
    Fp_set(&ANS->y,&A->y.x0.x0.x0);
    ANS->infinity=A->infinity;
}

void EFp_to_EFp12_karatsuba(EFp12 *ANS,EFp *A){
    Fp12_set_ui(&ANS->x,0);
    Fp12_set_ui(&ANS->y,0);
    Fp_set(&ANS->x.x0.x0.x0,&A->x);
    Fp_set(&ANS->y.x0.x0.x0,&A->y);
    ANS->infinity=A->infinity;
}

void EFp12_to_EFp_CVMA(EFp *ANS,EFp12 *A){
    Fp_set_neg(&ANS->x,&A->x.x0.x0.x0);
    Fp_set_neg(&ANS->y,&A->y.x0.x0.x0);
    ANS->infinity=A->infinity;
}

void EFp_to_EFp12_CVMA(EFp12 *ANS,EFp *A){
    Fp12_set_ui(&ANS->x,0);
    Fp12_set_ui(&ANS->y,0);
    Fp_set_neg(&ANS->x.x0.x0.x0,&A->x);
    Fp_set(&ANS->x.x0.x1.x0,&ANS->x.x0.x0.x0);
    Fp_set(&ANS->x.x0.x2.x0,&ANS->x.x0.x0.x0);
    Fp_set_neg(&ANS->y.x0.x0.x0,&A->y);
    Fp_set(&ANS->y.x0.x1.x0,&ANS->y.x0.x0.x0);
    Fp_set(&ANS->y.x0.x2.x0,&ANS->y.x0.x0.x0);
    ANS->infinity=A->infinity;
}

/*----------------------------------------------------------------------------*/
//Pseudo 8-sparse
void Pseudo_8_sparse_mapping(EFp *P,EFp2 *Q,Fp *L){	
	EFp_set(&TMP1_EFP,P);
	EFp2_set(&TMP1_EFP2,Q);
	Fp_mul(&TMP1_FP,&TMP1_EFP.x,&TMP1_EFP.y);
	Fp_inv(&TMP1_FP,&TMP1_FP);
	Fp_mul(&TMP2_FP,&TMP1_EFP.x,&TMP1_EFP.x);
	Fp_mul(&TMP2_FP,&TMP2_FP,&TMP1_FP);
	Fp_mul(&TMP3_FP,&TMP1_EFP.y,&TMP1_FP);
	Fp_mul(&TMP4_FP,&TMP2_FP,&TMP2_FP);
	Fp2_mul_mpz(&Q->x,&TMP1_EFP2.x,TMP4_FP.x0);
	Fp_mul(&TMP5_FP,&TMP2_FP,&TMP4_FP);
	Fp2_mul_mpz(&Q->y,&TMP1_EFP2.y,TMP5_FP.x0);
	Fp_mul(&P->x,&TMP4_FP,&TMP1_EFP.x);
	Fp_set(&P->y,&P->x);
	Fp_mul(L,&TMP3_FP,&TMP1_EFP.y);
	Fp_mul(L,L,L);
	Fp_mul(L,L,&TMP3_FP);
}

void Pseudo_8_sparse_mul_karatsuba(Fp12 *ANS,Fp12 *A,Fp12 *B){
	Fp2_mul(&TMP1_FP2,&A->x0.x0,&B->x1.x0);
	Fp2_mul(&TMP2_FP2,&A->x0.x1,&B->x1.x1);
	Fp2_add(&TMP3_FP2,&A->x0.x0,&A->x0.x1);
	Fp2_add(&TMP4_FP2,&B->x1.x0,&B->x1.x1);
	Fp2_mul(&TMP3_FP2,&TMP3_FP2,&TMP4_FP2);
	Fp2_sub(&TMP3_FP2,&TMP3_FP2,&TMP1_FP2);
	Fp2_sub(&TMP3_FP2,&TMP3_FP2,&TMP2_FP2);
	Fp2_add(&TMP1_FP12.x1.x0,&TMP1_FP2,&A->x1.x0);
	Fp2_add(&TMP1_FP12.x1.x1,&TMP3_FP2,&A->x1.x1);
	Fp2_add(&TMP1_FP12.x1.x2,&TMP2_FP2,&A->x1.x2);
	Fp2_mul(&TMP1_FP2,&A->x1.x0,&B->x1.x0);
	Fp2_mul(&TMP2_FP2,&A->x1.x1,&B->x1.x1);
	Fp2_add(&TMP3_FP2,&A->x1.x0,&A->x1.x1);
	Fp2_add(&TMP4_FP2,&B->x1.x0,&B->x1.x1);
	Fp2_mul(&TMP3_FP2,&TMP3_FP2,&TMP4_FP2);
	Fp2_sub(&TMP3_FP2,&TMP3_FP2,&TMP1_FP2);
	Fp2_sub(&TMP3_FP2,&TMP3_FP2,&TMP2_FP2);
	Fp2_mul(&TMP4_FP2,&A->x1.x2,&B->x1.x1);
	Fp2_add(&TMP4_FP2,&TMP4_FP2,&TMP4_FP2);
	Fp2_add(&TMP4_FP2,&TMP1_FP2,&TMP4_FP2);
	Fp2_mul_basis(&TMP1_FP12.x0.x0,&TMP4_FP2);
	Fp2_add(&TMP1_FP12.x0.x0,&TMP1_FP12.x0.x0,&A->x0.x0);
	Fp2_mul_basis(&TMP1_FP12.x0.x1,&TMP3_FP2);
	Fp2_add(&TMP1_FP12.x0.x1,&TMP1_FP12.x0.x1,&A->x0.x1);
	Fp2_mul(&TMP4_FP2,&A->x1.x2,&B->x1.x0);
	Fp2_add(&TMP4_FP2,&TMP2_FP2,&TMP4_FP2);
	Fp2_mul_basis(&TMP1_FP12.x0.x2,&TMP4_FP2);
	Fp2_add(&TMP1_FP12.x0.x2,&TMP1_FP12.x0.x2,&A->x0.x2);
	Fp2_mul(&TMP1_FP2,&A->x0.x2,&B->x1.x1);
	Fp2_add(&TMP1_FP2,&TMP1_FP2,&TMP1_FP2);
	Fp2_add(&TMP1_FP12.x1.x0,&TMP1_FP12.x1.x0,&TMP1_FP2);
	Fp2_mul(&TMP1_FP2,&A->x0.x2,&B->x1.x0);
	Fp2_add(&TMP1_FP12.x1.x2,&TMP1_FP12.x1.x2,&TMP1_FP2);
	Fp12_set(ANS,&TMP1_FP12);
}

void ff_ltt_karatsuba(Fp12 *f,EFp2 *T,EFp *P,Fp *L){
	EFp2_set(&TMP1_EFP2,T);
	Fp2_add(&TMP1_FP2,&TMP1_EFP2.y,&TMP1_EFP2.y);
	Fp2_inv(&TMP1_FP2,&TMP1_FP2);
	Fp2_sqr(&TMP2_FP2,&TMP1_EFP2.x);
	Fp2_add(&TMP3_FP2,&TMP2_FP2,&TMP2_FP2);
	Fp2_add(&TMP2_FP2,&TMP2_FP2,&TMP3_FP2);
	Fp2_mul(&TMP3_FP2,&TMP1_FP2,&TMP2_FP2);
	Fp2_add(&TMP4_FP2,&TMP1_EFP2.x,&TMP1_EFP2.x);
	Fp2_sqr(&T->x,&TMP3_FP2);
	Fp2_sub(&T->x,&T->x,&TMP4_FP2);
	Fp2_mul(&TMP5_FP2,&TMP3_FP2,&TMP1_EFP2.x);
	Fp2_sub(&TMP5_FP2,&TMP5_FP2,&TMP1_EFP2.y);
	Fp2_mul(&T->y,&TMP3_FP2,&T->x);
	Fp2_sub(&T->y,&TMP5_FP2,&T->y);
    
	Fp_set_ui(&TMP2_FP12.x0.x0.x0,1);
	Fp2_set_neg(&TMP2_FP12.x1.x1,&TMP3_FP2);
	Fp2_mul_mpz(&TMP2_FP12.x1.x0,&TMP5_FP2,L->x0);
    Fp2_mul_basis(&TMP2_FP12.x1.x0,&TMP2_FP12.x1.x0);
    Fp2_add(&TMP2_FP12.x1.x0,&TMP2_FP12.x1.x0,&TMP2_FP12.x1.x0);
    
	Fp12_sqr_karatsuba(f,f);
	Pseudo_8_sparse_mul_karatsuba(f,f,&TMP2_FP12);
}

void f_ltq_karatsuba(Fp12 *f,EFp2 *T,EFp2 *Q,EFp *P,Fp *L){
	EFp2_set(&TMP1_EFP2,T);
	Fp2_sub(&TMP1_FP2,&Q->x,&TMP1_EFP2.x);
	Fp2_inv(&TMP1_FP2,&TMP1_FP2);
	Fp2_sub(&TMP2_FP2,&Q->y,&TMP1_EFP2.y);
	Fp2_mul(&TMP3_FP2,&TMP1_FP2,&TMP2_FP2);
	Fp2_add(&TMP4_FP2,&TMP1_EFP2.x,&Q->x);
	Fp2_sqr(&T->x,&TMP3_FP2);
	Fp2_sub(&T->x,&T->x,&TMP4_FP2);
	Fp2_mul(&TMP5_FP2,&TMP3_FP2,&TMP1_EFP2.x);
	Fp2_sub(&TMP5_FP2,&TMP5_FP2,&TMP1_EFP2.y);
	Fp2_mul(&T->y,&TMP3_FP2,&T->x);
	Fp2_sub(&T->y,&TMP5_FP2,&T->y);
	
	Fp_set_ui(&TMP2_FP12.x0.x0.x0,1);
	Fp2_set_neg(&TMP2_FP12.x1.x1,&TMP3_FP2);
	Fp2_mul_mpz(&TMP2_FP12.x1.x0,&TMP5_FP2,L->x0);
    Fp2_mul_basis(&TMP2_FP12.x1.x0,&TMP2_FP12.x1.x0);
    Fp2_add(&TMP2_FP12.x1.x0,&TMP2_FP12.x1.x0,&TMP2_FP12.x1.x0);
	
	Pseudo_8_sparse_mul_karatsuba(f,f,&TMP2_FP12);
}

/*----------------------------------------------------------------------------*/
//miller
void Miller_algo_for_plain_ate_karatsuba(Fp12 *ANS,EFp12 *Q,EFp12 *P){
    gettimeofday(&tv_start,NULL);
    
    EFp12 Test;
	EFp12_init(&Test);
	EFp2 T;
	EFp2_init(&T);
	EFp2 mapped_Q;
	EFp2_init(&mapped_Q);
	EFp mapped_P;
	EFp_init(&mapped_P);
	Fp12 f;
	Fp12_init(&f);
	Fp L;
	Fp_init(&L);
	mpz_t loop;
	mpz_init(loop);
	mpz_sub_ui(loop,trace,1);
	int i,length;
	length=(int)mpz_sizeinbase(loop,2);
	char binary[length];
	mpz_get_str(binary,2,loop);
	
	EFp12_to_EFp_karatsuba(&mapped_P,P);
	EFp12_to_EFp2_karatsuba(&mapped_Q,Q);
	Pseudo_8_sparse_mapping(&mapped_P,&mapped_Q,&L);
	EFp2_set(&T,&mapped_Q);
	Fp_set_ui(&f.x0.x0.x0,1);
    
    for(i=1; binary[i]!='\0'; i++){
		ff_ltt_karatsuba(&f,&T,&mapped_P,&L);
		if(binary[i]=='1'){
			f_ltq_karatsuba(&f,&T,&mapped_Q,&mapped_P,&L);
		}
	}
	Fp12_set(ANS,&f);
	
	Fp12_clear(&f);
	EFp2_clear(&T);
	EFp2_clear(&mapped_Q);
	EFp_clear(&mapped_P);
	Fp_clear(&L);
	EFp12_clear(&Test);
	mpz_clear(loop);
    
    gettimeofday(&tv_end,NULL);
    MILLER_PLAINATE=timedifference_msec(tv_start,tv_end);
}

void Miller_algo_for_opt_ate_karatsuba(Fp12 *ANS,EFp12 *Q,EFp12 *P){
    gettimeofday(&tv_start,NULL);
    
    EFp12 Buf;
    EFp12_init(&Buf);
    EFp2 T;
    EFp2_init(&T);
    EFp2 mapped_Q,mapped_Q_neg,mapped_Q1,mapped_Q2_neg;
    EFp2_init(&mapped_Q);
    EFp2_init(&mapped_Q_neg);
    EFp2_init(&mapped_Q1);
    EFp2_init(&mapped_Q2_neg);
    EFp mapped_P;
    EFp_init(&mapped_P);
    Fp12 f;
    Fp12_init(&f);
    Fp L;
    Fp_init(&L);
    int i;
    
    EFp12_to_EFp_karatsuba(&mapped_P,P);
    EFp12_to_EFp2_karatsuba(&mapped_Q,Q);
    Pseudo_8_sparse_mapping(&mapped_P,&mapped_Q,&L);
    EFp2_set(&mapped_Q_neg,&mapped_Q);
    Fp2_set_neg(&mapped_Q_neg.y,&mapped_Q_neg.y);
    EFp2_set(&T,&mapped_Q);
    Fp_set_ui(&f.x0.x0.x0,1);
    
    for(i=X6_2_length-1; i>=0; i--){
        switch(X6_2_binary[i]){
            case 0:
                ff_ltt_karatsuba(&f,&T,&mapped_P,&L);
                break;
            case 1:
                ff_ltt_karatsuba(&f,&T,&mapped_P,&L);
                f_ltq_karatsuba(&f,&T,&mapped_Q,&mapped_P,&L);
                break;
            case -1:
                ff_ltt_karatsuba(&f,&T,&mapped_P,&L);
                f_ltq_karatsuba(&f,&T,&mapped_Q_neg,&mapped_P,&L);
                break;
            default:
                break;
        }
        
    }
    
    EFp2_skew_frobenius_map_p1_karatsuba(&mapped_Q1,&mapped_Q);
    EFp2_skew_frobenius_map_p2_karatsuba(&mapped_Q2_neg,&mapped_Q);
    EFp2_set_neg(&mapped_Q2_neg,&mapped_Q2_neg);
    f_ltq_karatsuba(&f,&T,&mapped_Q1,&mapped_P,&L);
    f_ltq_karatsuba(&f,&T,&mapped_Q2_neg,&mapped_P,&L);
    Fp12_set(ANS,&f);
    
    EFp12_clear(&Buf);
    Fp12_clear(&f);
    EFp2_clear(&T);
    EFp2_clear(&mapped_Q);
    EFp2_clear(&mapped_Q_neg);
    EFp2_clear(&mapped_Q1);
    EFp2_clear(&mapped_Q2_neg);
    EFp_clear(&mapped_P);
    Fp_clear(&L);
    
    gettimeofday(&tv_end,NULL);
    MILLER_OPTATE=timedifference_msec(tv_start,tv_end);
}

void Miller_algo_for_x_ate_karatsuba(Fp12 *ANS,EFp12 *Q,EFp12 *P){
    gettimeofday(&tv_start,NULL);
    
    EFp12 Buf;
    EFp12_init(&Buf);
    EFp2 T;
    EFp2_init(&T);
    EFp2 mapped_Q,mapped_Q_neg,mapped_Q1,mapped_Q2;
    EFp2_init(&mapped_Q);
    EFp2_init(&mapped_Q_neg);
    EFp2_init(&mapped_Q1);
    EFp2_init(&mapped_Q2);
    EFp mapped_P;
    EFp_init(&mapped_P);
    Fp12 f;
    Fp12_init(&f);
    Fp L;
    Fp_init(&L);
    int i;
    
    EFp12_to_EFp_karatsuba(&mapped_P,P);
    EFp12_to_EFp2_karatsuba(&mapped_Q,Q);
    Pseudo_8_sparse_mapping(&mapped_P,&mapped_Q,&L);
    EFp2_set(&mapped_Q_neg,&mapped_Q);
    Fp2_set_neg(&mapped_Q_neg.y,&mapped_Q_neg.y);
    EFp2_set(&T,&mapped_Q);
    Fp_set_ui(&f.x0.x0.x0,1);
    
    for(i=X_length-1; i>=0; i--){
        switch(X_binary[i]){
            case 0:
                ff_ltt_karatsuba(&f,&T,&mapped_P,&L);
                break;
            case 1:
                ff_ltt_karatsuba(&f,&T,&mapped_P,&L);
                f_ltq_karatsuba(&f,&T,&mapped_Q,&mapped_P,&L);
                break;
            case -1:
                ff_ltt_karatsuba(&f,&T,&mapped_P,&L);
                f_ltq_karatsuba(&f,&T,&mapped_Q_neg,&mapped_P,&L);
                break;
            default:
                break;
        }
    }
    
    Fp12_frobenius_map_p3_karatsuba(&Buf.x,&f);
    Fp12_mul_karatsuba(&f,&Buf.x,&f);
    EFp2_skew_frobenius_map_p3_karatsuba(&mapped_Q1,&T);
    EFp2_set(&mapped_Q2,&T);
    f_ltq_karatsuba(&f,&mapped_Q2,&mapped_Q1,&mapped_P,&L);
    Fp12_frobenius_map_p10_karatsuba(&Buf.x,&f);
    Fp12_mul_karatsuba(&f,&Buf.x,&f);
    EFp2_skew_frobenius_map_p10_karatsuba(&T,&mapped_Q2);
    f_ltq_karatsuba(&f,&T,&mapped_Q2,&mapped_P,&L);
    Fp12_set(ANS,&f);
    
    EFp12_clear(&Buf);
    Fp12_clear(&f);
    EFp2_clear(&T);
    EFp2_clear(&mapped_Q);
    EFp2_clear(&mapped_Q_neg);
    EFp2_clear(&mapped_Q1);
    EFp2_clear(&mapped_Q2);
    EFp_clear(&mapped_P);
    Fp_clear(&L);
    
    gettimeofday(&tv_end,NULL);
    MILLER_XATE=timedifference_msec(tv_start,tv_end);
}

/*----------------------------------------------------------------------------*/
//final exp
//CVMA
void Fp12_pow_X_CVMA(Fp12 *ANS,Fp12 *A){
    int i;
    Fp12 tmp,A_inv;
	Fp12_init(&tmp);
	Fp12_init(&A_inv);
	Fp12_frobenius_map_p6_CVMA(&A_inv,A);
	
	Fp12_set(&tmp,A);
	for(i=X_length-1; i>=0; i--){
		switch(X_binary[i]){
			case 0:
				Fp12_sqr_cyclotomic_CVMA(&tmp,&tmp);
				break;
			case 1:
				Fp12_sqr_cyclotomic_CVMA(&tmp,&tmp);
				Fp12_mul_CVMA(&tmp,&tmp,A);
				break;
			case -1:
				Fp12_sqr_cyclotomic_CVMA(&tmp,&tmp);
				Fp12_mul_CVMA(&tmp,&tmp,&A_inv);
				break;
			default:
				break;
		}
	}
	Fp12_set(ANS,&tmp);
	
	Fp12_clear(&tmp);
	Fp12_clear(&A_inv);
}

void Final_exp_easy_CVMA(Fp12 *ANS,Fp12 *f){
    gettimeofday(&tv_start,NULL);
    
    Fp12_frobenius_map_p6_CVMA(&TMP1_FP12,f);
    Fp12_inv_CVMA(&TMP2_FP12,f);
    Fp12_mul_CVMA(f,&TMP1_FP12,&TMP2_FP12);
    Fp12_frobenius_map_p2_CVMA(&TMP1_FP12,f);
    Fp12_mul_CVMA(ANS,&TMP1_FP12,f);
    
    gettimeofday(&tv_end,NULL);
    FINALEXP_EASY=timedifference_msec(tv_start,tv_end);
}

void Final_exp_hard_plain_CVMA(Fp12 *ANS,Fp12 *f){
    gettimeofday(&tv_start,NULL);
    
    mpz_t exp,buf;
    mpz_init(exp);
    mpz_init(buf);
    
    mpz_pow_ui(exp,prime,4);
    mpz_pow_ui(buf,prime,2);
    mpz_sub(exp,exp,buf);
    mpz_add_ui(exp,exp,1);
    mpz_tdiv_q(exp,exp,order);
    Fp12_pow_CVMA(ANS,f,exp);
    
    mpz_clear(exp);
    mpz_clear(buf);
    
    gettimeofday(&tv_end,NULL);
    FINALEXP_HARD_PLAIN=timedifference_msec(tv_start,tv_end);
}

void Final_exp_hard_optimal_CVMA(Fp12 *ANS,Fp12 *f){
    gettimeofday(&tv_start,NULL);
    
    Fp12_pow_X_CVMA(&TMP1_FP12,f);
    Fp12_frobenius_map_p6_CVMA(&TMP1_FP12,&TMP1_FP12);
    Fp12_sqr_cyclotomic_CVMA(&TMP1_FP12,&TMP1_FP12);
    Fp12_sqr_cyclotomic_CVMA(&TMP2_FP12,&TMP1_FP12);
    Fp12_mul_CVMA(&TMP2_FP12,&TMP1_FP12,&TMP2_FP12);
    Fp12_pow_X_CVMA(&TMP3_FP12,&TMP2_FP12);
    Fp12_frobenius_map_p6_CVMA(&TMP3_FP12,&TMP3_FP12);
    Fp12_frobenius_map_p6_CVMA(&TMP4_FP12,&TMP2_FP12);
    Fp12_mul_CVMA(&TMP2_FP12,&TMP3_FP12,&TMP4_FP12);
    Fp12_sqr_cyclotomic_CVMA(&TMP4_FP12,&TMP3_FP12);
    Fp12_pow_X_CVMA(&TMP5_FP12,&TMP4_FP12);
    Fp12_frobenius_map_p6_CVMA(&TMP5_FP12,&TMP5_FP12);
    Fp12_frobenius_map_p6_CVMA(&TMP5_FP12,&TMP5_FP12);
    Fp12_mul_CVMA(&TMP5_FP12,&TMP5_FP12,&TMP2_FP12);
    Fp12_mul_CVMA(&TMP4_FP12,&TMP5_FP12,&TMP1_FP12);
    Fp12_mul_CVMA(&TMP1_FP12,&TMP3_FP12,&TMP5_FP12);
    Fp12_mul_CVMA(&TMP1_FP12,&TMP1_FP12,f);
    Fp12_frobenius_map_p1_CVMA(&TMP3_FP12,&TMP4_FP12);
    Fp12_mul_CVMA(&TMP1_FP12,&TMP3_FP12,&TMP1_FP12);
    Fp12_frobenius_map_p2_CVMA(&TMP3_FP12,&TMP5_FP12);
    Fp12_mul_CVMA(&TMP1_FP12,&TMP3_FP12,&TMP1_FP12);
    Fp12_frobenius_map_p6_CVMA(&TMP3_FP12,f);
    Fp12_mul_CVMA(&TMP3_FP12,&TMP3_FP12,&TMP4_FP12);
    Fp12_frobenius_map_p3_CVMA(&TMP3_FP12,&TMP3_FP12);
    Fp12_mul_CVMA(ANS,&TMP3_FP12,&TMP1_FP12);
    
    gettimeofday(&tv_end,NULL);
    FINALEXP_HARD_OPT=timedifference_msec(tv_start,tv_end);
}

/*----------------------------------------------------------------------------*/
//pairing
void Plain_ate_pairing(Fp12 *ANS,EFp12 *P,EFp12 *Q){
    Miller_algo_for_plain_ate_karatsuba(ANS,P,Q);
    karatsuba_to_CVMA(ANS,ANS);
    Final_exp_easy_CVMA(ANS,ANS);
    Final_exp_hard_plain_CVMA(ANS,ANS);
}

void Opt_ate_pairing(Fp12 *ANS,EFp12 *P,EFp12 *Q){
    Miller_algo_for_opt_ate_karatsuba(ANS,P,Q);
    karatsuba_to_CVMA(ANS,ANS);
    Final_exp_easy_CVMA(ANS,ANS);
    Final_exp_hard_optimal_CVMA(ANS,ANS);
}

void X_ate_pairing(Fp12 *ANS,EFp12 *P,EFp12 *Q){
    Miller_algo_for_x_ate_karatsuba(ANS,P,Q);
    karatsuba_to_CVMA(ANS,ANS);
    Final_exp_easy_CVMA(ANS,ANS);
    Final_exp_hard_optimal_CVMA(ANS,ANS);
}

/*----------------------------------------------------------------------------*/
//JSF
void Joint_sparse_form(int **binary,mpz_t scalar[2],int *loop_length){
	int i,j;
	unsigned long int u;
	mpz_t mod_2,mod_4,mod_8;
	mpz_init(mod_2);
	mpz_init(mod_4);
	mpz_init(mod_8);
	mpz_t k[2];
	mpz_init(k[0]);
	mpz_init(k[1]);
    
	j=0;
	mpz_set(k[0],scalar[0]);
	mpz_set(k[1],scalar[1]);
	
	while(mpz_cmp_ui(k[0],0)>0 || mpz_cmp_ui(k[1],0)>0){
		for(i=0; i<2; i++){
			mpz_mod_ui(mod_2,k[i],2);
			if(mpz_cmp_ui(mod_2,0)==0){
				u=0;
			}else{
				mpz_mod_ui(mod_4,k[i],4);
				u=mpz_get_ui(mod_4);
				if(u==3){
					u=-1;
				}
				mpz_mod_ui(mod_8,k[i],8);
				mpz_mod_ui(mod_4,k[1-i],4);
				if((mpz_cmp_ui(mod_8,3)==0 || mpz_cmp_ui(mod_8,5)==0) && mpz_cmp_ui(mod_4,2)==0){
					u=-u;
				}
			}
			binary[i][j]=u;
		}
		for(i=0; i<2; i++){
			u=binary[i][j];
			switch (u){
				case 1:
					mpz_sub_ui(k[i],k[i],1);
					break;
				case -1:
					mpz_add_ui(k[i],k[i],1);
					break;
				default:
					break;
			}
			mpz_tdiv_q_ui(k[i],k[i],2);
		}
		j=j+1;
	}
	*loop_length=j-1;
	
	mpz_clear(mod_2);
	mpz_clear(mod_4);
	mpz_clear(mod_8);
	mpz_clear(k[0]);
	mpz_clear(k[1]);	
}

/*----------------------------------------------------------------------------*/
//G1 SCM
void EFp12_G1_SCM_plain_karatsuba(EFp12 *ANS,EFp12 *P,mpz_t scalar){
    gettimeofday(&tv_start,NULL);
    
    EFp tmp_P;
	EFp_init(&tmp_P);
	
	EFp12_to_EFp_karatsuba(&tmp_P,P);
	EFp_SCM(&tmp_P,&tmp_P,scalar);
	EFp_to_EFp12_karatsuba(ANS,&tmp_P);
	
	EFp_clear(&tmp_P);
	
    gettimeofday(&tv_end,NULL);
    G1SCM_PLAIN=timedifference_msec(tv_start,tv_end);
}

void EFp12_G1_SCM_2split_karatsuba(EFp12 *ANS,EFp12 *P,mpz_t scalar){
    gettimeofday(&tv_start,NULL);
    
    int i,length_s[2],loop_length;
	EFp next_tmp_P,tmp_P,skew_P;
	EFp_init(&next_tmp_P);
	EFp_init(&tmp_P);
	EFp_init(&skew_P);
	mpz_t s[2],buf,V1,V2,s1,s2,s3,s4,s5,CHECK;
	mpz_init(buf);
	mpz_init(V1);
	mpz_init(V2);
	mpz_init(s1);
	mpz_init(s2);
	mpz_init(s3);
	mpz_init(s4);
	mpz_init(s5);
	mpz_init(CHECK);
	for(i=0; i<2; i++){
		mpz_init(s[i]);
	}
	EFp table[4];
	for(i=0; i<4; i++){
		EFp_init(&table[i]);
	}
	
	mpz_mul(V1,X,X);
	mpz_mul_ui(V1,V1,6);
	mpz_mul_ui(buf,X,4);
	mpz_add(V1,V1,buf);
	mpz_add_ui(V1,V1,1);
	mpz_add(V2,X,X);
	mpz_add_ui(V2,V2,1);
	mpz_tdiv_qr(s1,s2,scalar,V1);
	mpz_mul(buf,V2,s1);
	mpz_tdiv_qr(s3,s4,buf,V1);
	mpz_mul(s5,V2,s3);
	mpz_add(s[1],s4,s5);
	mpz_mod(s[1],s[1],order);
	mpz_sub(s[0],s2,s5);
	mpz_mod(s[0],s[0],order);
	mpz_sub_ui(CHECK,order,1);
	mpz_tdiv_q_ui(CHECK,CHECK,2);
	
	EFp12_to_EFp_karatsuba(&tmp_P,P);
	EFp_skew_frobenius_map_p2(&skew_P,&tmp_P);
	if(mpz_cmp(s[0],CHECK)>0){
		mpz_sub(s[0],order,s[0]);
		EFp_set_neg(&tmp_P,&tmp_P);
	}
	table[0].infinity=1;                //00
	EFp_set(&table[1],&tmp_P);          //01
	EFp_set(&table[2],&skew_P);         //10
	EFp_ECA(&table[3],&tmp_P,&skew_P);  //11
    
	loop_length=0;
	for(i=0; i<2; i++){
		length_s[i]=(int)mpz_sizeinbase(s[i],2);
		if(loop_length<length_s[i]){
			loop_length=length_s[i];
		}
	}
	char binary_s[2][loop_length+1];
	char str[5],*e;
	int binary[loop_length+1];
	for(i=0; i<2; i++){
		if(length_s[i]==loop_length){
			mpz_get_str(binary_s[i],2,s[i]);
		}else{
			char binary_buf[loop_length+1];
			mpz_get_str(binary_buf,2,s[i]);
			memset(binary_s[i],'0',sizeof(binary_s[i]));
			memmove(binary_s[i]+loop_length-length_s[i],binary_buf,sizeof(binary_buf));
		}
	}
	for(i=0; i<loop_length; i++){
		sprintf(str,"%c%c",binary_s[1][i],binary_s[0][i]);
		binary[i]=(int)strtol(str,&e,2);
	}
	EFp_set(&next_tmp_P,&table[binary[0]]);
	
	for(i=1; i<loop_length; i++){
		EFp_ECD(&next_tmp_P,&next_tmp_P);
		EFp_ECA(&next_tmp_P,&next_tmp_P,&table[binary[i]]);
	}
	EFp_to_EFp12_karatsuba(ANS,&next_tmp_P);
	
	mpz_clear(buf);
	mpz_clear(V1);
	mpz_clear(V2);
	mpz_clear(s1);
	mpz_clear(s2);
	mpz_clear(s3);
	mpz_clear(s4);
	mpz_clear(s5);
	mpz_clear(CHECK);
	EFp_clear(&next_tmp_P);
	EFp_clear(&tmp_P);
	EFp_clear(&skew_P);
	for(i=0; i<2; i++){
		mpz_clear(s[i]);
	}
	for(i=0; i<4; i++){
		EFp_clear(&table[i]);
	}
	
    gettimeofday(&tv_end,NULL);
    G1SCM_2SPLIT=timedifference_msec(tv_start,tv_end);
}

void EFp12_G1_SCM_2split_JSF_karatsuba(EFp12 *ANS,EFp12 *P,mpz_t scalar){
    gettimeofday(&tv_start,NULL);
    
    int i,length_s[2],loop_length;
	EFp next_tmp_P,tmp_P,tmp_P_neg,skew_P,skew_P_neg;
	EFp_init(&next_tmp_P);
	EFp_init(&tmp_P);
	EFp_init(&tmp_P_neg);
	EFp_init(&skew_P);
	EFp_init(&skew_P_neg);
	mpz_t s[2],buf,V1,V2,s1,s2,s3,s4,s5,CHECK;
	mpz_init(buf);
	mpz_init(V1);
	mpz_init(V2);
	mpz_init(s1);
	mpz_init(s2);
	mpz_init(s3);
	mpz_init(s4);
	mpz_init(s5);
	mpz_init(CHECK);
	for(i=0; i<2; i++){
		mpz_init(s[i]);
	}
	EFp table[9];
	for(i=0; i<9; i++){
		EFp_init(&table[i]);
	}
	
	mpz_mul(V1,X,X);
	mpz_mul_ui(V1,V1,6);
	mpz_mul_ui(buf,X,4);
	mpz_add(V1,V1,buf);
	mpz_add_ui(V1,V1,1);
	mpz_add(V2,X,X);
	mpz_add_ui(V2,V2,1);
	mpz_tdiv_qr(s1,s2,scalar,V1);
	mpz_mul(buf,V2,s1);
	mpz_tdiv_qr(s3,s4,buf,V1);
	mpz_mul(s5,V2,s3);
    mpz_add(s[1],s4,s5);
	mpz_mod(s[1],s[1],order);
	mpz_sub(s[0],s2,s5);
	mpz_mod(s[0],s[0],order);
	mpz_sub_ui(CHECK,order,1);
	mpz_tdiv_q_ui(CHECK,CHECK,2);
	
	EFp12_to_EFp_karatsuba(&tmp_P,P);
	EFp_set_neg(&tmp_P_neg,&tmp_P);
	EFp_skew_frobenius_map_p2(&skew_P,&tmp_P);
	EFp_set_neg(&skew_P_neg,&skew_P);
	
	if(mpz_cmp(s[0],CHECK)>0){
		mpz_sub(s[0],order,s[0]);
		EFp_set_neg(&tmp_P,&tmp_P);
		EFp_set_neg(&tmp_P_neg,&tmp_P_neg);
	}
	
	table[0].infinity=1;                        //00
	EFp_set(&table[1],&tmp_P);                  //01
	EFp_set(&table[2],&skew_P);                 //10
	EFp_ECA(&table[3],&skew_P,&tmp_P);          //11
	EFp_set(&table[4],&tmp_P_neg);              //0-1
	EFp_set(&table[5],&skew_P_neg);             //-10
	EFp_ECA(&table[6],&skew_P_neg,&tmp_P_neg);  //-1-1
	EFp_ECA(&table[7],&skew_P,&tmp_P_neg);      //1-1
	EFp_ECA(&table[8],&skew_P_neg,&tmp_P);      //-11
	
	loop_length=0;
	for(i=0; i<2; i++){
		length_s[i]=(int)mpz_sizeinbase(s[i],2);
		if(loop_length<length_s[i]){
			loop_length=length_s[i];
		}
	}
	int JSF_length;
	int JSF_binary[2][loop_length+1];
	char check[5];
	for(i=0; i<loop_length; i++){
		JSF_binary[0][i]=0;
		JSF_binary[1][i]=0;
	}
	int *JSF_pointer[2];
	JSF_pointer[0]=JSF_binary[0];
	JSF_pointer[1]=JSF_binary[1];
	Joint_sparse_form(JSF_pointer,s,&JSF_length);
	int binary[JSF_length+1];
	for(i=JSF_length; i>=0; i--){
		if(JSF_binary[1][i]==0 && JSF_binary[0][i]==0) 		binary[i]=0;
		else if(JSF_binary[1][i]==0 && JSF_binary[0][i]==1) 	binary[i]=1;
		else if(JSF_binary[1][i]==1 && JSF_binary[0][i]==0) 	binary[i]=2;
		else if(JSF_binary[1][i]==1 && JSF_binary[0][i]==1)	binary[i]=3;
		else if(JSF_binary[1][i]==0 && JSF_binary[0][i]==-1)	binary[i]=4;
		else if(JSF_binary[1][i]==-1 && JSF_binary[0][i]==0)	binary[i]=5;
		else if(JSF_binary[1][i]==-1 && JSF_binary[0][i]==-1)	binary[i]=6;
		else if(JSF_binary[1][i]==1 && JSF_binary[0][i]==-1)	binary[i]=7;
		else if(JSF_binary[1][i]==-1 && JSF_binary[0][i]==1)	binary[i]=8;
	}
	EFp_set(&next_tmp_P,&table[binary[JSF_length]]);
    
	for(i=JSF_length-1; i>=0; i--){
		EFp_ECD(&next_tmp_P,&next_tmp_P);
		EFp_ECA(&next_tmp_P,&next_tmp_P,&table[binary[i]]);
	}
	EFp_to_EFp12_karatsuba(ANS,&next_tmp_P);
	
	mpz_clear(buf);
	mpz_clear(V1);
	mpz_clear(V2);
	mpz_clear(s1);
	mpz_clear(s2);
	mpz_clear(s3);
	mpz_clear(s4);
	mpz_clear(s5);
	mpz_clear(CHECK);
	EFp_clear(&next_tmp_P);
	EFp_clear(&tmp_P);
	EFp_clear(&tmp_P_neg);
	EFp_clear(&skew_P);
	EFp_clear(&skew_P_neg);
	for(i=0; i<2; i++){
		mpz_clear(s[i]);
	}
	for(i=0; i<9; i++){
		EFp_clear(&table[i]);
	}
	
    gettimeofday(&tv_end,NULL);
    G1SCM_2SPLIT_JSF=timedifference_msec(tv_start,tv_end);
}

/*----------------------------------------------------------------------------*/
//G2 SCM
void EFp12_G2_SCM_plain_karatsuba(EFp12 *ANS,EFp12 *Q,mpz_t scalar){
    gettimeofday(&tv_start,NULL);
    
    EFp2 tmp_Q;
	EFp2_init(&tmp_Q);
	
	EFp12_to_EFp2_karatsuba(&tmp_Q,Q);
	EFp2_SCM(&tmp_Q,&tmp_Q,scalar);
	EFp2_to_EFp12_karatsuba(ANS,&tmp_Q);
	
	EFp2_clear(&tmp_Q);
	
    gettimeofday(&tv_end,NULL);
    G2SCM_PLAIN=timedifference_msec(tv_start,tv_end);
}

void EFp12_G2_SCM_2split_karatsuba(EFp12 *ANS,EFp12 *Q,mpz_t scalar){
    gettimeofday(&tv_start,NULL);
    
    int i,length_s[2],loop_length;
	EFp2 next_twisted_Q,twisted_Q,skew_Q;
	EFp2_init(&next_twisted_Q);
	EFp2_init(&twisted_Q);
	EFp2_init(&skew_Q);
	mpz_t s[2],buf;
	mpz_init(buf);
	for(i=0; i<2; i++){
		mpz_init(s[i]);
	}
	EFp2 table[4];
	for(i=0; i<4; i++){
		EFp2_init(&table[i]);
	}
	
    mpz_sub_ui(buf,trace,1);
    mpz_tdiv_qr(s[1],s[0],scalar,buf);
    
	EFp12_to_EFp2_karatsuba(&twisted_Q,Q);
	EFp2_skew_frobenius_map_p1_karatsuba(&skew_Q,&twisted_Q);
	table[0].infinity=1;                    //00
	EFp2_set(&table[1],&twisted_Q);         //01
	EFp2_set(&table[2],&skew_Q);            //10
	EFp2_ECA(&table[3],&twisted_Q,&skew_Q); //11
    
	loop_length=0;
	for(i=0; i<2; i++){
		length_s[i]=(int)mpz_sizeinbase(s[i],2);
		if(loop_length<length_s[i]){
			loop_length=length_s[i];
		}
	}
	char binary_s[2][loop_length+1];
	char str[5],*e;
	int binary[loop_length+1];
	for(i=0; i<2; i++){
		if(length_s[i]==loop_length){
			mpz_get_str(binary_s[i],2,s[i]);
		}else{
			char binary_buf[loop_length+1];
			mpz_get_str(binary_buf,2,s[i]);
			memset(binary_s[i],'0',sizeof(binary_s[i]));
			memmove(binary_s[i]+loop_length-length_s[i],binary_buf,sizeof(binary_buf));
		}
	}
	for(i=0; i<loop_length; i++){
		sprintf(str,"%c%c",binary_s[1][i],binary_s[0][i]);
		binary[i]=(int)strtol(str,&e,2);
	}
	EFp2_set(&next_twisted_Q,&table[binary[0]]);
	
	for(i=1; i<loop_length; i++){
		EFp2_ECD(&next_twisted_Q,&next_twisted_Q);
		EFp2_ECA(&next_twisted_Q,&next_twisted_Q,&table[binary[i]]);
	}
	EFp2_to_EFp12_karatsuba(ANS,&next_twisted_Q);
	
	mpz_clear(buf);
	EFp2_clear(&next_twisted_Q);
	EFp2_clear(&twisted_Q);
	EFp2_clear(&skew_Q);
	for(i=0; i<2; i++){
		mpz_clear(s[i]);
	}
	for(i=0; i<4; i++){
		EFp2_clear(&table[i]);
	}
	
    gettimeofday(&tv_end,NULL);
    G2SCM_2SPLIT=timedifference_msec(tv_start,tv_end);
}

void EFp12_G2_SCM_2split_JSF_karatsuba(EFp12 *ANS,EFp12 *Q,mpz_t scalar){
    gettimeofday(&tv_start,NULL);
    
    int i,length_s[2],loop_length;
	EFp2 next_tmp_Q,tmp_Q,tmp_Q_neg,skew_Q,skew_Q_neg;
	EFp2_init(&next_tmp_Q);
	EFp2_init(&tmp_Q);
	EFp2_init(&tmp_Q_neg);
	EFp2_init(&skew_Q);
	EFp2_init(&skew_Q_neg);
	mpz_t s[2],buf;
	mpz_init(buf);
	for(i=0; i<2; i++){
		mpz_init(s[i]);
	}
	EFp2 table[9];
	for(i=0; i<9; i++){
		EFp2_init(&table[i]);
	}
	
    mpz_sub_ui(buf,trace,1);
    mpz_tdiv_qr(s[1],s[0],scalar,buf);
    
	EFp12_to_EFp2_karatsuba(&tmp_Q,Q);
	EFp2_set_neg(&tmp_Q_neg,&tmp_Q);
	EFp2_skew_frobenius_map_p1_karatsuba(&skew_Q,&tmp_Q);
	EFp2_set_neg(&skew_Q_neg,&skew_Q);
	table[0].infinity=1;                        //00
	EFp2_set(&table[1],&tmp_Q);                 //01
	EFp2_set(&table[2],&skew_Q);                //10
	EFp2_ECA(&table[3],&skew_Q,&tmp_Q);         //11
	EFp2_set(&table[4],&tmp_Q_neg);             //0-1
    EFp2_set(&table[5],&skew_Q_neg);            //-10
	EFp2_ECA(&table[6],&skew_Q_neg,&tmp_Q_neg); //-1-1
	EFp2_ECA(&table[7],&skew_Q,&tmp_Q_neg);     //1-1
	EFp2_ECA(&table[8],&skew_Q_neg,&tmp_Q);     //-11
	
	loop_length=0;
	for(i=0; i<2; i++){
		length_s[i]=(int)mpz_sizeinbase(s[i],2);
		if(loop_length<length_s[i]){
			loop_length=length_s[i];
		}
	}
	int JSF_length;
	int JSF_binary[2][loop_length+1];
	char check[5];
	for(i=0; i<loop_length; i++){
		JSF_binary[0][i]=0;
		JSF_binary[1][i]=0;
	}
	int *JSF_pointer[2];
	JSF_pointer[0]=JSF_binary[0];
	JSF_pointer[1]=JSF_binary[1];
	Joint_sparse_form(JSF_pointer,s,&JSF_length);
	int binary[JSF_length+1];
	for(i=JSF_length; i>=0; i--){
		if(JSF_binary[1][i]==0 && JSF_binary[0][i]==0) 		binary[i]=0;
		else if(JSF_binary[1][i]==0 && JSF_binary[0][i]==1) 	binary[i]=1;
		else if(JSF_binary[1][i]==1 && JSF_binary[0][i]==0) 	binary[i]=2;
		else if(JSF_binary[1][i]==1 && JSF_binary[0][i]==1)	binary[i]=3;
		else if(JSF_binary[1][i]==0 && JSF_binary[0][i]==-1)	binary[i]=4;
		else if(JSF_binary[1][i]==-1 && JSF_binary[0][i]==0)	binary[i]=5;
		else if(JSF_binary[1][i]==-1 && JSF_binary[0][i]==-1)	binary[i]=6;
		else if(JSF_binary[1][i]==1 && JSF_binary[0][i]==-1)	binary[i]=7;
		else if(JSF_binary[1][i]==-1 && JSF_binary[0][i]==1)	binary[i]=8;
	}
	EFp2_set(&next_tmp_Q,&table[binary[JSF_length]]);
    
	for(i=JSF_length-1; i>=0; i--){
		EFp2_ECD(&next_tmp_Q,&next_tmp_Q);
		EFp2_ECA(&next_tmp_Q,&next_tmp_Q,&table[binary[i]]);
	}
	EFp2_to_EFp12_karatsuba(ANS,&next_tmp_Q);
	
	mpz_clear(buf);
	EFp2_clear(&next_tmp_Q);
	EFp2_clear(&tmp_Q);
	EFp2_clear(&tmp_Q_neg);
	EFp2_clear(&skew_Q);
	EFp2_clear(&skew_Q_neg);
	for(i=0; i<2; i++){
		mpz_clear(s[i]);
	}
	for(i=0; i<9; i++){
		EFp2_clear(&table[i]);
	}
	
    gettimeofday(&tv_end,NULL);
    G2SCM_2SPLIT_JSF=timedifference_msec(tv_start,tv_end);
}

void EFp12_G2_SCM_4split_karatsuba(EFp12 *ANS,EFp12 *Q,mpz_t scalar){
    gettimeofday(&tv_start,NULL);
    
    int i,length_s[4],loop_length;
	EFp2 next_twisted_Q,twisted_Q,twisted_Q_6x,twisted_Q_6xx,twisted_Q_36xxx,skew_Q,skew_Q_neg,skew_Q_puls1,minus_skew_Q_puls1;
	EFp2_init(&next_twisted_Q);
	EFp2_init(&twisted_Q);
	EFp2_init(&twisted_Q_6x);
	EFp2_init(&twisted_Q_6xx);
	EFp2_init(&twisted_Q_36xxx);
	EFp2_init(&skew_Q);
	EFp2_init(&skew_Q_neg);
	EFp2_init(&skew_Q_puls1);
	EFp2_init(&minus_skew_Q_puls1);
	
	mpz_t buf,A,B,s[4];
	mpz_init(buf);
	mpz_init(A);
	mpz_init(B);
	for(i=0; i<4; i++){
		mpz_init(s[i]);
	}
	EFp2 table[16];
	for(i=0; i<16; i++){
		EFp2_init(&table[i]);
	}
    
    mpz_sub_ui(buf,trace,1);
    mpz_tdiv_qr(B,A,scalar,buf);
    mpz_mul_ui(buf,X,6);
    mpz_tdiv_qr(s[1],s[0],A,buf);
    mpz_tdiv_qr(s[3],s[2],B,buf);
	
	EFp12_to_EFp2_karatsuba(&twisted_Q,Q);
	EFp2_skew_frobenius_map_p1_karatsuba(&skew_Q,&twisted_Q);
	EFp2_set_neg(&skew_Q_neg,&skew_Q);
	EFp2_set(&twisted_Q_6xx,&skew_Q);
	EFp2_ECA(&skew_Q_puls1,&skew_Q,&twisted_Q);
	EFp2_ECA(&minus_skew_Q_puls1,&skew_Q_neg,&twisted_Q);
	EFp2_skew_frobenius_map_p3_karatsuba(&minus_skew_Q_puls1,&minus_skew_Q_puls1);
	EFp2_ECA(&twisted_Q_6x,&skew_Q_puls1,&minus_skew_Q_puls1);
	EFp2_set_neg(&twisted_Q_6x,&twisted_Q_6x);
	EFp2_skew_frobenius_map_p1_karatsuba(&twisted_Q_36xxx,&twisted_Q_6x);
	table[0].infinity=1;                                    //0000
	EFp2_set(&table[1],&twisted_Q);                         //0001
	EFp2_set(&table[2],&twisted_Q_6x);                      //0010
	EFp2_ECA(&table[3],&twisted_Q_6x,&twisted_Q);           //0011
	EFp2_set(&table[4],&twisted_Q_6xx);                     //0100
	EFp2_ECA(&table[5],&twisted_Q_6xx,&twisted_Q);          //0101
	EFp2_ECA(&table[6],&twisted_Q_6xx,&twisted_Q_6x);       //0110
	EFp2_ECA(&table[7],&table[6],&twisted_Q);               //0111
	EFp2_set(&table[8],&twisted_Q_36xxx);                   //1000
	EFp2_ECA(&table[9],&twisted_Q_36xxx,&twisted_Q);        //1001
	EFp2_ECA(&table[10],&twisted_Q_36xxx,&twisted_Q_6x);    //1010
	EFp2_ECA(&table[11],&twisted_Q_36xxx,&table[3]);        //1011
	EFp2_ECA(&table[12],&twisted_Q_36xxx,&twisted_Q_6xx);   //1100
	EFp2_ECA(&table[13],&table[12],&twisted_Q);	            //1101
	EFp2_ECA(&table[14],&table[12],&twisted_Q_6x);          //1110
	EFp2_ECA(&table[15],&table[14],&twisted_Q);             //1111
	
	loop_length=0;
	for(i=0; i<4; i++){
		length_s[i]=(int)mpz_sizeinbase(s[i],2);
		if(loop_length<length_s[i]){
			loop_length=length_s[i];
		}
	}
	char binary_s[4][loop_length+1];
	char str[5],*e;
	int binary[loop_length+1];
	for(i=0; i<4; i++){
		if(length_s[i]==loop_length){
			mpz_get_str(binary_s[i],2,s[i]);
		}else{
			char binary_buf[loop_length+1];
			mpz_get_str(binary_buf,2,s[i]);
			memset(binary_s[i],'0',sizeof(binary_s[i]));
			memmove(binary_s[i]+loop_length-length_s[i],binary_buf,sizeof(binary_buf));
		}
	}
	for(i=0; i<loop_length; i++){
		sprintf(str,"%c%c%c%c",binary_s[3][i],binary_s[2][i],binary_s[1][i],binary_s[0][i]);
		binary[i]=(int)strtol(str,&e,2);
	}
	EFp2_set(&next_twisted_Q,&table[binary[0]]);
	
	for(i=1; i<loop_length; i++){
		EFp2_ECD(&next_twisted_Q,&next_twisted_Q);
		EFp2_ECA(&next_twisted_Q,&next_twisted_Q,&table[binary[i]]);
	}
	EFp2_to_EFp12_karatsuba(ANS,&next_twisted_Q);
	
	mpz_clear(buf);
	mpz_clear(A);
	mpz_clear(B);
	EFp2_clear(&next_twisted_Q);
	EFp2_clear(&twisted_Q);
	EFp2_clear(&twisted_Q_6x);
	EFp2_clear(&twisted_Q_6xx);
	EFp2_clear(&twisted_Q_36xxx);
	EFp2_clear(&skew_Q);
	EFp2_clear(&skew_Q_neg);
	EFp2_clear(&skew_Q_puls1);
	EFp2_clear(&minus_skew_Q_puls1);
	for(i=0; i<4; i++){
		mpz_clear(s[i]);
	}
	for(i=0; i<16; i++){
		EFp2_clear(&table[i]);
	}
	
    gettimeofday(&tv_end,NULL);
    G2SCM_4SPLIT=timedifference_msec(tv_start,tv_end);
}

/*----------------------------------------------------------------------------*/
//G3 EXP
void Fp12_G3_EXP_plain_CVMA(Fp12 *ANS,Fp12 *A,mpz_t scalar){
    gettimeofday(&tv_start,NULL);
    
    int i,length;
	length=(int)mpz_sizeinbase(scalar,2);
	char binary[length];
	mpz_get_str(binary,2,scalar);
	Fp12 tmp,buf;
	Fp12_init(&tmp);
	Fp12_init(&buf);
	Fp12_set(&tmp,A);
	Fp12_set(&buf,&tmp);
	
	for(i=1; binary[i]!='\0'; i++){
		Fp12_sqr_cyclotomic_CVMA(&buf,&buf);
		if(binary[i]=='1'){
			Fp12_mul_CVMA(&buf,&tmp,&buf);
		}
	}
	
	Fp12_set(ANS,&buf);
	Fp12_clear(&tmp);
	Fp12_clear(&buf);
	
    gettimeofday(&tv_end,NULL);
    G3EXP_PLAIN=timedifference_msec(tv_start,tv_end);
}

void Fp12_G3_EXP_2split_CVMA(Fp12 *ANS,Fp12 *A,mpz_t scalar){
    gettimeofday(&tv_start,NULL);
    
    int i,length_s[2],loop_length;
	Fp12 Buf,next_f,f,frobenius_f;
	Fp12_init(&Buf);
	Fp12_init(&next_f);
	Fp12_init(&f);
	Fp12_init(&frobenius_f);
	mpz_t s[2],buf;
	mpz_init(buf);
	for(i=0; i<2; i++){
		mpz_init(s[i]);
	}
	Fp12 table[4];
	for(i=0; i<4; i++){
		Fp12_init(&table[i]);
	}
	
    mpz_sub_ui(buf,trace,1);
    mpz_tdiv_qr(s[1],s[0],scalar,buf);
    
	Fp12_set(&f,A);
    Fp12_frobenius_map_p1_CVMA(&frobenius_f,&f);
    Fp_set_ui(&table[0].x0.x0.x0,1);                    //00
	Fp_set_neg(&table[0].x0.x0.x0,&table[0].x0.x0.x0);
	Fp_set(&table[0].x0.x1.x0,&table[0].x0.x0.x0);
	Fp_set(&table[0].x0.x2.x0,&table[0].x0.x0.x0);
	Fp12_set(&table[1],&f);                             //01
	Fp12_set(&table[2],&frobenius_f);                   //10
	Fp12_mul_CVMA(&table[3],&f,&frobenius_f);           //11
	
	loop_length=0;
	for(i=0; i<2; i++){
		length_s[i]=(int)mpz_sizeinbase(s[i],2);
		if(loop_length<length_s[i]){
			loop_length=length_s[i];
		}
	}
	char binary_s[2][loop_length+1];
	char str[5],*e;
	int binary[loop_length+1];
	for(i=0; i<2; i++){
		if(length_s[i]==loop_length){
			mpz_get_str(binary_s[i],2,s[i]);
		}else{
			char binary_buf[loop_length+1];
			mpz_get_str(binary_buf,2,s[i]);
			memset(binary_s[i],'0',sizeof(binary_s[i]));
			memmove(binary_s[i]+loop_length-length_s[i],binary_buf,sizeof(binary_buf));
		}
	}
	for(i=0; i<loop_length; i++){
		sprintf(str,"%c%c",binary_s[1][i],binary_s[0][i]);
		binary[i]=(int)strtol(str,&e,2);
	}
	Fp12_set(&next_f,&table[binary[0]]);
	
	for(i=1; i<loop_length; i++){
		Fp12_sqr_cyclotomic_CVMA(&next_f,&next_f);
		Fp12_mul_CVMA(&next_f,&next_f,&table[binary[i]]);
	}
	Fp12_set(ANS,&next_f);
	
	mpz_clear(buf);
	Fp12_clear(&Buf);
	Fp12_clear(&next_f);
	Fp12_clear(&f);
	Fp12_clear(&frobenius_f);
	for(i=0; i<2; i++){
		mpz_clear(s[i]);
	}
	for(i=0; i<4; i++){
		Fp12_clear(&table[i]);
	}
	
    gettimeofday(&tv_end,NULL);
    G3EXP_2SPLIT=timedifference_msec(tv_start,tv_end);
}

void Fp12_G3_EXP_2split_JSF_CVMA(Fp12 *ANS,Fp12 *A,mpz_t scalar){
    gettimeofday(&tv_start,NULL);
    
    int i,length_s[2],loop_length;
	Fp12 next_f,f,f_inv,frobenius_f,frobenius_f_inv;
	Fp12_init(&next_f);
	Fp12_init(&f);
	Fp12_init(&f_inv);
	Fp12_init(&frobenius_f);
	Fp12_init(&frobenius_f_inv);
	mpz_t s[2],buf;
	mpz_init(buf);
	for(i=0; i<2; i++){
		mpz_init(s[i]);
	}
	Fp12 table[9];
	for(i=0; i<9; i++){
		Fp12_init(&table[i]);
	}
    
    mpz_sub_ui(buf,trace,1);
    mpz_tdiv_qr(s[1],s[0],scalar,buf);
	
	Fp12_set(&f,A);
	Fp12_frobenius_map_p6_CVMA(&f_inv,&f);
	Fp12_frobenius_map_p1_CVMA(&frobenius_f,&f);
	Fp12_frobenius_map_p6_CVMA(&frobenius_f_inv,&frobenius_f);
	Fp_set_ui(&table[0].x0.x0.x0,1);                    //00
	Fp_set_neg(&table[0].x0.x0.x0,&table[0].x0.x0.x0);
	Fp_set(&table[0].x0.x1.x0,&table[0].x0.x0.x0);
	Fp_set(&table[0].x0.x2.x0,&table[0].x0.x0.x0);
	Fp12_set(&table[1],&f);                             //01
	Fp12_set(&table[2],&frobenius_f);                   //10
	Fp12_mul_CVMA(&table[3],&frobenius_f,&f);           //11
	Fp12_set(&table[4],&f_inv);                         //0-1
	Fp12_set(&table[5],&frobenius_f_inv);               //-10
	Fp12_mul_CVMA(&table[6],&frobenius_f_inv,&f_inv);   //-1-1
	Fp12_mul_CVMA(&table[7],&frobenius_f,&f_inv);       //1-1
	Fp12_mul_CVMA(&table[8],&frobenius_f_inv,&f);       //-11
    
	loop_length=0;
	for(i=0; i<2; i++){
		length_s[i]=(int)mpz_sizeinbase(s[i],2);
		if(loop_length<length_s[i]){
			loop_length=length_s[i];
		}
	}
	int JSF_length;
	int JSF_binary[2][loop_length+1];
	char check[5];
	for(i=0; i<loop_length; i++){
		JSF_binary[0][i]=0;
		JSF_binary[1][i]=0;
	}
	int *JSF_pointer[2];
	JSF_pointer[0]=JSF_binary[0];
	JSF_pointer[1]=JSF_binary[1];
	Joint_sparse_form(JSF_pointer,s,&JSF_length);
	int binary[JSF_length+1];
	for(i=JSF_length; i>=0; i--){
		if(JSF_binary[1][i]==0 && JSF_binary[0][i]==0) 		binary[i]=0;
		else if(JSF_binary[1][i]==0 && JSF_binary[0][i]==1) 	binary[i]=1;
		else if(JSF_binary[1][i]==1 && JSF_binary[0][i]==0) 	binary[i]=2;
		else if(JSF_binary[1][i]==1 && JSF_binary[0][i]==1)	binary[i]=3;
		else if(JSF_binary[1][i]==0 && JSF_binary[0][i]==-1)	binary[i]=4;
		else if(JSF_binary[1][i]==-1 && JSF_binary[0][i]==0)	binary[i]=5;
		else if(JSF_binary[1][i]==-1 && JSF_binary[0][i]==-1)	binary[i]=6;
		else if(JSF_binary[1][i]==1 && JSF_binary[0][i]==-1)	binary[i]=7;
		else if(JSF_binary[1][i]==-1 && JSF_binary[0][i]==1)	binary[i]=8;
	}
	Fp12_set(&next_f,&table[binary[JSF_length]]);
    
	for(i=JSF_length-1; i>=0; i--){
		Fp12_sqr_cyclotomic_CVMA(&next_f,&next_f);
		Fp12_mul_CVMA(&next_f,&next_f,&table[binary[i]]);
	}
	Fp12_set(ANS,&next_f);
	
	mpz_clear(buf);
	Fp12_clear(&next_f);
	Fp12_clear(&f);
	Fp12_clear(&f_inv);
	Fp12_clear(&frobenius_f);
	Fp12_clear(&frobenius_f_inv);
	for(i=0; i<2; i++){
		mpz_clear(s[i]);
	}
	for(i=0; i<9; i++){
		Fp12_clear(&table[i]);
	}
	
    gettimeofday(&tv_end,NULL);
    G3EXP_2SPLIT_JSF=timedifference_msec(tv_start,tv_end);
}

void Fp12_G3_EXP_4split_CVMA(Fp12 *ANS,Fp12 *A,mpz_t scalar){
    gettimeofday(&tv_start,NULL);
    
    int i,length_s[4],loop_length;
	Fp12 Buf;
	Fp12_init(&Buf);
	Fp12 next_f,f,f_6x,f_6xx,f_36xxx,frobenius_f,frobenius_f_inv,frobenius_f_f,frobenius_f_inv_f;
	Fp12_init(&next_f);
	Fp12_init(&f);
	Fp12_init(&f_6x);
	Fp12_init(&f_6xx);
	Fp12_init(&f_36xxx);
	Fp12_init(&frobenius_f);
	Fp12_init(&frobenius_f_inv);
	Fp12_init(&frobenius_f_f);
	Fp12_init(&frobenius_f_inv_f);
	mpz_t buf,a,b,s[4];
	mpz_init(buf);
	mpz_init(a);
	mpz_init(b);
	for(i=0; i<4; i++){
		mpz_init(s[i]);
	}
	Fp12 table[16];
	for(i=0; i<16; i++){
		Fp12_init(&table[i]);
	}
    
    mpz_sub_ui(buf,trace,1);
    mpz_tdiv_qr(b,a,scalar,buf);
    mpz_mul_ui(buf,X,6);
    mpz_tdiv_qr(s[1],s[0],a,buf);
    mpz_tdiv_qr(s[3],s[2],b,buf);
	
	Fp12_set(&f,A);
	Fp12_frobenius_map_p1_CVMA(&frobenius_f,&f);
	Fp12_set(&f_6xx,&frobenius_f);
	Fp12_frobenius_map_p6_CVMA(&frobenius_f_inv,&frobenius_f);
	Fp12_mul_CVMA(&frobenius_f_f,&frobenius_f,&f);
	Fp12_mul_CVMA(&frobenius_f_inv_f,&frobenius_f_inv,&f);
	Fp12_frobenius_map_p3_CVMA(&frobenius_f_inv_f,&frobenius_f_inv_f);
	Fp12_mul_CVMA(&f_6x,&frobenius_f_f,&frobenius_f_inv_f);
	Fp12_frobenius_map_p6_CVMA(&f_6x,&f_6x);
	Fp12_frobenius_map_p1_CVMA(&f_36xxx,&f_6x);
	Fp_set_ui(&table[0].x0.x0.x0,1);                    //0000
	Fp_set_neg(&table[0].x0.x0.x0,&table[0].x0.x0.x0);
	Fp_set(&table[0].x0.x1.x0,&table[0].x0.x0.x0);
	Fp_set(&table[0].x0.x2.x0,&table[0].x0.x0.x0);
	Fp12_set(&table[1],&f);                             //0001
	Fp12_set(&table[2],&f_6x);                          //0010
	Fp12_mul_CVMA(&table[3],&f_6x,&f);                  //0011
	Fp12_set(&table[4],&f_6xx);                         //0100
	Fp12_mul_CVMA(&table[5],&f_6xx,&f);                 //0101
	Fp12_mul_CVMA(&table[6],&f_6xx,&f_6x);              //0110
	Fp12_mul_CVMA(&table[7],&table[6],&f);              //0111
	Fp12_set(&table[8],&f_36xxx);                       //1000
	Fp12_mul_CVMA(&table[9],&f_36xxx,&f);               //1001
	Fp12_mul_CVMA(&table[10],&f_36xxx,&f_6x);           //1010
	Fp12_mul_CVMA(&table[11],&f_36xxx,&table[3]);       //1011
	Fp12_mul_CVMA(&table[12],&f_36xxx,&f_6xx);          //1100
	Fp12_mul_CVMA(&table[13],&table[12],&f);            //1101
	Fp12_mul_CVMA(&table[14],&table[12],&f_6x);         //1110
	Fp12_mul_CVMA(&table[15],&table[14],&f);            //1111
	
	loop_length=0;
	for(i=0; i<4; i++){
		length_s[i]=(int)mpz_sizeinbase(s[i],2);
		if(loop_length<length_s[i]){
			loop_length=length_s[i];
		}
	}
	char binary_s[4][loop_length+1];
	char str[5],*e;
	int binary[loop_length+1];
	for(i=0; i<4; i++){
		if(length_s[i]==loop_length){
			mpz_get_str(binary_s[i],2,s[i]);
		}else{
			char binary_buf[loop_length+1];
			mpz_get_str(binary_buf,2,s[i]);
			memset(binary_s[i],'0',sizeof(binary_s[i]));
			memmove(binary_s[i]+loop_length-length_s[i],binary_buf,sizeof(binary_buf));
		}
	}
	for(i=0; i<loop_length; i++){
		sprintf(str,"%c%c%c%c",binary_s[3][i],binary_s[2][i],binary_s[1][i],binary_s[0][i]);
		binary[i]=(int)strtol(str,&e,2);
	}
	Fp12_set(&next_f,&table[binary[0]]);
	
	for(i=1; i<loop_length; i++){
		Fp12_sqr_cyclotomic_CVMA(&next_f,&next_f);
		Fp12_mul_CVMA(&next_f,&next_f,&table[binary[i]]);
	}
	Fp12_set(ANS,&next_f);
	
	mpz_clear(buf);
	mpz_clear(a);
	mpz_clear(b);
	Fp12_clear(&Buf);
	Fp12_clear(&next_f);
	Fp12_clear(&f);
	Fp12_clear(&f_6x);
	Fp12_clear(&f_6xx);
	Fp12_clear(&f_36xxx);
	Fp12_clear(&frobenius_f);
	Fp12_clear(&frobenius_f_inv);
	Fp12_clear(&frobenius_f_f);
	Fp12_clear(&frobenius_f_inv_f);
	for(i=0; i<4; i++){
		mpz_clear(s[i]);
	}
	for(i=0; i<16; i++){
		Fp12_clear(&table[i]);
	}
	
    gettimeofday(&tv_end,NULL);
    G3EXP_4SPLIT=timedifference_msec(tv_start,tv_end);
}

/*----------------------------------------------------------------------------*/
//settings
void BN12_init(){
    init_parameters();
    generate_X();
    if(generate_prime()==1 && generate_order()==1){
        generate_trace();
        weil();
        get_epsilon();
        get_Two_inv();
        set_basis();
        set_frobenius_constant_karatsuba();
        set_frobenius_constant_CVMA();
        set_curve_parameter();
        get_basis_conversion_matrix();
    }else{
        BN12_clear();
        printf("error : prime\nexit\n");
        exit(1);
    }
}

void init_parameters(){
    int i,j;

    mpz_init(X);
    mpz_init(prime);
    mpz_init(order);
    mpz_init(trace);
    mpz_init(EFp_total);
    mpz_init(EFp12_total);
    mpz_init(curve_b);
    
    for(i=0; i<X_length+1; i++){
        X_binary[i]=0;
    }
    for(i=0; i<X6_2_length+1; i++){
        X6_2_binary[i]=0;
    }
    
    mpz_init(epsilon1);
    mpz_init(epsilon2);
    mpz_init(Two_inv);
    Fp2_init(&Alpha_1);
    Fp2_init(&Alpha_1_inv);
    
    for(i=0; i<12; i++){
        Fp2_init(&frobenius_constant_CVMA[i]);
        for(j=0; j<6; j++){
            Fp2_init(&frobenius_constant_karatsuba[i][j]);
        }
        for(j=0; j<2; j++){
            Fp2_init(&skew_frobenius_constant_karatsuba[i][j]);
        }
    }
    
    Matrix_init(&Matrix_karatsuba_to_CVMA);
    Matrix_init(&Matrix_CVMA_to_karatsuba);
    
    Fp_init(&TMP1_FP);
    Fp_init(&TMP2_FP);
    Fp_init(&TMP3_FP);
    Fp_init(&TMP4_FP);
    Fp_init(&TMP5_FP);
    Fp2_init(&TMP1_FP2);
    Fp2_init(&TMP2_FP2);
    Fp2_init(&TMP3_FP2);
    Fp2_init(&TMP4_FP2);
    Fp2_init(&TMP5_FP2);
    Fp2_init(&TMP6_FP2);
    Fp2_init(&TMP7_FP2);
    Fp2_init(&TMP8_FP2);
    Fp2_init(&TMP9_FP2);
    Fp2_init(&TMP10_FP2);
    Fp2_init(&TMP11_FP2);
    Fp2_init(&TMP12_FP2);
    Fp2_init(&TMP13_FP2);
    Fp2_init(&TMP14_FP2);
    Fp6_init(&TMP1_FP6);
    Fp6_init(&TMP2_FP6);
    Fp6_init(&TMP3_FP6);
    Fp6_init(&TMP4_FP6);
    Fp12_init(&TMP1_FP12);
    Fp12_init(&TMP2_FP12);
    Fp12_init(&TMP3_FP12);
    EFp_init(&TMP1_EFP);
    EFp_init(&TMP2_EFP);
    EFp2_init(&TMP1_EFP2);
    EFp2_init(&TMP2_EFP2);
    EFp6_init(&TMP1_EFP6);
    EFp6_init(&TMP2_EFP6);
    EFp12_init(&TMP1_EFP12);
    EFp12_init(&TMP2_EFP12);
}

void generate_X(){
    int i;
    mpz_t buf;
    mpz_init(buf);
    
    X_binary[114]=1;
    X_binary[20]=1;
    X_binary[18]=1;
    X_binary[12]=-1;
    X_binary[6]=1;
    X_binary[0]=-1;
    
    X6_2_binary[116]=1;
    X6_2_binary[115]=1;
    X6_2_binary[22]=1;
    X6_2_binary[21]=1;
    X6_2_binary[20]=1;
    X6_2_binary[19]=1;
    X6_2_binary[14]=-1;
    X6_2_binary[13]=-1;
    X6_2_binary[8]=1;
    X6_2_binary[7]=1;
    X6_2_binary[2]=-1;
    
    mpz_set_ui(X,0);
    for(i=X_length; i>=0; i--){
        if(X_binary[i]==1){
            mpz_ui_pow_ui(buf,2,i);
            mpz_add(X,X,buf);
        }else if(X_binary[i]==-1){
            mpz_ui_pow_ui(buf,2,i);
            mpz_sub(X,X,buf);
        }
    }
    
    mpz_clear(buf);
}

int generate_prime(){
    mpz_t buf,result;
    mpz_init(buf);
    mpz_init(result);
    
    mpz_pow_ui(buf,X,4);
    mpz_mul_ui(buf,buf,36);
    mpz_set(result,buf);
    mpz_pow_ui(buf,X,3);
    mpz_mul_ui(buf,buf,36);
    mpz_add(result,result,buf);
    mpz_pow_ui(buf,X,2);
    mpz_mul_ui(buf,buf,24);
    mpz_add(result,result,buf);
    mpz_mul_ui(buf,X,6);
    mpz_add(result,result,buf);
    mpz_add_ui(result,result,1);
    
    if(mpz_probab_prime_p(result,25)==0){
        mpz_clear(buf);
        mpz_clear(result);
        return 0;
    }else{
        mpz_set(prime,result);
        mpz_clear(buf);
        mpz_clear(result);
        return 1;
    }
}

int generate_order(){
    mpz_t buf,result;
    mpz_init(buf);
    mpz_init(result);
    
    mpz_pow_ui(buf,X,4);
    mpz_mul_ui(buf,buf,36);
    mpz_set(result,buf);
    mpz_pow_ui(buf,X,3);
    mpz_mul_ui(buf,buf,36);
    mpz_add(result,result,buf);
    mpz_pow_ui(buf,X,2);
    mpz_mul_ui(buf,buf,18);
    mpz_add(result,result,buf);
    mpz_mul_ui(buf,X,6);
    mpz_add(result,result,buf);
    mpz_add_ui(result,result,1);
    
    if(mpz_probab_prime_p(result,25)==0){
        mpz_clear(buf);
        mpz_clear(result);
        return 0;
    }else{
        mpz_set(order,result);
        mpz_clear(buf);
        mpz_clear(result);
        return 1;
    }
}

void generate_trace(){
    mpz_t buf;
    mpz_init(buf);
    
    mpz_pow_ui(buf,X,2);
    mpz_mul_ui(buf,buf,6);
    mpz_add_ui(trace,buf,1);
    
    mpz_clear(buf);
}

void weil(){
    mpz_t t2,t6,t12,p2,p6,buf;
    mpz_init(t2);
    mpz_init(t6);
    mpz_init(t12);
    mpz_init(p2);
    mpz_init(p6);
    mpz_init(buf);
    
    mpz_add_ui(buf,prime,1);
    mpz_sub(EFp_total,buf,trace);
    
    mpz_pow_ui(t2,trace,2);
    mpz_mul_ui(buf,prime,2);
    mpz_sub(t2,t2,buf);
    mpz_pow_ui(p2,prime,2);
    mpz_pow_ui(t6,t2,3);
    mpz_mul(buf,t2,p2);
    mpz_mul_ui(buf,buf,3);
    mpz_sub(t6,t6,buf);
    mpz_pow_ui(p6,p2,3);
    mpz_pow_ui(t12,t6,2);
    mpz_mul_ui(buf,p6,2);
    mpz_sub(t12,t12,buf);
    
    mpz_pow_ui(buf,p6,2);
    mpz_sub(buf,buf,t12);
    mpz_add_ui(EFp12_total,buf,1);
    
    mpz_clear(t2);
    mpz_clear(t6);
    mpz_clear(t12);
    mpz_clear(p2);
    mpz_clear(p6);
    mpz_clear(buf);
}

void get_epsilon(){
    Fp inv,buf,result1,result2;
    Fp_init(&inv);
    Fp_init(&buf);
    Fp_init(&result1);
    Fp_init(&result2);
    
    Fp_set_ui(&buf,2);
    Fp_inv(&inv,&buf);
    mpz_sub_ui(buf.x0,prime,3);
    Fp_sqrt(&buf,&buf);
    Fp_sub_ui(&buf,&buf,1);
    Fp_mul(&result1,&buf,&inv);
    Fp_mul(&result2,&result1,&result1);
    mpz_set(epsilon1,result1.x0);
    mpz_set(epsilon2,result2.x0);
    
    Fp_clear(&inv);
    Fp_clear(&buf);
    Fp_clear(&result1);
    Fp_clear(&result2);
}

void get_Two_inv(){
    mpz_set_ui(Two_inv,2);
    mpz_invert(Two_inv,Two_inv,prime);
}

void set_basis(){
    Fp2_set_ui(&Alpha_1,1);
    Fp2_inv(&Alpha_1_inv,&Alpha_1);
}

void set_frobenius_constant_karatsuba(){
    mpz_t exp,p1,p2,p3,p4,p6,p8,p10;
    mpz_init(exp);
    mpz_init(p1);
    mpz_init(p2);
    mpz_init(p3);
    mpz_init(p4);
    mpz_init(p6);
    mpz_init(p8);
    mpz_init(p10);
    Fp2 tmp_Fp2;
    Fp2_init(&tmp_Fp2);
    Fp two_constant1,two_constant2;
    Fp_init(&two_constant1);
    Fp_init(&two_constant2);
    
    mpz_set(p1,prime);
    mpz_mul(p2,p1,p1);
    mpz_mul(p3,p2,p1);
    mpz_mul(p4,p3,p1);
    mpz_mul(p6,p4,p2);
    mpz_mul(p8,p6,p2);
    mpz_mul(p10,p8,p2);
    
    mpz_sub_ui(exp,p1,1);
    mpz_tdiv_q_ui(exp,exp,3);
    Fp_set_ui(&two_constant1,2);
    Fp_pow(&two_constant1,&two_constant1,exp);
    Fp_mul(&two_constant2,&two_constant1,&two_constant1);
    
    //p1
    mpz_sub_ui(exp,p1,1);
    mpz_tdiv_q_ui(exp,exp,2);
    Fp2_pow(&tmp_Fp2,&Alpha_1,exp);
    Fp_set_ui(&frobenius_constant_karatsuba[f_p1][0].x0,1);
    Fp_set(&frobenius_constant_karatsuba[f_p1][1].x0,&two_constant1);
    Fp_set(&frobenius_constant_karatsuba[f_p1][2].x0,&two_constant2);
    Fp2_set(&frobenius_constant_karatsuba[f_p1][3],&tmp_Fp2);
    Fp2_mul_mpz(&frobenius_constant_karatsuba[f_p1][4],&tmp_Fp2,two_constant1.x0);
    Fp2_mul_mpz(&frobenius_constant_karatsuba[f_p1][5],&tmp_Fp2,two_constant2.x0);
    //p2
    mpz_sub_ui(exp,p2,1);
    mpz_tdiv_q_ui(exp,exp,2);
    Fp2_pow(&tmp_Fp2,&Alpha_1,exp);
    Fp_set_ui(&frobenius_constant_karatsuba[f_p2][0].x0,1);
    Fp_set(&frobenius_constant_karatsuba[f_p2][1].x0,&two_constant2);
    Fp_set(&frobenius_constant_karatsuba[f_p2][2].x0,&two_constant1);
    Fp2_set(&frobenius_constant_karatsuba[f_p2][3],&tmp_Fp2);
    Fp2_mul_mpz(&frobenius_constant_karatsuba[f_p2][4],&tmp_Fp2,two_constant2.x0);
    Fp2_mul_mpz(&frobenius_constant_karatsuba[f_p2][5],&tmp_Fp2,two_constant1.x0);
    //p3
    mpz_sub_ui(exp,p3,1);
    mpz_tdiv_q_ui(exp,exp,2);
    Fp2_pow(&tmp_Fp2,&Alpha_1,exp);
    Fp_set_ui(&frobenius_constant_karatsuba[f_p3][0].x0,1);
    Fp_set_ui(&frobenius_constant_karatsuba[f_p3][1].x0,1);
    Fp_set_ui(&frobenius_constant_karatsuba[f_p3][2].x0,1);
    Fp2_set(&frobenius_constant_karatsuba[f_p3][3],&tmp_Fp2);
    Fp2_set(&frobenius_constant_karatsuba[f_p3][4],&tmp_Fp2);
    Fp2_set(&frobenius_constant_karatsuba[f_p3][5],&tmp_Fp2);
    //p4
    mpz_sub_ui(exp,p4,1);
    mpz_tdiv_q_ui(exp,exp,2);
    Fp2_pow(&tmp_Fp2,&Alpha_1,exp);
    Fp_set_ui(&frobenius_constant_karatsuba[f_p4][0].x0,1);
    Fp_set(&frobenius_constant_karatsuba[f_p4][1].x0,&two_constant1);
    Fp_set(&frobenius_constant_karatsuba[f_p4][2].x0,&two_constant2);
    Fp2_set(&frobenius_constant_karatsuba[f_p4][3],&tmp_Fp2);
    Fp2_mul_mpz(&frobenius_constant_karatsuba[f_p4][4],&tmp_Fp2,two_constant1.x0);
    Fp2_mul_mpz(&frobenius_constant_karatsuba[f_p4][5],&tmp_Fp2,two_constant2.x0);
    //p6
    mpz_sub_ui(exp,p6,1);
    mpz_tdiv_q_ui(exp,exp,2);
    Fp2_pow(&tmp_Fp2,&Alpha_1,exp);
    Fp_set_ui(&frobenius_constant_karatsuba[f_p6][0].x0,1);
    Fp_set_ui(&frobenius_constant_karatsuba[f_p6][1].x0,1);
    Fp_set_ui(&frobenius_constant_karatsuba[f_p6][2].x0,1);
    Fp2_set(&frobenius_constant_karatsuba[f_p6][3],&tmp_Fp2);
    Fp2_set(&frobenius_constant_karatsuba[f_p6][4],&tmp_Fp2);
    Fp2_set(&frobenius_constant_karatsuba[f_p6][5],&tmp_Fp2);
    //p8
    mpz_sub_ui(exp,p8,1);
    mpz_tdiv_q_ui(exp,exp,2);
    Fp2_pow(&tmp_Fp2,&Alpha_1,exp);
    Fp_set_ui(&frobenius_constant_karatsuba[f_p8][0].x0,1);
    Fp_set(&frobenius_constant_karatsuba[f_p8][1].x0,&two_constant2);
    Fp_set(&frobenius_constant_karatsuba[f_p8][2].x0,&two_constant1);
    Fp2_set(&frobenius_constant_karatsuba[f_p8][3],&tmp_Fp2);
    Fp2_mul_mpz(&frobenius_constant_karatsuba[f_p8][4],&tmp_Fp2,two_constant2.x0);
    Fp2_mul_mpz(&frobenius_constant_karatsuba[f_p8][5],&tmp_Fp2,two_constant1.x0);
    //p10
    mpz_sub_ui(exp,p10,1);
    mpz_tdiv_q_ui(exp,exp,2);
    Fp2_pow(&tmp_Fp2,&Alpha_1,exp);
    Fp_set_ui(&frobenius_constant_karatsuba[f_p10][0].x0,1);
    Fp_set(&frobenius_constant_karatsuba[f_p10][1].x0,&two_constant1);
    Fp_set(&frobenius_constant_karatsuba[f_p10][2].x0,&two_constant2);
    Fp2_set(&frobenius_constant_karatsuba[f_p10][3],&tmp_Fp2);
    Fp2_mul_mpz(&frobenius_constant_karatsuba[f_p10][4],&tmp_Fp2,two_constant1.x0);
    Fp2_mul_mpz(&frobenius_constant_karatsuba[f_p10][5],&tmp_Fp2,two_constant2.x0);
    
    //skew_frobenius_1
	mpz_sub_ui(exp,prime,1);
	Fp2_pow(&tmp_Fp2,&Alpha_1,exp);
	Fp2_mul(&skew_frobenius_constant_karatsuba[f_p1][0],&tmp_Fp2,&frobenius_constant_karatsuba[f_p1][2]);
	Fp2_mul(&skew_frobenius_constant_karatsuba[f_p1][1],&tmp_Fp2,&frobenius_constant_karatsuba[f_p1][3]);
	
	//skew_frobenius_2
	mpz_sub_ui(exp,p2,1);
	Fp2_pow(&tmp_Fp2,&Alpha_1,exp);
	Fp2_mul(&skew_frobenius_constant_karatsuba[f_p2][0],&tmp_Fp2,&frobenius_constant_karatsuba[f_p2][2]);
	Fp2_mul(&skew_frobenius_constant_karatsuba[f_p2][1],&tmp_Fp2,&frobenius_constant_karatsuba[f_p2][3]);
	
	//skew_frobenius_3
	mpz_sub_ui(exp,p3,1);
	Fp2_pow(&tmp_Fp2,&Alpha_1,exp);
	Fp2_mul(&skew_frobenius_constant_karatsuba[f_p3][0],&tmp_Fp2,&frobenius_constant_karatsuba[f_p3][2]);
	Fp2_mul(&skew_frobenius_constant_karatsuba[f_p3][1],&tmp_Fp2,&frobenius_constant_karatsuba[f_p3][3]);
	
	//skew_frobenius_10
	mpz_sub_ui(exp,p10,1);
	Fp2_pow(&tmp_Fp2,&Alpha_1,exp);
	Fp2_mul(&skew_frobenius_constant_karatsuba[f_p10][0],&tmp_Fp2,&frobenius_constant_karatsuba[f_p10][2]);
	Fp2_mul(&skew_frobenius_constant_karatsuba[f_p10][1],&tmp_Fp2,&frobenius_constant_karatsuba[f_p10][3]);
    
    Fp2_clear(&tmp_Fp2);
    Fp_clear(&two_constant1);
    Fp_clear(&two_constant2);
    mpz_clear(exp);
    mpz_clear(p1);
    mpz_clear(p2);
    mpz_clear(p3);
    mpz_clear(p4);
    mpz_clear(p6);
    mpz_clear(p8);
    mpz_clear(p10);
}

void set_frobenius_constant_CVMA(){
    mpz_t exp,p1,p2,p3,p4,p6,p8,p10;
    mpz_init(exp);
    mpz_init(p1);
    mpz_init(p2);
    mpz_init(p3);
    mpz_init(p4);
    mpz_init(p6);
    mpz_init(p8);
    mpz_init(p10);
    Fp2 tmp_Fp2;
    Fp2_init(&tmp_Fp2);
    Fp6 tmp_Fp6;
    Fp6_init(&tmp_Fp6);
    
    mpz_set(p1,prime);
    mpz_mul(p2,p1,p1);
    mpz_mul(p3,p2,p1);
    mpz_mul(p4,p3,p1);
    mpz_mul(p6,p4,p2);
    mpz_mul(p8,p6,p2);
    mpz_mul(p10,p8,p2);
    
    //p1
    mpz_sub_ui(exp,p1,1);
    mpz_tdiv_q_ui(exp,exp,2);
    Fp2_pow(&frobenius_constant_CVMA[f_p1],&Alpha_1,exp);
    //p2
    mpz_sub_ui(exp,p2,1);
    mpz_tdiv_q_ui(exp,exp,2);
    Fp2_pow(&frobenius_constant_CVMA[f_p2],&Alpha_1,exp);
    //p3
    mpz_sub_ui(exp,p3,1);
    mpz_tdiv_q_ui(exp,exp,2);
    Fp2_pow(&frobenius_constant_CVMA[f_p3],&Alpha_1,exp);
    //p4
    mpz_sub_ui(exp,p4,1);
    mpz_tdiv_q_ui(exp,exp,2);
    Fp2_pow(&frobenius_constant_CVMA[f_p4],&Alpha_1,exp);
    //p6
    mpz_sub_ui(exp,p6,1);
    mpz_tdiv_q_ui(exp,exp,2);
    Fp2_pow(&frobenius_constant_CVMA[f_p6],&Alpha_1,exp);
    //p8
    mpz_sub_ui(exp,p8,1);
    mpz_tdiv_q_ui(exp,exp,2);
    Fp2_pow(&frobenius_constant_CVMA[f_p8],&Alpha_1,exp);
    //p10
    mpz_sub_ui(exp,p10,1);
    mpz_tdiv_q_ui(exp,exp,2);
    Fp2_pow(&frobenius_constant_CVMA[f_p10],&Alpha_1,exp);
    
    Fp2_clear(&tmp_Fp2);
    Fp6_clear(&tmp_Fp6);
    mpz_clear(exp);
    mpz_clear(p1);
    mpz_clear(p2);
    mpz_clear(p3);
    mpz_clear(p4);
    mpz_clear(p6);
    mpz_clear(p8);
    mpz_clear(p10);
}

void set_curve_parameter(){
    mpz_set_ui(curve_b,2);
}

void BN12_print_parameters(){
    printf("====================================================================================\n");
    printf("BN12\n\n");
    gmp_printf("parameters\n");
    gmp_printf("X     (%dbit length) : %Zd \n",(int)mpz_sizeinbase(X,2),X);
    gmp_printf("prime (%dbit length) : %Zd \n",(int)mpz_sizeinbase(prime,2),prime);
    gmp_printf("order (%dbit length) : %Zd \n",(int)mpz_sizeinbase(order,2),order);
    gmp_printf("trace (%dbit length) : %Zd \n",(int)mpz_sizeinbase(trace,2),trace);
    
    gmp_printf("\nelliptic curve\n");
    gmp_printf("E:y^2=x^3+2\n");
    
    gmp_printf("\nmodulo polynomial (karatsuba)\n");
    gmp_printf("Fp2 = Fp[alpha]/(alpha^2+1)\n");
    gmp_printf("Fp6 = Fp2[omega]/(omega^3-2)\n");
    gmp_printf("Fp12= Fp6[beta]/(beta^2-(alpha+1))\n");
    
    gmp_printf("\nmodulo polynomial (CVMA)\n");
    gmp_printf("Fp2 = Fp[alpha]/(alpha^2+1)\n");
    gmp_printf("Fp6 = Fp2[omega]/(omega^6+omega^5+omega^4+omega^3+omega^2+omega+1)\n");
    gmp_printf("Fp12= Fp6[beta]/(beta^2+(alpha+1))\n");
}

void BN12_clear(){
    int i,j;
    
    mpz_clear(X);
    mpz_clear(prime);
    mpz_clear(order);
    mpz_clear(trace);
    mpz_clear(EFp_total);
    mpz_clear(EFp12_total);
    mpz_clear(curve_b);
    mpz_clear(epsilon1);
    mpz_clear(epsilon2);
    mpz_clear(Two_inv);
    Fp2_clear(&Alpha_1);
    Fp2_clear(&Alpha_1_inv);
    
    for(i=0; i<12; i++){
        Fp2_clear(&frobenius_constant_CVMA[i]);
        for(j=0; j<6; j++){
            Fp2_clear(&frobenius_constant_karatsuba[i][j]);
        }
        for(j=0; j<2; j++){
            Fp2_clear(&skew_frobenius_constant_karatsuba[i][j]);
        }
    }
    
    Matrix_clear(&Matrix_karatsuba_to_CVMA);
    Matrix_clear(&Matrix_CVMA_to_karatsuba);
    
    Fp_clear(&TMP1_FP);
    Fp_clear(&TMP2_FP);
    Fp_clear(&TMP3_FP);
    Fp_clear(&TMP4_FP);
    Fp_clear(&TMP5_FP);
    Fp2_clear(&TMP1_FP2);
    Fp2_clear(&TMP2_FP2);
    Fp2_clear(&TMP3_FP2);
    Fp2_clear(&TMP4_FP2);
    Fp2_clear(&TMP5_FP2);
    Fp2_clear(&TMP6_FP2);
    Fp2_clear(&TMP7_FP2);
    Fp2_clear(&TMP8_FP2);
    Fp2_clear(&TMP9_FP2);
    Fp2_clear(&TMP10_FP2);
    Fp2_clear(&TMP11_FP2);
    Fp2_clear(&TMP12_FP2);
    Fp2_clear(&TMP13_FP2);
    Fp2_clear(&TMP14_FP2);
    Fp6_clear(&TMP1_FP6);
    Fp6_clear(&TMP2_FP6);
    Fp6_clear(&TMP3_FP6);
    Fp6_clear(&TMP4_FP6);
    Fp12_clear(&TMP1_FP12);
    Fp12_clear(&TMP2_FP12);
    Fp12_clear(&TMP3_FP12);
    EFp_clear(&TMP1_EFP);
    EFp_clear(&TMP2_EFP);
    EFp2_clear(&TMP1_EFP2);
    EFp2_clear(&TMP2_EFP2);
    EFp6_clear(&TMP1_EFP6);
    EFp6_clear(&TMP2_EFP6);
    EFp12_clear(&TMP1_EFP12);
    EFp12_clear(&TMP2_EFP12);
}

void get_basis_conversion_matrix(){
    int i,j;
    Fp6 tmp1,tmp2,tmp3;
    Fp6_init(&tmp1);
    Fp6_init(&tmp2);
    Fp6_init(&tmp3);
    MATRIX tmp_karatsuba,tmp_karatsuba_inv,tmp_CVMA,tmp_CVMA_inv,test;
    Matrix_init(&tmp_karatsuba);
    Matrix_init(&tmp_karatsuba_inv);
    Matrix_init(&tmp_CVMA);
    Matrix_init(&tmp_CVMA_inv);
    Matrix_init(&test);
    mpz_t exp;
    mpz_init(exp);
    mpz_pow_ui(exp,prime,12);
    mpz_sub_ui(exp,exp,1);
    mpz_tdiv_q_ui(exp,exp,7);
    
    //karatsuba
    mpz_set_str(tmp1.x0.x0.x0,"3626936454119470773250005253042749723227359392902167519969643608576607899049219594514785306544179722478682698644649040622360947873115419192",10);
    mpz_set_str(tmp1.x1.x0.x0,"8004388770406465558461022206574004656018235354390130437084847809138151064715822062832032864402025736760200971431015401041382927143090864",10);
    mpz_set_str(tmp1.x2.x0.x0,"3166414805863787863474668763893866983421198258388435458393153737179600344815135451334905411516494232181883762907294420810687765201950499789",10);
    Fp6_sqr_karatsuba(&tmp2,&tmp1);
    Fp6_mul_karatsuba(&tmp3,&tmp2,&tmp1);
    
    Fp_set(&tmp_karatsuba.data[0][0],&tmp1.x0.x0);
    Fp_set(&tmp_karatsuba.data[0][1],&tmp1.x1.x0);
    Fp_set(&tmp_karatsuba.data[0][2],&tmp1.x2.x0);
    Fp_set(&tmp_karatsuba.data[1][0],&tmp2.x0.x0);
    Fp_set(&tmp_karatsuba.data[1][1],&tmp2.x1.x0);
    Fp_set(&tmp_karatsuba.data[1][2],&tmp2.x2.x0);
    Fp_set(&tmp_karatsuba.data[2][0],&tmp3.x0.x0);
    Fp_set(&tmp_karatsuba.data[2][1],&tmp3.x1.x0);
    Fp_set(&tmp_karatsuba.data[2][2],&tmp3.x2.x0);
    //CVMA
    mpz_set_str(tmp1.x0.x0.x0,"5503613203322570470031868916911126013223080721880281088192322836989779752610277231471441539726623644446745085675019092845269779616806097622",10);
    mpz_set_str(tmp1.x1.x0.x0,"6101079443215441786667218566046361919496302661527639029494168169072006857558875181701765359716169813667778753481566651625471231172643893173",10);
    mpz_set_str(tmp1.x2.x0.x0,"4308680723536827836761169618640654200676636842585565205588632172825325542713081331010793899747531306004677750061923975284866876505130506521",10);
    Fp6_sqr_CVMA(&tmp2,&tmp1);
    Fp6_mul_CVMA(&tmp3,&tmp2,&tmp1);
    Fp_set(&tmp_CVMA.data[0][0],&tmp1.x0.x0);
    Fp_set(&tmp_CVMA.data[0][1],&tmp1.x1.x0);
    Fp_set(&tmp_CVMA.data[0][2],&tmp1.x2.x0);
    Fp_set(&tmp_CVMA.data[1][0],&tmp2.x0.x0);
    Fp_set(&tmp_CVMA.data[1][1],&tmp2.x1.x0);
    Fp_set(&tmp_CVMA.data[1][2],&tmp2.x2.x0);
    Fp_set(&tmp_CVMA.data[2][0],&tmp3.x0.x0);
    Fp_set(&tmp_CVMA.data[2][1],&tmp3.x1.x0);
    Fp_set(&tmp_CVMA.data[2][2],&tmp3.x2.x0);
    
    Matrix_inv(&tmp_karatsuba_inv,&tmp_karatsuba);
    Matrix_inv(&tmp_CVMA_inv,&tmp_CVMA);
    Matrix_mul(&Matrix_karatsuba_to_CVMA,&tmp_karatsuba_inv,&tmp_CVMA);
    Matrix_mul(&Matrix_CVMA_to_karatsuba,&tmp_CVMA_inv,&tmp_karatsuba);
    Matrix_mul(&test,&Matrix_CVMA_to_karatsuba,&Matrix_karatsuba_to_CVMA);
    
    Matrix_clear(&tmp_karatsuba);
    Matrix_clear(&tmp_karatsuba_inv);
    Matrix_clear(&tmp_CVMA);
    Matrix_clear(&tmp_CVMA_inv);
    Matrix_clear(&test);
    Fp6_clear(&tmp1);
    Fp6_clear(&tmp2);
    Fp6_clear(&tmp3);
    mpz_clear(exp);
    
}

/*----------------------------------------------------------------------------*/
//Matrix arithmetic
void Matrix_init(MATRIX *M){
    int i,j;
    for(i=0; i<3; i++){
        for(j=0; j<3; j++){
            Fp_init(&M->data[i][j]);
        }
    }
}

void Matrix_set(MATRIX *ANS,MATRIX *M){
    int i,j;
    for(i=0; i<3; i++){
        for(j=0; j<3; j++){
            Fp_set(&ANS->data[i][j],&M->data[i][j]);
        }
    }
}

void Matrix_clear(MATRIX *M){
    int i,j;
    for(i=0; i<3; i++){
        for(j=0; j<3; j++){
            Fp_clear(&M->data[i][j]);
        }
    }
}

void Matrix_printf(MATRIX *M){
    int i,j;
    for(i=0; i<3; i++){
        for(j=0; j<3; j++){
            Fp_printf(&M->data[i][j],"");
            printf("\n");
        }
    }
}

int Matrix_cmp(MATRIX *M,MATRIX *N){
    int i,j;
    for(i=0; i<3; i++){
        for(j=0; j<3; j++){
            if(Fp_cmp(&M->data[i][j],&N->data[i][j])!=0) return 1;
        }
    }
    return 0;
}

void Matrix_inv(MATRIX *ANS,MATRIX *A){
    int i,j,k;
    Fp tmp0,tmp1;
    Fp_init(&tmp0);
    Fp_init(&tmp1);
    MATRIX tmp_A,ans;
    Matrix_init(&tmp_A);
    Matrix_init(&ans);
    
    Matrix_set(&tmp_A,A);
    for(i=0; i<3; i++){
        for(j=0; j<3; j++){
            if(i==j){
                Fp_set_ui(&ans.data[i][j],1);
            }
        }
    }
    
    for(i=0; i<3; i++){
        Fp_inv(&tmp0,&tmp_A.data[i][i]);
        for(j=0; j<3; j++){
            Fp_mul(&tmp_A.data[i][j],&tmp_A.data[i][j],&tmp0);
            Fp_mul(&ans.data[i][j],&ans.data[i][j],&tmp0);
        }
        for(j=0; j<3; j++){
            if(i!=j){
                Fp_set(&tmp0,&tmp_A.data[j][i]);
                for(k=0; k<3; k++){
                    Fp_mul(&tmp1,&tmp_A.data[i][k],&tmp0);
                    Fp_sub(&tmp_A.data[j][k],&tmp_A.data[j][k],&tmp1);
                    Fp_mul(&tmp1,&ans.data[i][k],&tmp0);
                    Fp_sub(&ans.data[j][k],&ans.data[j][k],&tmp1);
                }
            }
        }
    }
    Matrix_set(ANS,&ans);
    
    Matrix_clear(&tmp_A);
    Matrix_clear(&ans);
    Fp_clear(&tmp0);
    Fp_clear(&tmp1);
}

void Matrix_mul(MATRIX *ANS,MATRIX *A,MATRIX *B){
    int i,j,k;
    Fp tmp;
    Fp_init(&tmp);
    MATRIX ans;
    Matrix_init(&ans);
    
    for(i=0; i<3; i++){
        for(j=0; j<3; j++){
            for(k=0; k<3; k++){
                Fp_mul(&tmp,&A->data[i][k],&B->data[k][j]);
                Fp_add(&ans.data[i][j],&ans.data[i][j],&tmp);
            }
        }
    }
    Matrix_set(ANS,&ans);
    
    Fp_clear(&tmp);
    Matrix_clear(&ans);
}

/*============================================================================*/
/* Test                                                                       */
/*============================================================================*/
//time
float timedifference_msec(struct timeval t0, struct timeval t1){
    return (t1.tv_sec - t0.tv_sec) * 1000.0f + (t1.tv_usec - t0.tv_usec) / 1000.0f;
}

float timedifference_usec(struct timeval t0, struct timeval t1){
    return (t1.tv_sec - t0.tv_sec) + (t1.tv_usec - t0.tv_usec);
}

/*----------------------------------------------------------------------------*/
//cost
void Init_mpz_Cost(struct mpz_Cost *cost){
    cost->mpz_mul=0;
    cost->mpz_mul_ui=0;
    cost->mpz_sqr=0;
    cost->mpz_add=0;
    cost->mpz_add_ui=0;
    cost->mpz_invert=0;
}

void Print_mpz_Cost(struct mpz_Cost *cost,char *str){
    printf("%s",str);
    printf("mpz_mul,mpz_mul_ui,mpz_sqr,mpz_add,mpz_add_ui,mpz_invert\n");
    printf("%ld,",cost->mpz_mul);
    printf("%ld,",cost->mpz_mul_ui);
    printf("%ld,",cost->mpz_sqr);
    printf("%ld,",cost->mpz_add);
    printf("%ld,",cost->mpz_add_ui);
    printf("%ld",cost->mpz_invert);
    printf("\n");
}

void Init_Fp_Cost(struct Fp_Cost *cost){
    cost->Fp_mul=0;
    cost->Fp_mul_mpz=0;
    cost->Fp_mul_ui=0;
    cost->Fp_sqr=0;
    cost->Fp_basis=0;
    cost->Fp_add=0;
    cost->Fp_add_mpz=0;
    cost->Fp_add_ui=0;
    cost->Fp_inv=0;
    cost->Fp_neg=0;
}

void Print_Fp_Cost(struct Fp_Cost *cost,char *str){
    printf("%s",str);
    printf("Fp_mul,Fp_mul_mpz,Fp_mul_ui,Fp_sqr,Fp_basis,Fp_add,Fp_add_mpz,Fp_add_ui,Fp_inv,Fp_neg\n");
    printf("%ld,",cost->Fp_mul);
    printf("%ld,",cost->Fp_mul_mpz);
    //printf("%ld,",cost->Fp_mul_ui);
    printf("%ld,",cost->Fp_sqr);
    //printf("%ld,",cost->Fp_basis);
    printf("%ld,",cost->Fp_add);
    //printf("%ld,",cost->Fp_add_mpz);
    printf("%ld,",cost->Fp_add_ui);
    printf("%ld,",cost->Fp_inv);
    printf("%ld",cost->Fp_neg);
    printf("\n");
}

/*----------------------------------------------------------------------------*/
//test
void test_plain_ate_pairing(){
    printf("====================================================================================\n");
    printf("Plain-ate pairing\n\n");
    EFp12 P,Q,s1P,s2P,s1Q,s2Q;
    EFp12_init(&P);
    EFp12_init(&Q);
    EFp12_init(&s1P);
    EFp12_init(&s2P);
    EFp12_init(&s1Q);
    EFp12_init(&s2Q);
    Fp12 Z,test1,test2,test3;
    Fp12_init(&Z);
    Fp12_init(&test1);
    Fp12_init(&test2);
    Fp12_init(&test3);
    mpz_t s12,s1,s2;
    mpz_init(s12);
    mpz_init(s1);
    mpz_init(s2);
    
    gmp_randstate_t state;
    gmp_randinit_default(state);
    gmp_randseed_ui(state,(unsigned long)time(NULL));
    mpz_urandomm(s1,state,order);
    mpz_urandomm(s2,state,order);
    mpz_mul(s12,s1,s2);
    mpz_mod(s12,s12,order);
    
    EFp12_generate_G1_karatsuba(&P);
    EFp12_printf(&P,"P\n");
    printf("\n\n");
    EFp12_generate_G2_karatsuba(&Q);
    EFp12_printf(&Q,"Q\n");
    printf("\n\n");
    EFp12_SCM_karatsuba(&s1P,&P,s1);
    EFp12_SCM_karatsuba(&s2P,&P,s2);
    EFp12_SCM_karatsuba(&s1Q,&Q,s1);
    EFp12_SCM_karatsuba(&s2Q,&Q,s2);
    
    printf("------------------------------------------------------------------------------------\n");
    printf("test1\n\n");
    printf("plain_ate(Q,P)^s1*s2\n");
    Plain_ate_pairing(&Z,&Q,&P);
    Fp12_pow_CVMA(&test1,&Z,s12);
    printf("Miller's Algo. (Plain-ate) : %.2f[ms]\n",MILLER_PLAINATE);
    printf("Final Exp. (easy)  : %.2f[ms]\n",FINALEXP_EASY);
    printf("Final Exp. (hard)  : %.2f[ms]\n",FINALEXP_HARD_OPT);
    printf("Final Exp. (total) : %.2f[ms]\n",FINALEXP_EASY+FINALEXP_HARD_OPT);
    Fp12_printf(&test1,"");
    printf("\n\n");
    
    printf("------------------------------------------------------------------------------------\n");
    printf("test2\n\n");
    printf("plain_ate([s2]Q,[s1]P)\n");
    Plain_ate_pairing(&test2,&s2Q,&s1P);
    printf("Miller's Algo. (Plain-ate) : %.2f[ms]\n",MILLER_PLAINATE);
    printf("Final Exp. (easy)  : %.2f[ms]\n",FINALEXP_EASY);
    printf("Final Exp. (hard)  : %.2f[ms]\n",FINALEXP_HARD_OPT);
    printf("Final Exp. (total) : %.2f[ms]\n",FINALEXP_EASY+FINALEXP_HARD_OPT);
    Fp12_printf(&test2,"");
    printf("\n\n");
    
    printf("------------------------------------------------------------------------------------\n");
    printf("test3\n\n");
    printf("plain_ate([s1]Q,[s2]P)\n");
    Plain_ate_pairing(&test3,&s1Q,&s2P);
    printf("Miller's Algo. (Plain-ate) : %.2f[ms]\n",MILLER_PLAINATE);
    printf("Final Exp. (easy)  : %.2f[ms]\n",FINALEXP_EASY);
    printf("Final Exp. (hard)  : %.2f[ms]\n",FINALEXP_HARD_OPT);
    printf("Final Exp. (total) : %.2f[ms]\n",FINALEXP_EASY+FINALEXP_HARD_OPT);
    Fp12_printf(&test3,"");
    printf("\n\n");
    
    printf("bilinear test\n");
    if(Fp12_cmp(&test1,&test2)==0 && Fp12_cmp(&test1,&test3)==0){
        printf("success\n\n");
    }else{
        printf("failed\n\n");
    }
    
    mpz_clear(s12);
    mpz_clear(s1);
    mpz_clear(s2);
    EFp12_clear(&P);
    EFp12_clear(&Q);
    EFp12_clear(&s1P);
    EFp12_clear(&s2P);
    EFp12_clear(&s1Q);
    EFp12_clear(&s2Q);
    Fp12_clear(&Z);
    Fp12_clear(&test1);
    Fp12_clear(&test2);
    Fp12_clear(&test3);
}

void test_opt_ate_pairing(){
    printf("====================================================================================\n");
    printf("Opt-ate pairing\n\n");
    EFp12 P,Q,s1P,s2P,s1Q,s2Q;
    EFp12_init(&P);
    EFp12_init(&Q);
    EFp12_init(&s1P);
    EFp12_init(&s2P);
    EFp12_init(&s1Q);
    EFp12_init(&s2Q);
    Fp12 Z,test1,test2,test3;
    Fp12_init(&Z);
    Fp12_init(&test1);
    Fp12_init(&test2);
    Fp12_init(&test3);
    mpz_t s12,s1,s2;
    mpz_init(s12);
    mpz_init(s1);
    mpz_init(s2);
    
    gmp_randstate_t state;
    gmp_randinit_default(state);
    gmp_randseed_ui(state,(unsigned long)time(NULL));
    mpz_urandomm(s1,state,order);
    mpz_urandomm(s2,state,order);
    mpz_mul(s12,s1,s2);
    mpz_mod(s12,s12,order);
    
    EFp12_generate_G1_karatsuba(&P);
    EFp12_printf(&P,"P\n");
    printf("\n\n");
    EFp12_generate_G2_karatsuba(&Q);
    EFp12_printf(&Q,"Q\n");
    printf("\n\n");
    EFp12_SCM_karatsuba(&s1P,&P,s1);
    EFp12_SCM_karatsuba(&s2P,&P,s2);
    EFp12_SCM_karatsuba(&s1Q,&Q,s1);
    EFp12_SCM_karatsuba(&s2Q,&Q,s2);
    
    printf("------------------------------------------------------------------------------------\n");
    printf("test1\n\n");
    printf("opt_ate(Q,P)^s1*s2\n");
    Opt_ate_pairing(&Z,&Q,&P);
    Fp12_pow_CVMA(&test1,&Z,s12);
    printf("Miller's Algo. (Opt-ate) : %.2f[ms]\n",MILLER_OPTATE);
    printf("Final Exp. (easy)  : %.2f[ms]\n",FINALEXP_EASY);
    printf("Final Exp. (hard)  : %.2f[ms]\n",FINALEXP_HARD_OPT);
    printf("Final Exp. (total) : %.2f[ms]\n",FINALEXP_EASY+FINALEXP_HARD_OPT);
    Fp12_printf(&test1,"");
    printf("\n\n");
    
    printf("------------------------------------------------------------------------------------\n");
    printf("test2\n\n");
    printf("opt_ate([s2]Q,[s1]P)\n");
    Opt_ate_pairing(&test2,&s2Q,&s1P);
    printf("Miller's Algo. (Opt-ate) : %.2f[ms]\n",MILLER_OPTATE);
    printf("Final Exp. (easy)  : %.2f[ms]\n",FINALEXP_EASY);
    printf("Final Exp. (hard)  : %.2f[ms]\n",FINALEXP_HARD_OPT);
    printf("Final Exp. (total) : %.2f[ms]\n",FINALEXP_EASY+FINALEXP_HARD_OPT);
    Fp12_printf(&test2,"");
    printf("\n\n");
    
    printf("------------------------------------------------------------------------------------\n");
    printf("test3\n\n");
    printf("opt_ate([s1]Q,[s2]P)\n");
    Opt_ate_pairing(&test3,&s1Q,&s2P);
    printf("Miller's Algo. (Opt-ate) : %.2f[ms]\n",MILLER_OPTATE);
    printf("Final Exp. (easy)  : %.2f[ms]\n",FINALEXP_EASY);
    printf("Final Exp. (hard)  : %.2f[ms]\n",FINALEXP_HARD_OPT);
    printf("Final Exp. (total) : %.2f[ms]\n",FINALEXP_EASY+FINALEXP_HARD_OPT);
    Fp12_printf(&test3,"");
    printf("\n\n");
    
    printf("bilinear test\n");
    if(Fp12_cmp(&test1,&test2)==0 && Fp12_cmp(&test1,&test3)==0){
        printf("success\n\n");
    }else{
        printf("failed\n\n");
    }
    
    mpz_clear(s12);
    mpz_clear(s1);
    mpz_clear(s2);
    EFp12_clear(&P);
    EFp12_clear(&Q);
    EFp12_clear(&s1P);
    EFp12_clear(&s2P);
    EFp12_clear(&s1Q);
    EFp12_clear(&s2Q);
    Fp12_clear(&Z);
    Fp12_clear(&test1);
    Fp12_clear(&test2);
    Fp12_clear(&test3);
}

void test_x_ate_pairing(){
    printf("====================================================================================\n");
    printf("X-ate pairing\n\n");
    EFp12 P,Q,s1P,s2P,s1Q,s2Q;
    EFp12_init(&P);
    EFp12_init(&Q);
    EFp12_init(&s1P);
    EFp12_init(&s2P);
    EFp12_init(&s1Q);
    EFp12_init(&s2Q);
    Fp12 Z,test1,test2,test3;
    Fp12_init(&Z);
    Fp12_init(&test1);
    Fp12_init(&test2);
    Fp12_init(&test3);
    mpz_t s12,s1,s2;
    mpz_init(s12);
    mpz_init(s1);
    mpz_init(s2);
    
    gmp_randstate_t state;
    gmp_randinit_default(state);
    gmp_randseed_ui(state,(unsigned long)time(NULL));
    mpz_urandomm(s1,state,order);
    mpz_urandomm(s2,state,order);
    mpz_mul(s12,s1,s2);
    mpz_mod(s12,s12,order);
    
    EFp12_generate_G1_karatsuba(&P);
    EFp12_printf(&P,"P\n");
    printf("\n\n");
    EFp12_generate_G2_karatsuba(&Q);
    EFp12_printf(&Q,"Q\n");
    printf("\n\n");
    EFp12_SCM_karatsuba(&s1P,&P,s1);
    EFp12_SCM_karatsuba(&s2P,&P,s2);
    EFp12_SCM_karatsuba(&s1Q,&Q,s1);
    EFp12_SCM_karatsuba(&s2Q,&Q,s2);
    
    printf("------------------------------------------------------------------------------------\n");
    printf("test1\n\n");
    printf("x_ate(Q,P)^s1*s2\n");
    X_ate_pairing(&Z,&Q,&P);
    Fp12_pow_CVMA(&test1,&Z,s12);
    printf("Miller's Algo. (X-ate) : %.2f[ms]\n",MILLER_XATE);
    printf("Final Exp. (easy)  : %.2f[ms]\n",FINALEXP_EASY);
    printf("Final Exp. (hard)  : %.2f[ms]\n",FINALEXP_HARD_OPT);
    printf("Final Exp. (total) : %.2f[ms]\n",FINALEXP_EASY+FINALEXP_HARD_OPT);
    Fp12_printf(&test1,"");
    printf("\n\n");
    
    printf("------------------------------------------------------------------------------------\n");
    printf("test2\n\n");
    printf("x_ate([s2]Q,[s1]P)\n");
    X_ate_pairing(&test2,&s2Q,&s1P);
    printf("Miller's Algo. (X-ate) : %.2f[ms]\n",MILLER_XATE);
    printf("Final Exp. (easy)  : %.2f[ms]\n",FINALEXP_EASY);
    printf("Final Exp. (hard)  : %.2f[ms]\n",FINALEXP_HARD_OPT);
    printf("Final Exp. (total) : %.2f[ms]\n",FINALEXP_EASY+FINALEXP_HARD_OPT);
    Fp12_printf(&test2,"");
    printf("\n\n");
    
    printf("------------------------------------------------------------------------------------\n");
    printf("test3\n\n");
    printf("x_ate([s1]Q,[s2]P)\n");
    X_ate_pairing(&test3,&s1Q,&s2P);
    printf("Miller's Algo. (X-ate) : %.2f[ms]\n",MILLER_XATE);
    printf("Final Exp. (easy)  : %.2f[ms]\n",FINALEXP_EASY);
    printf("Final Exp. (hard)  : %.2f[ms]\n",FINALEXP_HARD_OPT);
    printf("Final Exp. (total) : %.2f[ms]\n",FINALEXP_EASY+FINALEXP_HARD_OPT);
    Fp12_printf(&test3,"");
    printf("\n\n");
    
    printf("bilinear test\n");
    if(Fp12_cmp(&test1,&test2)==0 && Fp12_cmp(&test1,&test3)==0){
        printf("success\n\n");
    }else{
        printf("failed\n\n");
    }
    
    mpz_clear(s12);
    mpz_clear(s1);
    mpz_clear(s2);
    EFp12_clear(&P);
    EFp12_clear(&Q);
    EFp12_clear(&s1P);
    EFp12_clear(&s2P);
    EFp12_clear(&s1Q);
    EFp12_clear(&s2Q);
    Fp12_clear(&Z);
    Fp12_clear(&test1);
    Fp12_clear(&test2);
    Fp12_clear(&test3);
}

void test_G1_SCM(){
    printf("====================================================================================\n");
    printf("G1 SCM\n\n");
    EFp12 P,test1,test2,test3;
    EFp12_init(&P);
    EFp12_init(&test1);
    EFp12_init(&test2);
    EFp12_init(&test3);
    mpz_t scalar;
    mpz_init(scalar);
    
    gmp_randstate_t state;
    gmp_randinit_default(state);
    gmp_randseed_ui(state,(unsigned long)time(NULL));
    mpz_urandomm(scalar,state,order);
    printf("scalar : %d-bit\n",(int)mpz_sizeinbase(scalar,2));
    gmp_printf("%Zd\n",scalar);
    
    printf("generating rational point P in G1\n\n");
    EFp12_generate_G1_karatsuba(&P);
    
    printf("------------------------------------------------------------------------------------\n");
    printf("test1\n\n");
    printf("plain\n");
    EFp12_G1_SCM_plain_karatsuba(&test1,&P,scalar);
    printf("G1 SCM (plain) : %.2f[ms]\n",G1SCM_PLAIN);
    EFp12_printf(&test1,"");
    printf("\n\n");
    
    printf("------------------------------------------------------------------------------------\n");
    printf("test2\n\n");
    printf("2split\n");
    EFp12_G1_SCM_2split_karatsuba(&test2,&P,scalar);
    printf("G1 SCM (2split) : %.2f[ms]\n",G1SCM_2SPLIT);
    EFp12_printf(&test2,"");
    printf("\n\n");
    
    printf("------------------------------------------------------------------------------------\n");
    printf("test3\n\n");
    printf("2split-JSF\n");
    EFp12_G1_SCM_2split_JSF_karatsuba(&test3,&P,scalar);
    printf("G1 SCM (2split-JSF) : %.2f[ms]\n",G1SCM_2SPLIT_JSF);
    EFp12_printf(&test3,"");
    printf("\n\n");
    
    if(Fp12_cmp(&test1.x,&test2.x)==0 && Fp12_cmp(&test1.y,&test2.y)==0
       && Fp12_cmp(&test1.x,&test3.x)==0 && Fp12_cmp(&test1.y,&test3.y)==0){
        printf("success\n\n");
    }else{
        printf("failed\n\n");
    }
    
    mpz_clear(scalar);
    EFp12_clear(&P);
    EFp12_clear(&test1);
    EFp12_clear(&test2);
    EFp12_clear(&test3);
}

void test_G2_SCM(){
    printf("====================================================================================\n");
    printf("G2 SCM\n\n");
    EFp12 Q,test1,test2,test3,test4;
    EFp12_init(&Q);
    EFp12_init(&test1);
    EFp12_init(&test2);
    EFp12_init(&test3);
    EFp12_init(&test4);
    mpz_t scalar;
    mpz_init(scalar);
    
    gmp_randstate_t state;
    gmp_randinit_default(state);
    gmp_randseed_ui(state,(unsigned long)time(NULL));
    mpz_urandomm(scalar,state,order);
    printf("scalar : %d-bit\n",(int)mpz_sizeinbase(scalar,2));
    
    printf("generating rational point Q in G2\n\n");
    EFp12_generate_G2_karatsuba(&Q);
    
    printf("------------------------------------------------------------------------------------\n");
    printf("test1\n\n");
    printf("plain\n");
    EFp12_G2_SCM_plain_karatsuba(&test1,&Q,scalar);
    printf("G2 SCM (plain) : %.2f[ms]\n",G2SCM_PLAIN);
    EFp12_printf(&test1,"");
    printf("\n\n");
    
    printf("------------------------------------------------------------------------------------\n");
    printf("test2\n\n");
    printf("2split\n");
    EFp12_G2_SCM_2split_karatsuba(&test2,&Q,scalar);
    printf("G2 SCM (2split) : %.2f[ms]\n",G2SCM_2SPLIT);
    EFp12_printf(&test2,"");
    printf("\n\n");
    
    printf("------------------------------------------------------------------------------------\n");
    printf("test3\n\n");
    printf("2split-JSF\n");
    EFp12_G2_SCM_2split_JSF_karatsuba(&test3,&Q,scalar);
    printf("G2 SCM (2split-JSF) : %.2f[ms]\n",G2SCM_2SPLIT_JSF);
    EFp12_printf(&test3,"");
    printf("\n\n");
    
    printf("------------------------------------------------------------------------------------\n");
    printf("test4\n\n");
    printf("4split\n");
    EFp12_G2_SCM_4split_karatsuba(&test4,&Q,scalar);
    printf("G2 SCM (4split) : %.2f[ms]\n",G2SCM_4SPLIT);
    EFp12_printf(&test4,"");
    printf("\n\n");
    
    if(Fp12_cmp(&test1.x,&test2.x)==0 && Fp12_cmp(&test1.y,&test2.y)==0
       && Fp12_cmp(&test1.x,&test3.x)==0 && Fp12_cmp(&test1.y,&test3.y)==0
       && Fp12_cmp(&test1.x,&test4.x)==0 && Fp12_cmp(&test1.y,&test4.y)==0){
        printf("success\n\n");
    }else{
        printf("failed\n\n");
    }
    
    mpz_clear(scalar);
    EFp12_clear(&Q);
    EFp12_clear(&test1);
    EFp12_clear(&test2);
    EFp12_clear(&test3);
    EFp12_clear(&test4);
}

void test_G3_EXP(){
    printf("====================================================================================\n");
    printf("G3 Exp.\n\n");
    EFp12 P,Q;
    EFp12_init(&P);
    EFp12_init(&Q);
    Fp12 Z,test1,test2,test3,test4;
    Fp12_init(&Z);
    Fp12_init(&test1);
    Fp12_init(&test2);
    Fp12_init(&test3);
    Fp12_init(&test4);
    mpz_t scalar;
    mpz_init(scalar);
    
    gmp_randstate_t state;
    gmp_randinit_default(state);
    gmp_randseed_ui(state,(unsigned long)time(NULL));
    mpz_urandomm(scalar,state,order);
    printf("scalar : %d-bit\n",(int)mpz_sizeinbase(scalar,2));
    
    printf("generating rational point P in G1\n\n");
    EFp12_generate_G1_karatsuba(&P);
    printf("generating rational point Q in G2\n\n");
    EFp12_generate_G2_karatsuba(&Q);
    printf("x-ate(Q,P)\n");
    Opt_ate_pairing(&Z,&Q,&P);
    
    printf("------------------------------------------------------------------------------------\n");
    printf("test1\n\n");
    printf("plain\n");
    Fp12_G3_EXP_plain_CVMA(&test1,&Z,scalar);
    printf("G3 SCM (plain) : %.2f[ms]\n",G3EXP_PLAIN);
    Fp12_printf(&test1,"");
    printf("\n\n");
    
    printf("------------------------------------------------------------------------------------\n");
    printf("test2\n\n");
    printf("2split\n");
    Fp12_G3_EXP_2split_CVMA(&test2,&Z,scalar);
    printf("G3 SCM (2split) : %.2f[ms]\n",G3EXP_2SPLIT);
    Fp12_printf(&test2,"");
    printf("\n\n");
    
    printf("------------------------------------------------------------------------------------\n");
    printf("test3\n\n");
    printf("2split-JSF\n");
    Fp12_G3_EXP_2split_JSF_CVMA(&test3,&Z,scalar);
    printf("G3 SCM (2split-JSF) : %.2f[ms]\n",G3EXP_2SPLIT_JSF);
    Fp12_printf(&test3,"");
    printf("\n\n");
    
    printf("------------------------------------------------------------------------------------\n");
    printf("test4\n\n");
    printf("4split\n");
    Fp12_G3_EXP_4split_CVMA(&test4,&Z,scalar);
    printf("G3 SCM (4split) : %.2f[ms]\n",G3EXP_4SPLIT);
    Fp12_printf(&test4,"");
    printf("\n\n");
    
    if(Fp12_cmp(&test1,&test2)==0 && Fp12_cmp(&test1,&test3)==0 && Fp12_cmp(&test1,&test4)==0){
        printf("success\n\n");
    }else{
        printf("failed\n\n");
    }
    
    mpz_clear(scalar);
    EFp12_clear(&P);
    EFp12_clear(&Q);
    Fp12_clear(&Z);
    Fp12_clear(&test1);
    Fp12_clear(&test2);
    Fp12_clear(&test3);
    Fp12_clear(&test4);
}

void test_basis_conversion(){
    printf("====================================================================================\n");
    printf("Basis conversion\n\n");
    Fp12 tmp,test1,test2;
    Fp12_init(&tmp);
    Fp12_init(&test1);
    Fp12_init(&test2);
    mpz_t exp;
    mpz_init(exp);
    
    gmp_randstate_t state;
    gmp_randinit_default(state);
    gmp_randseed_ui(state,(unsigned long)time(NULL));
    Fp12_set_random(&tmp,state);
    mpz_urandomm(exp,state,order);
    
    printf("test1\n");
    Fp12_pow_karatsuba(&test1,&tmp,exp);
    Fp12_printf(&test1,"");
    printf("\n\n");
    
    printf("test2\n");
    karatsuba_to_CVMA(&test2,&tmp);
    Fp12_pow_CVMA(&test2,&test2,exp);
    CVMA_to_karatsuba(&test2,&test2);
    Fp12_printf(&test2,"");
    printf("\n\n");
    
    if(Fp12_cmp(&test1,&test2)==0){
        printf("success\n");
    }else{
        printf("failed\n");
    }
    
    Fp12_clear(&tmp);
    Fp12_clear(&test1);
    Fp12_clear(&test2);
    mpz_clear(exp);
}

void computation_time(){
    printf("====================================================================================\n");
    printf("Computation time\n\n");
    int count;
    float AVE_BASIS_CONVERSION_CTOK=0,AVE_BASIS_CONVERSION_KTOC=0;
    float AVE_MILLER_PLAINATE=0,AVE_MILLER_XATE=0,AVE_MILLER_OPTATE=0;
    float AVE_FINALEXP_EASY=0,AVE_FINALEXP_HARD_PLAIN=0,AVE_FINALEXP_HARD_OPT=0;
    float AVE_G1SCM_PLAIN=0,AVE_G1SCM_2SPLIT=0,AVE_G1SCM_2SPLIT_JSF=0;
    float AVE_G2SCM_PLAIN=0,AVE_G2SCM_2SPLIT=0,AVE_G2SCM_2SPLIT_JSF=0,AVE_G2SCM_4SPLIT=0;
    float AVE_G3EXP_PLAIN=0,AVE_G3EXP_2SPLIT=0,AVE_G3EXP_2SPLIT_JSF=0,AVE_G3EXP_4SPLIT=0;
    EFp12 P,Q,tmp_P,tmp_Q;
    EFp12_init(&P);
    EFp12_init(&Q);
    EFp12_init(&tmp_P);
    EFp12_init(&tmp_Q);
    Fp12 Z,tmp_Z;
    Fp12_init(&Z);
    Fp12_init(&tmp_Z);
    gmp_randstate_t state;
    gmp_randinit_default(state);
    gmp_randseed_ui(state,(unsigned long)time(NULL));
    mpz_t scalar;
    mpz_init(scalar);
    
    printf("generating rational point P in G1\n\n");
    EFp12_generate_G1_karatsuba(&P);
    EFp12_printf(&P,"P\n");
    printf("\n\n");
    
    printf("generating rational point Q in G2\n\n");
    EFp12_generate_G2_karatsuba(&Q);
    EFp12_printf(&Q,"Q\n");
    printf("\n\n");
    
    printf("Computing...\n");
    count=0;
    while(count<100){
        printf("%d\n",count);
        
        //Miller's Alg. (T'_1)
        Plain_ate_pairing(&Z,&Q,&P);
        AVE_MILLER_PLAINATE=AVE_MILLER_PLAINATE+MILLER_PLAINATE;
        Opt_ate_pairing(&Z,&Q,&P);
        AVE_MILLER_OPTATE=AVE_MILLER_OPTATE+MILLER_OPTATE;
        X_ate_pairing(&Z,&Q,&P);
        AVE_MILLER_XATE=AVE_MILLER_XATE+MILLER_XATE;
        
        //Basis conversion (T'_1 to T'_2)
        karatsuba_to_CVMA(&Z,&Z);
        AVE_BASIS_CONVERSION_KTOC=AVE_BASIS_CONVERSION_KTOC+BASIS_CONVERSION_KTOC;
        
        //Final exp (T'_2)
        Final_exp_easy_CVMA(&Z,&Z);
        Final_exp_hard_optimal_CVMA(&Z,&Z);
        AVE_FINALEXP_EASY=AVE_FINALEXP_EASY+FINALEXP_EASY;
        AVE_FINALEXP_HARD_OPT=AVE_FINALEXP_HARD_OPT+FINALEXP_HARD_OPT;
        
        mpz_urandomm(scalar,state,order);
        
        //G1_SCM (T'_1)
        EFp12_G1_SCM_plain_karatsuba(&tmp_P,&P,scalar);
        EFp12_G1_SCM_2split_karatsuba(&tmp_P,&P,scalar);
        EFp12_G1_SCM_2split_JSF_karatsuba(&tmp_P,&P,scalar);
        AVE_G1SCM_PLAIN=AVE_G1SCM_PLAIN+G1SCM_PLAIN;
        AVE_G1SCM_2SPLIT=AVE_G1SCM_2SPLIT+G1SCM_2SPLIT;
        AVE_G1SCM_2SPLIT_JSF=AVE_G1SCM_2SPLIT_JSF+G1SCM_2SPLIT_JSF;
        
        //G2_SCM (T'_1)
        EFp12_G2_SCM_plain_karatsuba(&tmp_Q,&Q,scalar);
        EFp12_G2_SCM_2split_karatsuba(&tmp_Q,&Q,scalar);
        EFp12_G2_SCM_2split_JSF_karatsuba(&tmp_Q,&Q,scalar);
        EFp12_G2_SCM_4split_karatsuba(&tmp_Q,&Q,scalar);
        AVE_G2SCM_PLAIN=AVE_G2SCM_PLAIN+G2SCM_PLAIN;
        AVE_G2SCM_2SPLIT=AVE_G2SCM_2SPLIT+G2SCM_2SPLIT;
        AVE_G2SCM_2SPLIT_JSF=AVE_G2SCM_2SPLIT_JSF+G2SCM_2SPLIT_JSF;
        AVE_G2SCM_4SPLIT=AVE_G2SCM_4SPLIT+G2SCM_4SPLIT;
        
        //G3_EXP (T'_2)
        Fp12_G3_EXP_plain_CVMA(&tmp_Z,&Z,scalar);
        Fp12_G3_EXP_2split_CVMA(&tmp_Z,&Z,scalar);
        Fp12_G3_EXP_2split_JSF_CVMA(&tmp_Z,&Z,scalar);
        Fp12_G3_EXP_4split_CVMA(&tmp_Z,&Z,scalar);
        AVE_G3EXP_PLAIN=AVE_G3EXP_PLAIN+G3EXP_PLAIN;
        AVE_G3EXP_2SPLIT=AVE_G3EXP_2SPLIT+G3EXP_2SPLIT;
        AVE_G3EXP_2SPLIT_JSF=AVE_G3EXP_2SPLIT_JSF+G3EXP_2SPLIT_JSF;
        AVE_G3EXP_4SPLIT=AVE_G3EXP_4SPLIT+G3EXP_4SPLIT;
        
        //Basis conversion (T'_2 tp T'_1)
        CVMA_to_karatsuba(&Z,&Z);
        AVE_BASIS_CONVERSION_CTOK=AVE_BASIS_CONVERSION_CTOK+BASIS_CONVERSION_CTOK;
        
        count++;
    }
    printf("\n");
    
    printf("Miller's Algo. (plain-ate) : %.2f[ms]\n",AVE_MILLER_PLAINATE/100);
    printf("Miller's Algo. (opt-ate)   : %.2f[ms]\n",AVE_MILLER_OPTATE/100);
    printf("Miller's Algo. (x-ate)     : %.2f[ms]\n",AVE_MILLER_XATE/100);
    printf("\n");
    printf("Final Exp. Plain   (easy) : %.2f[ms]\n",AVE_FINALEXP_EASY/100);
    printf("Final Exp. Optimal (hard) : %.2f[ms]\n",AVE_FINALEXP_HARD_OPT/100);
    printf("Karatsuba to CVMA         : %.2f[ms]\n",AVE_BASIS_CONVERSION_KTOC/100);
    printf("CVMA to karatsuba         : %.2f[ms]\n",AVE_BASIS_CONVERSION_CTOK/100);
    printf("G1 SCM (plain)            : %.2f[ms]\n",AVE_G1SCM_PLAIN/100);
    printf("G1 SCM (2split)           : %.2f[ms]\n",AVE_G1SCM_2SPLIT/100);
    printf("G1 SCM (2split-JSF)       : %.2f[ms]\n",AVE_G1SCM_2SPLIT_JSF/100);
    printf("\n");
    printf("G2 SCM (plain)            : %.2f[ms]\n",AVE_G2SCM_PLAIN/100);
    printf("G2 SCM (2split)           : %.2f[ms]\n",AVE_G2SCM_2SPLIT/100);
    printf("G2 SCM (2split-JSF)       : %.2f[ms]\n",AVE_G2SCM_2SPLIT_JSF/100);
    printf("G2 SCM (4split)           : %.2f[ms]\n",AVE_G2SCM_4SPLIT/100);
    printf("\n");
    printf("G3 EXP (plain)            : %.2f[ms]\n",AVE_G3EXP_PLAIN/100);
    printf("G3 EXP (2split)           : %.2f[ms]\n",AVE_G3EXP_2SPLIT/100);
    printf("G3 EXP (2split-JSF)       : %.2f[ms]\n",AVE_G3EXP_2SPLIT_JSF/100);
    printf("G3 EXP (4split)           : %.2f[ms]\n",AVE_G3EXP_4SPLIT/100);
    printf("\n");
    
    EFp12_clear(&P);
    EFp12_clear(&Q);
    EFp12_clear(&tmp_P);
    EFp12_clear(&tmp_Q);
    Fp12_clear(&Z);
    Fp12_clear(&tmp_Z);
    mpz_clear(scalar);
}

void computation_cost(){
    printf("====================================================================================\n");
    printf("Pairing operations\n\n");
    EFp12 P,Q,tmp_P,tmp_Q;
    EFp12_init(&P);
    EFp12_init(&Q);
    EFp12_init(&tmp_P);
    EFp12_init(&tmp_Q);
    Fp12 Z,tmp_Z;
    Fp12_init(&Z);
    Fp12_init(&tmp_Z);
    mpz_t scalar;
    mpz_init(scalar);
    
    //gmp_randstate_t state;
    //gmp_randinit_default(state);
    //gmp_randseed_ui(state,(unsigned long)time(NULL));
    //mpz_urandomm(scalar,state,order);
    
    mpz_set_str(scalar,"6433987872172996767370742739789192855206084819751810317146156999801714004530205817441825065255166612237926117694268872766909697908042480267",10);
    printf("scalar : %d-bit\n",(int)mpz_sizeinbase(scalar,2));
    
    printf("generating rational point P in G1\n\n");
    EFp12_generate_G1_karatsuba(&P);
    printf("generating rational point Q in G2\n\n");
    EFp12_generate_G2_karatsuba(&Q);
    
    //Miller's Alg. (T'_1)
    printf("Miller's Alg. (T'_1)\n");
    Init_Fp_Cost(&Fp_cost);
    Miller_algo_for_plain_ate_karatsuba(&Z,&P,&Q);
    Print_Fp_Cost(&Fp_cost,"Plain-ate\n");
    printf("\n");
    Init_Fp_Cost(&Fp_cost);
    Miller_algo_for_opt_ate_karatsuba(&Z,&P,&Q);
    Print_Fp_Cost(&Fp_cost,"Opt-ate\n");
    printf("\n");
    Init_Fp_Cost(&Fp_cost);
    Miller_algo_for_x_ate_karatsuba(&Z,&P,&Q);
    Print_Fp_Cost(&Fp_cost,"X-ate\n");
    printf("\n");
    
    //Basis conversion (T'_1 to T'_2)
    printf("Basis conversion (T'_1 to T'_2)\n");
    Init_Fp_Cost(&Fp_cost);
    karatsuba_to_CVMA(&Z,&Z);
    Print_Fp_Cost(&Fp_cost,"");
    printf("\n");

    //Final Exp. (T'_2)
    printf("Final Exp. (T'_2)\n");
    Init_Fp_Cost(&Fp_cost);
    Final_exp_easy_CVMA(&Z,&Z);
    Print_Fp_Cost(&Fp_cost,"Final exp (easy)\n");
    Init_Fp_Cost(&Fp_cost);
    Final_exp_hard_plain_CVMA(&Z,&Z);
    Print_Fp_Cost(&Fp_cost,"Final exp Plain (hard)\n");
    Init_Fp_Cost(&Fp_cost);
    Final_exp_hard_optimal_CVMA(&Z,&Z);
    Print_Fp_Cost(&Fp_cost,"Final exp Opt (hard)\n");
    printf("\n");
    
    //G1_SCM
    printf("G1 SCM (T'_1)\n");
    Init_Fp_Cost(&Fp_cost);
    EFp12_G1_SCM_plain_karatsuba(&tmp_P,&P,scalar);
    Print_Fp_Cost(&Fp_cost,"plain\n");
    Init_Fp_Cost(&Fp_cost);
    EFp12_G1_SCM_2split_karatsuba(&tmp_P,&P,scalar);
    Print_Fp_Cost(&Fp_cost,"2split\n");
    Init_Fp_Cost(&Fp_cost);
    EFp12_G1_SCM_2split_JSF_karatsuba(&tmp_P,&P,scalar);
    Print_Fp_Cost(&Fp_cost,"2split-JSF\n");
    printf("\n");
    
    //G2_SCM
    printf("G2 SCM (T'_2)\n");
    Init_Fp_Cost(&Fp_cost);
    EFp12_G2_SCM_plain_karatsuba(&tmp_Q,&Q,scalar);
    Print_Fp_Cost(&Fp_cost,"plain\n");
    Init_Fp_Cost(&Fp_cost);
    EFp12_G2_SCM_2split_karatsuba(&tmp_Q,&Q,scalar);
    Print_Fp_Cost(&Fp_cost,"2split\n");
    Init_Fp_Cost(&Fp_cost);
    EFp12_G2_SCM_2split_JSF_karatsuba(&tmp_Q,&Q,scalar);
    Print_Fp_Cost(&Fp_cost,"2split-JSF\n");
    Init_Fp_Cost(&Fp_cost);
    EFp12_G2_SCM_4split_karatsuba(&tmp_Q,&Q,scalar);
    Print_Fp_Cost(&Fp_cost,"4split\n");
    printf("\n");
    
    //G3_EXP
    printf("G3 EXP (T'_3)\n");
    Init_Fp_Cost(&Fp_cost);
    Fp12_G3_EXP_plain_CVMA(&tmp_Z,&Z,scalar);
    Print_Fp_Cost(&Fp_cost,"plain\n");
    Init_Fp_Cost(&Fp_cost);
    Fp12_G3_EXP_2split_CVMA(&tmp_Z,&Z,scalar);
    Print_Fp_Cost(&Fp_cost,"2split\n");
    Init_Fp_Cost(&Fp_cost);
    Fp12_G3_EXP_2split_JSF_CVMA(&tmp_Z,&Z,scalar);
    Print_Fp_Cost(&Fp_cost,"2split-JSF\n");
    Init_Fp_Cost(&Fp_cost);
    Fp12_G3_EXP_4split_CVMA(&tmp_Z,&Z,scalar);
    Print_Fp_Cost(&Fp_cost,"4split\n");
    printf("\n");
    
    //Basis conversion (T'_2 to T'_1)
    printf("Basis conversion (T'_2 to T'_1)\n");
    Init_Fp_Cost(&Fp_cost);
    CVMA_to_karatsuba(&Z,&Z);
    Print_Fp_Cost(&Fp_cost,"");
    printf("\n");
    
    EFp12_clear(&P);
    EFp12_clear(&Q);
    EFp12_clear(&tmp_P);
    EFp12_clear(&tmp_Q);
    Fp12_clear(&Z);
    Fp12_clear(&tmp_Z);
    mpz_clear(scalar);
}

