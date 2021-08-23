/********************************************************************************************
* Supersingular Isogeny Key Encapsulation Library
*
* Abstract: ephemeral supersingular isogeny Diffie-Hellman key exchange (SIDH)
*********************************************************************************************/ 

#include "random/random.h"
#include <stdio.h>

static void init_basis(digit_t *gen, f2elm_t XP, f2elm_t XQ, f2elm_t XR)
{ // Initialization of basis points
    
    fpcopy(gen,                  XP[0]);
    fpcopy(gen +   NWORDS_FIELD, XP[1]);
    fpcopy(gen + 2*NWORDS_FIELD, XQ[0]);
    fpcopy(gen + 3*NWORDS_FIELD, XQ[1]);
    fpcopy(gen + 4*NWORDS_FIELD, XR[0]);
    fpcopy(gen + 5*NWORDS_FIELD, XR[1]);
}


void random_mod_order_A(unsigned char* random_digits)
{  // Generation of Alice's secret key  
   // Outputs random value in [0, 2^eA - 1]

    randombytes(random_digits, SECRETKEY_A_BYTES);
    random_digits[SECRETKEY_A_BYTES-1] &= MASK_ALICE;    // Masking last byte 
}


void random_mod_order_B(unsigned char* random_digits)
{  // Generation of Bob's secret key  
   // Outputs random value in [0, 2^Floor(Log(2, oB)) - 1]

    randombytes(random_digits, SECRETKEY_B_BYTES);
    random_digits[SECRETKEY_B_BYTES-1] &= MASK_BOB;     // Masking last byte 
}




// Alice using 3-isogeny
int EphemeralKeyGeneration_A(const unsigned char* PrivateKeyA, unsigned char* PublicKeyA)
{ // Alice's ephemeral public key generation
  // Input:  a private key PrivateKeyA in the range [0, 2^Floor(Log(2,oA)) - 1]. 
  // Output: the public key PublicKeyA consisting of 3 elements in GF(p^2) which are encoded by removing leading 0 bytes.
    point_proj_t R, phiP = {0}, phiQ = {0}, phiR = {0}, pts[MAX_INT_POINTS_ALICE];
    f2elm_t XPA, XQA, XRA, coeff[3], A24plus = {0}, A24minus = {0}, A = {0};
    unsigned int i, row, m, index = 0, pts_index[MAX_INT_POINTS_ALICE], npts = 0, ii = 0;
    digit_t SecretKeyA[NWORDS_ORDER] = {0};

    // Initialize basis points
    init_basis((digit_t*)A_gen, XPA, XQA, XRA);
    init_basis((digit_t*)B_gen, phiP->X, phiQ->X, phiR->X);
    fpcopy((digit_t*)&Montgomery_one, (phiP->Z)[0]);
    fpcopy((digit_t*)&Montgomery_one, (phiQ->Z)[0]);
    fpcopy((digit_t*)&Montgomery_one, (phiR->Z)[0]);

    // Initialize constants: A24minus = A-2C, A24plus = A+2C, where A=6, C=1
    fpcopy((digit_t*)&Montgomery_one, A24plus[0]);
    mp2_add(A24plus, A24plus, A24plus);
    mp2_add(A24plus, A24plus, A24minus);
    mp2_add(A24plus, A24minus, A);
    mp2_add(A24minus, A24minus, A24plus);

    // Retrieve kernel point
    decode_to_digits(PrivateKeyA, SecretKeyA, SECRETKEY_A_BYTES, NWORDS_ORDER);
    LADDER3PT(XPA, XQA, XRA, SecretKeyA, ALICE, R, A);
    
    // Traverse tree
    index = 0;  
    for (row = 1; row < MAX_Alice; row++) {
        while (index < MAX_Alice-row) {
            fp2copy(R->X, pts[npts]->X);
            fp2copy(R->Z, pts[npts]->Z);
            pts_index[npts] = index;
            npts += 1;
            m = strat_Alice_Huff[MAX_Alice-index-row];
            xTPLe(R, R, A24minus, A24plus, (int)m);
            index += m;
        } 
        get_3_isog(R, A24minus, A24plus, coeff);

        for (i = 0; i < npts; i++) {
            eval_3_isog(pts[i], coeff);
        }     
        eval_3_isog(phiP, coeff);
        eval_3_isog(phiQ, coeff);
        eval_3_isog(phiR, coeff);

        fp2copy(pts[npts-1]->X, R->X); 
        fp2copy(pts[npts-1]->Z, R->Z);
        index = pts_index[npts-1];
        npts -= 1;
    }
    
    get_3_isog(R, A24minus, A24plus, coeff);
    eval_3_isog(phiP, coeff);
    eval_3_isog(phiQ, coeff);
    eval_3_isog(phiR, coeff);


    inv_3_way(phiP->Z, phiQ->Z, phiR->Z);
    fp2mul_mont(phiP->X, phiP->Z, phiP->X);
    fp2mul_mont(phiQ->X, phiQ->Z, phiQ->X);
    fp2mul_mont(phiR->X, phiR->Z, phiR->X);


    // Format public key
    fp2_encode(phiP->X, PublicKeyA);
    fp2_encode(phiQ->X, PublicKeyA + FP2_ENCODED_BYTES);
    fp2_encode(phiR->X, PublicKeyA + 2*FP2_ENCODED_BYTES);

    return 0;
}

// Bob using 5-isogney
int EphemeralKeyGeneration_B(const unsigned char* PrivateKeyB, unsigned char* PublicKeyB)
{ // Bob's ephemeral public key generation
  // Input:  a private key PrivateKeyB in the range [0, 2^Floor(Log(2,oB)) - 1]. 
  // Output: the public key PublicKeyB consisting of 3 elements in GF(p^2) which are encoded by removing leading 0 bytes.
    point_proj_t R, R2, phiP = {0}, phiQ = {0}, phiR = {0}, pts[MAX_INT_POINTS_BOB], P2={0};
    f2elm_t XPB, XQB, XRB, coeff[3], A24plus = {0}, C24 = {0}, A = {0};
    unsigned int i, row, m, index = 0, pts_index[MAX_INT_POINTS_BOB], npts = 0, ii = 0;
    digit_t SecretKeyB[NWORDS_ORDER] = {0};

    f2elm_t t0, t1, t2, atmp, ctmp;
    // Initialize basis points
    init_basis((digit_t*)B_gen, XPB, XQB, XRB);
    init_basis((digit_t*)A_gen, phiP->X, phiQ->X, phiR->X);
    fpcopy((digit_t*)&Montgomery_one, (phiP->Z)[0]);
    fpcopy((digit_t*)&Montgomery_one, (phiQ->Z)[0]);
    fpcopy((digit_t*)&Montgomery_one, (phiR->Z)[0]);
    fpcopy((digit_t*)&Montgomery_one, (P2->Z)[0]);
    fpcopy((digit_t*)&Mont_P2, (P2->X)[0]);
    // Initialize constants: A24plus = A+2C, C24 = 4C, where A=6, C=1
    fpcopy((digit_t*)&Montgomery_one, A24plus[0]); //1
    mp2_add(A24plus, A24plus, A24plus); //2
    mp2_add(A24plus, A24plus, C24); //4
    mp2_add(A24plus, C24, A); //6
    mp2_add(C24, C24, A24plus); //

    // Retrieve kernel point
    decode_to_digits(PrivateKeyB, SecretKeyB, SECRETKEY_B_BYTES, NWORDS_ORDER);
    LADDER3PT(XPB, XQB, XRB, SecretKeyB, BOB, R, A);       


    // Traverse tree
    index = 0;        
    for (row = 1; row < MAX_Bob; row++) {
        while (index < MAX_Bob-row) {
            fp2copy(R->X, pts[npts]->X);
            fp2copy(R->Z, pts[npts]->Z);
            pts_index[npts] = index;
            npts += 1;
            m = strat_Bob_Huff[MAX_Bob-index-row];
            x5Pe(R, R, A24plus, C24, (int)(m));
            index += m;
        }
        // kernel points R, [2]R  

        xDBL(R, R2, A24plus, C24);



        for (i = 0; i < npts; i++) {
            eval_5_isog(pts[i], R, R2);
        }
        eval_5_isog(phiP, R, R2);
        eval_5_isog(phiQ, R, R2);
        eval_5_isog(phiR, R, R2);


        // coefficient : evaluate isogeny at 2-torsion, recover
        eval_5_isog(P2, R, R2);
        get_5_isog(P2, A24plus, C24);  

        fp2copy(pts[npts-1]->X, R->X); 
        fp2copy(pts[npts-1]->Z, R->Z);
        index = pts_index[npts-1];
        npts -= 1;
    }

    xDBL(R, R2, A24plus, C24);
    eval_5_isog(phiP, R, R2);
    eval_5_isog(phiQ, R, R2);
    eval_5_isog(phiR, R, R2);

    // last step doesnt need to compute the coefficeint


    inv_3_way(phiP->Z, phiQ->Z, phiR->Z);
    fp2mul_mont(phiP->X, phiP->Z, phiP->X);
    fp2mul_mont(phiQ->X, phiQ->Z, phiQ->X);
    fp2mul_mont(phiR->X, phiR->Z, phiR->X);


    // Format public key                   
    fp2_encode(phiP->X, PublicKeyB);
    fp2_encode(phiQ->X, PublicKeyB + FP2_ENCODED_BYTES);
    fp2_encode(phiR->X, PublicKeyB + 2*FP2_ENCODED_BYTES);

    return 0;
}



int EphemeralSecretAgreement_A(const unsigned char* PrivateKeyA, const unsigned char* PublicKeyB, unsigned char* SharedSecretA)
{ // Alice's ephemeral shared secret computation
  // It produces a shared secret key SharedSecretA using her secret key PrivateKeyA and Bob's public key PublicKeyB
  // Inputs: Alice's PrivateKeyA is an integer in the range [0, 2^Floor(Log(2,oA)) - 1]. 
  //         Bob's PublicKeyB consists of 3 elements in GF(p^2) encoded by removing leading 0 bytes.
  // Output: a shared secret SharedSecretA that consists of one element in GF(p^2) encoded by removing leading 0 bytes.  
    point_proj_t R, pts[MAX_INT_POINTS_ALICE];
    f2elm_t coeff[3], PKA[3], jinv;
    f2elm_t A24plus = {0}, A24minus = {0}, A = {0};
    unsigned int i, row, m, index = 0, pts_index[MAX_INT_POINTS_ALICE], npts = 0, ii = 0;
    digit_t SecretKeyA[NWORDS_ORDER] = {0};
      
    // Initialize images of Alice's basis
    fp2_decode(PublicKeyB, PKA[0]);
    fp2_decode(PublicKeyB + FP2_ENCODED_BYTES, PKA[1]);
    fp2_decode(PublicKeyB + 2*FP2_ENCODED_BYTES, PKA[2]);

    // Initialize constants: A24plus = A+2C, A24minus = A-2C, where C=1
    get_A(PKA[0], PKA[1], PKA[2], A);
    mp_add((digit_t*)&Montgomery_one, (digit_t*)&Montgomery_one, A24minus[0], NWORDS_FIELD);
    mp2_add(A, A24minus, A24plus);
    mp2_sub_p2(A, A24minus, A24minus);

    // Retrieve kernel point
    decode_to_digits(PrivateKeyA, SecretKeyA, SECRETKEY_A_BYTES, NWORDS_ORDER);
    LADDER3PT(PKA[0], PKA[1], PKA[2], SecretKeyA, ALICE, R, A);


    
    // Traverse tree
    index = 0;  
    for (row = 1; row < MAX_Alice; row++) {
        while (index < MAX_Alice-row) {
            fp2copy(R->X, pts[npts]->X);
            fp2copy(R->Z, pts[npts]->Z);
            pts_index[npts] = index;
            npts += 1;
            m = strat_Alice_Huff[MAX_Alice-index-row];
            xTPLe(R, R, A24minus, A24plus, (int)m);
            index += m;
        }
        get_3_isog(R, A24minus, A24plus, coeff);

        for (i = 0; i < npts; i++) {
            eval_3_isog(pts[i], coeff);
        } 

        fp2copy(pts[npts-1]->X, R->X); 
        fp2copy(pts[npts-1]->Z, R->Z);
        index = pts_index[npts-1];
        npts -= 1;
    }
     
    get_3_isog(R, A24minus, A24plus, coeff);  



    fp2add(A24plus, A24minus, A);                 
    fp2add(A, A, A);
    fp2sub(A24plus, A24minus, A24plus);                   
    j_inv(A, A24plus, jinv);
    fp2_encode(jinv, SharedSecretA);    // Format shared secret

    return 0;
}



int EphemeralSecretAgreement_B(const unsigned char* PrivateKeyB, const unsigned char* PublicKeyA, unsigned char* SharedSecretB)
{ // Bob's ephemeral shared secret computation
  // It produces a shared secret key SharedSecretB using his secret key PrivateKeyB and Alice's public key PublicKeyA
  // Inputs: Bob's PrivateKeyA is an integer in the range [0, oB-1]. 
  //         Alice's PublicKeyA consists of 3 elements in GF(p^2) encoded by removing leading 0 bytes.
  // Output: a shared secret SharedSecretB that consists of one element in GF(p^2) encoded by removing leading 0 bytes.  
    point_proj_t R, R2, pts[MAX_INT_POINTS_BOB], P2={0}, t0={0};
    f2elm_t coeff[3], PKA[3], jinv;
    f2elm_t A24plus = {0}, C24 = {0}, A = {0};
    unsigned int i, row, m, index = 0, pts_index[MAX_INT_POINTS_BOB], npts = 0, ii = 0;
    digit_t SecretKeyB[NWORDS_ORDER] = {0};
      
    // Initialize images of Bob's basis
    fp2_decode(PublicKeyA, PKA[0]);
    fp2_decode(PublicKeyA + FP2_ENCODED_BYTES, PKA[1]);
    fp2_decode(PublicKeyA + 2*FP2_ENCODED_BYTES, PKA[2]);
    fpcopy((digit_t*)&Montgomery_one, (t0->Z)[0]);
    // Initialize constants: A24plus = A+2C, C24 = 4C, where C=1
    get_A(PKA[0], PKA[1], PKA[2], A);
    get_2torsion(A, P2);



    mp_add((digit_t*)&Montgomery_one, (digit_t*)&Montgomery_one, C24[0], NWORDS_FIELD);
    mp2_add(A, C24, A24plus);
    mp_add(C24[0], C24[0], C24[0], NWORDS_FIELD);

    // Retrieve kernel point
    decode_to_digits(PrivateKeyB, SecretKeyB, SECRETKEY_B_BYTES, NWORDS_ORDER);
    LADDER3PT(PKA[0], PKA[1], PKA[2], SecretKeyB, BOB, R, A);    


    // Traverse tree
    index = 0;        
    for (row = 1; row < MAX_Bob; row++) {
        while (index < MAX_Bob-row) {
            fp2copy(R->X, pts[npts]->X);
            fp2copy(R->Z, pts[npts]->Z);
            pts_index[npts] = index;
            npts += 1;
            m = strat_Bob_Huff[MAX_Bob-index-row];
            x5Pe(R, R, A24plus, C24, (int)(m));
            index += m;
        }
        // kernel points R, [2]R      
        xDBL(R, R2, A24plus, C24);      

        for (i = 0; i < npts; i++) {
            eval_5_isog(pts[i], R, R2);

        }

        // coeffiecient
        eval_5_isog(P2, R, R2);
        get_5_isog(P2, A24plus, C24);  

        fp2copy(pts[npts-1]->X, R->X); 
        fp2copy(pts[npts-1]->Z, R->Z);
        index = pts_index[npts-1];
        npts -= 1;
    }

    xDBL(R, R2, A24plus, C24);    
    eval_5_isog(P2, R, R2);
    get_5_isog(P2, A24plus, C24);  

 


    mp2_add(A24plus, A24plus, A24plus);                                                
    fp2sub(A24plus, C24, A24plus); 
    fp2add(A24plus, A24plus, A24plus);                    
    j_inv(A24plus, C24, jinv);
    fp2_encode(jinv, SharedSecretB);    // Format shared secret

    return 0;
}

// HUFF SIDH
int EphemeralKeyGeneration_A_Huff(const unsigned char* PrivateKeyA, unsigned char* PublicKeyA)
{ // Alice's ephemeral public key generation
  // Input:  a private key PrivateKeyA in the range [0, 2^Floor(Log(2,oA)) - 1]. 
  // Output: the public key PublicKeyA consisting of 3 elements in GF(p^2) which are encoded by removing leading 0 bytes.
    point_proj_t R, phiP = {0}, phiQ = {0}, phiR = {0}, pts[MAX_INT_POINTS_ALICE];
    f2elm_t XPA, XQA, XRA, coeff[3], A24plus = {0}, A24minus = {0}, A = {0};
    unsigned int i, row, m, index = 0, pts_index[MAX_INT_POINTS_ALICE], npts = 0, ii = 0;
    digit_t SecretKeyA[NWORDS_ORDER] = {0};

    // Initialize basis points
    init_basis((digit_t*)A_gen_Huff, XPA, XQA, XRA);
    init_basis((digit_t*)B_gen_Huff, phiP->X, phiQ->X, phiR->X);
    fpcopy((digit_t*)&Montgomery_one, (phiP->Z)[0]);
    fpcopy((digit_t*)&Montgomery_one, (phiQ->Z)[0]);
    fpcopy((digit_t*)&Montgomery_one, (phiR->Z)[0]);

    // Initialize constants: A=c+1/c-2 where c=3+sqrt{8} in GF(p^2) A24minus = (C-D)^2, A24plus = (C+D)^2
    fpcopy((digit_t*)&Montgomery_one, A[0]); // A =1
    fp2add(A, A, A); // A =2
    fp2add(A, A, A); // A =4

    fp2copy(A, A24minus); //4
    fp2copy(A, A24plus); //4
    fp2add(A24plus, A24plus, A24plus); // 8


    // Retrieve kernel point
    decode_to_digits(PrivateKeyA, SecretKeyA, SECRETKEY_A_BYTES, NWORDS_ORDER);
    LADDER3PT_Huff(XPA, XQA, XRA, SecretKeyA, ALICE, R, A);
    
    // Traverse tree
    index = 0;  
    for (row = 1; row < MAX_Alice; row++) {
        while (index < MAX_Alice-row) {
            fp2copy(R->X, pts[npts]->X);
            fp2copy(R->Z, pts[npts]->Z);
            pts_index[npts] = index;
            npts += 1;
            m = strat_Alice_Huff[MAX_Alice-index-row];
            xTPLe_Huff(R, R, A24minus, A24plus, (int)m);
            index += m;
        } 
        get_3_isog_Huff(R, A24minus, A24plus, coeff);

        for (i = 0; i < npts; i++) {
            eval_3_isog_Huff(pts[i], coeff);
        }     
        eval_3_isog_Huff(phiP, coeff);
        eval_3_isog_Huff(phiQ, coeff);
        eval_3_isog_Huff(phiR, coeff);

        fp2copy(pts[npts-1]->X, R->X); 
        fp2copy(pts[npts-1]->Z, R->Z);
        index = pts_index[npts-1];
        npts -= 1;
    }
    
    get_3_isog_Huff(R, A24minus, A24plus, coeff);
    eval_3_isog_Huff(phiP, coeff);
    eval_3_isog_Huff(phiQ, coeff);
    eval_3_isog_Huff(phiR, coeff);

    inv_3_way(phiP->Z, phiQ->Z, phiR->Z);
    fp2mul_mont(phiP->X, phiP->Z, phiP->X);
    fp2mul_mont(phiQ->X, phiQ->Z, phiQ->X);
    fp2mul_mont(phiR->X, phiR->Z, phiR->X);


    // Format public key
    fp2_encode(phiP->X, PublicKeyA);
    fp2_encode(phiQ->X, PublicKeyA + FP2_ENCODED_BYTES);
    fp2_encode(phiR->X, PublicKeyA + 2*FP2_ENCODED_BYTES);

    return 0;
}

// Bob using 5-isogney
int EphemeralKeyGeneration_B_Huff(const unsigned char* PrivateKeyB, unsigned char* PublicKeyB)
{ // Bob's ephemeral public key generation
  // Input:  a private key PrivateKeyB in the range [0, 2^Floor(Log(2,oB)) - 1]. 
  // Output: the public key PublicKeyB consisting of 3 elements in GF(p^2) which are encoded by removing leading 0 bytes.
    point_proj_t R, R2, phiP = {0}, phiQ = {0}, phiR = {0}, pts[MAX_INT_POINTS_BOB];
    f2elm_t XPB, XQB, XRB, coeff[4], CmDsq = {0}, CD4 = {0}, A = {0}, C={0}, D={0};
    unsigned int i, row, m, index = 0, pts_index[MAX_INT_POINTS_BOB], npts = 0, ii = 0;
    digit_t SecretKeyB[NWORDS_ORDER] = {0};

    f2elm_t t0, t1, t2, atmp, ctmp;
    // Initialize basis points
    init_basis((digit_t*)B_gen_Huff, XPB, XQB, XRB);
    init_basis((digit_t*)A_gen_Huff, phiP->X, phiQ->X, phiR->X);
    fpcopy((digit_t*)&Montgomery_one, (phiP->Z)[0]);
    fpcopy((digit_t*)&Montgomery_one, (phiQ->Z)[0]);
    fpcopy((digit_t*)&Montgomery_one, (phiR->Z)[0]);
    fpcopy((digit_t*)&Montgomery_one, D[0]); // D=1
    fpcopy((digit_t*)&Huff_C, C[0]); //C=c

    // Initialize constants: A=c+1/c-2 where c=3+sqrt{8} in GF(p^2) A24minus = (C-D)^2, A24plus = (C+D)^2
    // c+1/c-2 = 4

    fpcopy((digit_t*)&Montgomery_one, A[0]); // A =1
    fp2add(A, A, A); // A =2
    fp2add(A, A, A); // A =4
    fpcopy((digit_t*)&Montgomery_one, CmDsq[0]); // A =1
    fpcopy((digit_t*)&Montgomery_one, CD4[0]); // A =1

    // Retrieve kernel point
    decode_to_digits(PrivateKeyB, SecretKeyB, SECRETKEY_B_BYTES, NWORDS_ORDER);
    LADDER3PT_Huff(XPB, XQB, XRB, SecretKeyB, BOB, R, A);       


    // Traverse tree
    index = 0;        
    for (row = 1; row < MAX_Bob; row++) {
        while (index < MAX_Bob-row) {
            fp2copy(R->X, pts[npts]->X);
            fp2copy(R->Z, pts[npts]->Z);
            pts_index[npts] = index;
            npts += 1;
            m = strat_Bob_Huff[MAX_Bob-index-row];
            x5Pe_Huff(R, R, CmDsq, CD4, (int)(m));
            index += m;
        }
        // kernel points R, [2]R  

        xDBL_Huff(R, R2, CmDsq, CD4);


        get_5_isog_huff(R, R2, C, D, CmDsq, CD4, coeff);
        for (i = 0; i < npts; i++) {
            eval_5_isog_Huff(pts[i], coeff);
        }
        eval_5_isog_Huff(phiP, coeff);
        eval_5_isog_Huff(phiQ, coeff);
        eval_5_isog_Huff(phiR, coeff);

        fp2copy(pts[npts-1]->X, R->X); 
        fp2copy(pts[npts-1]->Z, R->Z);
        index = pts_index[npts-1];
        npts -= 1;
    }

    xDBL_Huff(R, R2, CmDsq, CD4);
    fp2add(R->X, R->Z, coeff[0]);                 // coeff[0] = Xi+Zi
    fp2sub(R->X, R->Z, coeff[1]);                 // coeff[1] = Xi-Zi   
    fp2add(R2->X, R2->Z, coeff[2]);                 // coeff[2] = X2i+Z2i
    fp2sub(R2->X, R2->Z, coeff[3]);                 // coeff[3] = X2i-Z2i   

    eval_5_isog_Huff(phiP, coeff);
    eval_5_isog_Huff(phiQ, coeff);
    eval_5_isog_Huff(phiR, coeff);


    // last step doesnt need to compute the coefficeint

    inv_3_way(phiP->Z, phiQ->Z, phiR->Z);
    fp2mul_mont(phiP->X, phiP->Z, phiP->X);
    fp2mul_mont(phiQ->X, phiQ->Z, phiQ->X);
    fp2mul_mont(phiR->X, phiR->Z, phiR->X);

                
    // Format public key                   
    fp2_encode(phiP->X, PublicKeyB);
    fp2_encode(phiQ->X, PublicKeyB + FP2_ENCODED_BYTES);
    fp2_encode(phiR->X, PublicKeyB + 2*FP2_ENCODED_BYTES);

    return 0;
}



int EphemeralSecretAgreement_A_Huff(const unsigned char* PrivateKeyA, const unsigned char* PublicKeyB, unsigned char* SharedSecretA)
{ // Alice's ephemeral shared secret computation
  // It produces a shared secret key SharedSecretA using her secret key PrivateKeyA and Bob's public key PublicKeyB
  // Inputs: Alice's PrivateKeyA is an integer in the range [0, 2^Floor(Log(2,oA)) - 1]. 
  //         Bob's PublicKeyB consists of 3 elements in GF(p^2) encoded by removing leading 0 bytes.
  // Output: a shared secret SharedSecretA that consists of one element in GF(p^2) encoded by removing leading 0 bytes.  
    point_proj_t R, pts[MAX_INT_POINTS_ALICE];
    f2elm_t coeff[3], PKA[3], jinv;
    f2elm_t A24plus = {0}, A24minus = {0}, A = {0}, CD4={0};
    unsigned int i, row, m, index = 0, pts_index[MAX_INT_POINTS_ALICE], npts = 0, ii = 0;
    digit_t SecretKeyA[NWORDS_ORDER] = {0};
      
    // Initialize images of Alice's basis
    fp2_decode(PublicKeyB, PKA[0]);
    fp2_decode(PublicKeyB + FP2_ENCODED_BYTES, PKA[1]);
    fp2_decode(PublicKeyB + 2*FP2_ENCODED_BYTES, PKA[2]);

    // Initialize constants: A24plus = A+2C, A24minus = A-2C, where C=1
    get_A_Huff(PKA[0], PKA[1], PKA[2], A);
    fp2copy(A, A24minus);
    fpcopy((digit_t*)&Montgomery_one, A24plus[0]);
    fp2add(A24plus, A24plus, A24plus); // 2
    fp2add(A24plus, A24plus, A24plus); // 4
    fp2add(A24minus, A24plus, A24plus);


    // Retrieve kernel point
    decode_to_digits(PrivateKeyA, SecretKeyA, SECRETKEY_A_BYTES, NWORDS_ORDER);
    LADDER3PT_Huff(PKA[0], PKA[1], PKA[2], SecretKeyA, ALICE, R, A);


    
    // Traverse tree
    index = 0;  
    for (row = 1; row < MAX_Alice; row++) {
        while (index < MAX_Alice-row) {
            fp2copy(R->X, pts[npts]->X);
            fp2copy(R->Z, pts[npts]->Z);
            pts_index[npts] = index;
            npts += 1;
            m = strat_Alice_Huff[MAX_Alice-index-row];
            xTPLe_Huff(R, R, A24minus, A24plus, (int)m);
            index += m;
        }
        get_3_isog_Huff(R, A24minus, A24plus, coeff);

        for (i = 0; i < npts; i++) {
            eval_3_isog_Huff(pts[i], coeff);
        } 

        fp2copy(pts[npts-1]->X, R->X); 
        fp2copy(pts[npts-1]->Z, R->Z);
        index = pts_index[npts-1];
        npts -= 1;
    }
     
    get_3_isog_Huff(R, A24minus, A24plus, coeff);  

    fp2sub(A24plus, A24minus, CD4); //                 
                
    j_inv_Huff(A24plus, A24minus, CD4, jinv);
    fp2_encode(jinv, SharedSecretA);    // Format shared secret

    return 0;
}



int EphemeralSecretAgreement_B_Huff(const unsigned char* PrivateKeyB, const unsigned char* PublicKeyA, unsigned char* SharedSecretB)
{ // Bob's ephemeral shared secret computation
  // It produces a shared secret key SharedSecretB using his secret key PrivateKeyB and Alice's public key PublicKeyA
  // Inputs: Bob's PrivateKeyB is an integer in the range [0, 2^Floor(Log(2,oB)) - 1]. 
  //         Alice's PublicKeyA consists of 3 elements in GF(p^2) encoded by removing leading 0 bytes.
  // Output: a shared secret SharedSecretB that consists of one element in GF(p^2) encoded by removing leading 0 bytes.  
    point_proj_t R, R2, pts[MAX_INT_POINTS_BOB], P2={0}, t0={0};
    f2elm_t coeff[4], PKA[3], jinv;
    f2elm_t CmDsq = {0}, CpDsq={0}, CD4 = {0}, A = {0}, C={0}, D={0};
    unsigned int i, row, m, index = 0, pts_index[MAX_INT_POINTS_BOB], npts = 0, ii = 0;
    digit_t SecretKeyB[NWORDS_ORDER] = {0};
      
    // Initialize images of Bob's basis
    fp2_decode(PublicKeyA, PKA[0]);
    fp2_decode(PublicKeyA + FP2_ENCODED_BYTES, PKA[1]);
    fp2_decode(PublicKeyA + 2*FP2_ENCODED_BYTES, PKA[2]);
    fpcopy((digit_t*)&Montgomery_one, D[0]);
    // Initialize constants: A24plus = A+2C, C24 = 4C, where C=1
    get_A_Huff(PKA[0], PKA[1], PKA[2], A); // A = c+1/c-2
    get_C(A, C);



    fp2copy(A, CmDsq); //CmDsq = c+1/c-2
    fpcopy((digit_t*)&Montgomery_one,CD4[0]);
    fp2add(CD4, CD4, CD4); //2
    fp2add(CD4, CD4, CD4); //4

    // Retrieve kernel point
    decode_to_digits(PrivateKeyB, SecretKeyB, SECRETKEY_B_BYTES, NWORDS_ORDER);
    LADDER3PT_Huff(PKA[0], PKA[1], PKA[2], SecretKeyB, BOB, R, A);    


    // Traverse tree
    index = 0;        
    for (row = 1; row < MAX_Bob; row++) {
        while (index < MAX_Bob-row) {
            fp2copy(R->X, pts[npts]->X);
            fp2copy(R->Z, pts[npts]->Z);
            pts_index[npts] = index;
            npts += 1;
            m = strat_Bob_Huff[MAX_Bob-index-row];
            x5Pe_Huff(R, R, CmDsq, CD4, (int)(m));
            index += m;
        }
        // kernel points R, [2]R      
        xDBL_Huff(R, R2, CmDsq, CD4);      

        get_5_isog_huff(R, R2, C, D, CmDsq, CD4, coeff);
        for (i = 0; i < npts; i++) {
            eval_5_isog_Huff(pts[i], coeff);

        }

        // coeffiecient


        fp2copy(pts[npts-1]->X, R->X); 
        fp2copy(pts[npts-1]->Z, R->Z);
        index = pts_index[npts-1];
        npts -= 1;
    }

    xDBL_Huff(R, R2, CmDsq, CD4);    
    get_5_isog_huff(R, R2, C, D, CmDsq, CD4, coeff);

 
    fp2add(CmDsq, CD4, CpDsq);     

                  
    j_inv_Huff(CpDsq, CmDsq, CD4, jinv);
    fp2_encode(jinv, SharedSecretB);    // Format shared secret

    return 0;
}
