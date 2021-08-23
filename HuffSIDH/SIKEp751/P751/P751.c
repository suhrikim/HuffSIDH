/********************************************************************************************
* Supersingular Isogeny Key Encapsulation Library
*
* Abstract: supersingular isogeny parameters and generation of functions for P751
*********************************************************************************************/  

#include "api.h" 
#include "P751_internal.h"
#include "../internal.h"


// Encoding of field elements, elements over Z_order, elements over GF(p^2) and elliptic curve points:
// --------------------------------------------------------------------------------------------------
// Elements over GF(p) and Z_order are encoded with the least significant octet (and digit) located at the leftmost position (i.e., little endian format). 
// Elements (a+b*i) over GF(p^2), where a and b are defined over GF(p), are encoded as {a, b}, with a in the least significant position.
// Elliptic curve points P = (x,y) are encoded as {x, y}, with x in the least significant position. 
// Internally, the number of digits used to represent all these elements is obtained by approximating the number of bits to the immediately greater multiple of 32.
// For example, a 751-bit field element is represented with Ceil(751 / 64) = 12 64-bit digits or Ceil(751 / 32) = 24 32-bit digits.

//
// Curve isogeny system "SIDHp751". Base curve: Montgomery curve By^2 = Cx^3 + Ax^2 + Cx defined over GF(p751^2), where A=6, B=1, C=1 and p751 = 2^372*3^239-1
//
         
const uint64_t p751[NWORDS64_FIELD]              = { 0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF, 0xEEAFFFFFFFFFFFFF,
                                                     0xE3EC968549F878A8, 0xDA959B1A13F7CC76, 0x084E9867D6EBE876, 0x8562B5045CB25748, 0x0E12909F97BADC66, 0x00006FE5D541F71C };
const uint64_t p751x2[NWORDS64_FIELD]            = { 0xFFFFFFFFFFFFFFFE, 0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF, 0xDD5FFFFFFFFFFFFF, 
                                                     0xC7D92D0A93F0F151, 0xB52B363427EF98ED, 0x109D30CFADD7D0ED, 0x0AC56A08B964AE90, 0x1C25213F2F75B8CD, 0x0000DFCBAA83EE38 }; 
const uint64_t p751x4[NWORDS64_FIELD]            = { 0xFFFFFFFFFFFFFFFC, 0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF, 0xBABFFFFFFFFFFFFF, 
                                                     0x8FB25A1527E1E2A3, 0x6A566C684FDF31DB, 0x213A619F5BAFA1DB, 0x158AD41172C95D20, 0x384A427E5EEB719A, 0x0001BF975507DC70 }; 
const uint64_t p751p1[NWORDS64_FIELD]            = { 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0xEEB0000000000000,
                                                     0xE3EC968549F878A8, 0xDA959B1A13F7CC76, 0x084E9867D6EBE876, 0x8562B5045CB25748, 0x0E12909F97BADC66, 0x00006FE5D541F71C };   
const uint64_t p751x16p[2*NWORDS64_FIELD]        = { 0x0000000000000010, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x2A00000000000000, 
                                                     0x826D2F56C0F0EAE2, 0xAD4C9CBD81067123, 0xF62CF3052282F124, 0x53A95F7469B516FE, 0x3DADEC0D08A4732F, 0x58AD934557C11C7E, 
                                                     0x7F731B89B2DA43F2, 0x51AE9F5F5F6AFF3B, 0xD74319A6C9BCA375, 0x5BAB790796CF84D4, 0xA421554FE2E49CA8, 0x20AD617C8DF437CF, 
                                                     0x3AB06E7A12F5FF7B, 0x70A25E037E40347E, 0x51F1D323FB4C1151, 0xAE0D99AA4835FED9, 0xDF5429960D2536B6, 0x000000030E91D466 };
// Order of Alice's subgroup
const uint64_t Alice_order[NWORDS64_ORDER]       = { 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0010000000000000 }; 
// Order of Bob's subgroup
const uint64_t Bob_order[NWORDS64_ORDER]         = { 0xC968549F878A8EEB, 0x59B1A13F7CC76E3E, 0xE9867D6EBE876DA9, 0x2B5045CB25748084, 0x2909F97BADC66856, 0x06FE5D541F71C0E1 };
// Alice's generator values {XPA0 + XPA1*i, XQA0 + xQA1*i, XRA0 + XRA1*i} in GF(p751^2), expressed in Montgomery representation
const uint64_t A_gen[6 * NWORDS64_FIELD]         = { 0x884F46B74000BAA8, 0xBA52630F939DEC20, 0xC16FB97BA714A04D, 0x082536745B1AB3DB, 0x1117157F446F9E82, 0xD2F27D621A018490,
                                                     0x6B24AB523D544BCD, 0x9307D6AA2EA85C94, 0xE1A096729528F20F, 0x896446F868F3255C, 0x2401D996B1BFF8A5, 0x00000EF8786A5C0A,   // XPA0
                                                     0xAEB78B3B96F59394, 0xAB26681E29C90B74, 0xE520AC30FDC4ACF1, 0x870AAAE3A4B8111B, 0xF875BDB738D64EFF, 0x50109A7ECD7ED6BC,
                                                     0x4CC64848FF0C56FB, 0xE617CB6C519102C9, 0x9C74B3835921E609, 0xC91DDAE4A35A7146, 0x7FC82A155C1B9129, 0x0000214FA6B980B3,   // XPA1
                                                     0x0F93CC38680A8CA9, 0x762E733822E7FED7, 0xE549F005AC0ADB67, 0x94A71FDD2C43A4ED, 0xD48645C2B04721C5, 0x432DA1FE4D4CA4DC,
                                                     0xBC99655FAA7A80E8, 0xB2C6D502BCFD4823, 0xEE92F40CA2EC8BDB, 0x7B074132EFB6D16C, 0x3340B46FA38A7633, 0x0000215749657F6C,   // XQA0
                                                     0xECFF375BF3079F4C, 0xFBFE74B043E80EF3, 0x17376CBE3C5C7AD1, 0xC06327A7E29CDBF2, 0x2111649C438BF3D4, 0xC1F9298261BA2E97,
                                                     0x1F9FECE869CFD1C2, 0x01A39B4FC9346D62, 0x147CD1D3E82A3C9F, 0xDE84E9D249E533EE, 0x1C48A5ADFB7C578D, 0x000061ACA0B82E1D,   // XQA1
                                                     0x1600C525D41059F1, 0xA596899A0A1D83F7, 0x6BFDEED6D2B23F35, 0x5C7E707270C23910, 0x276CA1A4E8369411, 0xB193651A602925A0,
                                                     0x243D239F1CA1F04A, 0x543DC6DA457860AD, 0xCDA590F325181DE9, 0xD3AB7ACFDA80B395, 0x6C97468580FDDF7B, 0x0000352A3E5C4C77,   // XRA0
                                                     0x9B794F9FD1CC3EE8, 0xDB32E40A9B2FD23E, 0x26192A2542E42B67, 0xA18E94FCA045BCE7, 0x96DC1BC38E7CDA2D, 0x9A1D91B752487DE2,
                                                     0xCC63763987436DA3, 0x1316717AACCC551D, 0xC4C368A4632AFE72, 0x4B6EA85C9CCD5710, 0x7A12CAD582C7BC9A, 0x00001C7E240149BF }; // XRA1
// Bob's generator values {XPB0, XQB0, XRB0 + XRB1*i} in GF(p751^2), expressed in Montgomery representation
const uint64_t B_gen[6 * NWORDS64_FIELD]         = { 0x85691AAF4015F88C, 0x7478C5B8C36E9631, 0x7EF2A185DE4DD6E2, 0x943BBEE46BEB9DC7, 0x1A3EC62798792D22, 0x791BC4B084B31D69,
                                                     0x03DBE6522CEA17C4, 0x04749AA65D665D83, 0x3D52B5C45EF450F3, 0x0B4219848E36947D, 0xA4CF7070466BDE27, 0x0000334B1FA6D193,   // XPB0
                                                     0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000,
                                                     0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000,   // XPB1
                                                     0x8E7CB3FA53211340, 0xD67CE54F7A05EEE0, 0xFDDC2C8BCE46FC38, 0x08587FAE3110DF1E, 0xD6B8246FA22B058B, 0x4DAC3ACC905A5DBD,
                                                     0x51D0BF2FADCED3E8, 0xE5A2406DF6484425, 0x907F177584F671B8, 0x4738A2FFCCED051C, 0x2B0067B4177E4853, 0x00002806AC948D3D,   // XQB0
                                                     0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000,
                                                     0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000,   // XQB1
                                                     0xB56457016D1D6D1C, 0x03DECCB38F39C491, 0xDFB910AC8A559452, 0xA9D0F17D1FF24883, 0x8562BBAF515C248C, 0x249B2A6DDB1CB67D,
                                                     0x3131AF96FB46835C, 0xE10258398480C3E1, 0xEAB5E2B872D4FAB1, 0xB71E63875FAEB1DF, 0xF8384D4F13757CF6, 0x0000361EC9B09912,   // XRB0
                                                     0x58C967899ED16EF4, 0x81998376DC622A4B, 0x3D1C1DCFE0B12681, 0x9347DEBB953E1730, 0x9ABB344D3A82C2D7, 0xE4881BD2820552B2,
                                                     0x0037247923D90266, 0x2E3156EDB157E5A5, 0xF86A46A7506823F7, 0x8FE5523A7B7F1CFC, 0xFA3CFFA38372F67B, 0x0000692DCE85FFBD }; // XRB1


// Huff initial curve c
const uint64_t Huff_c[NWORDS64_FIELD]            = { 0x2AFD75A913F3D5E7, 0x2918FBA06F88C9AB, 0xA4AC4DC7CB526F05, 0x2D19E9391A607300, 0x07A79E2B34091B54 ,0x3AD809DCB42F1792, 
                                                     0xD46179328BD6402A, 0x1AFA73541E2C4F3F, 0xF602D73ACE9BDBD8, 0xD77AC58F6BAB7004, 0x4689D97F6793B3B3, 0x00004F26B00E42B7 };  

// Alice's generator values {XPA0 + XPA1*i, XQA0 + xQA1*i, XRA0 + XRA1*i} in GF(p751^2), expressed in Montgomery representation
const uint64_t A_gen_Huff[6 * NWORDS64_FIELD]    = { 0x958e4f66d4e30696, 0xa695e95da63b81d0, 0x846ab67831ebf470, 0x843b2db567d409af, 0x49ce84badb09692e, 0x0e8ba1d891278a4d,
                                                     0x6167cd229663beec, 0x77debb8107317f33, 0x2722972c9a522dd1, 0x7f24f70661d391c7, 0x595ccc2d58a83bef, 0x00005636ee3b71e1,   // XPA0
                                                     0x4fac3aa34730fe8a, 0x7cc67e782c71c806, 0x0877effc615c34f3, 0xc35ef339210bf057, 0x800cc1b9de3194ea, 0x2734c77eca924eaf,
                                                     0xd92ff46450f85363, 0x135fe698faf7e881, 0xdea2df73a1651607, 0x38fee25c42399106, 0x545ea488b8ea82ee, 0x00000e929459eb47,   // XPA1
                                                     0x5ef041e17a450b8c, 0x8b365cbb68639ee6, 0xa6e4af362e25cffc, 0x4eda3ef32bbd13db, 0x74aa780b22161dd8, 0x03d8defaf7e6ae18,
                                                     0x69e893b71cf1e45e, 0x77533d6a6d02fab9, 0x0b9d3b6d65e219de, 0x6bc00825a4c926ac, 0x4a48a7940febeb79, 0x000027869f6d2a3b,   // XQA0
                                                     0xe0a4fd73b64fdb72, 0x06470511e69a3cbe, 0x941f2a5c7ceaf65e, 0xb8ced52908604fd9, 0x64311db03d6ebe40, 0x7ebd0be5bdb78c28,
                                                     0x3f203c03de5f5029, 0x22fb7a36d2149c71, 0xf595187f4f22a1a0, 0x561ffe52a7c4be4d, 0x1f87b19ea8184e57, 0x000030b2fc41faa4,   // XQA1
                                                     0x5e5a6a1b19d3cb27, 0xa48402b5473d67fe, 0x3bc54e179a0f358a, 0x761c1d65571bd436, 0x66681f5f46ceebc9, 0x2f240847d85de10e,
                                                     0xf4fffaac33b8f114, 0x34fdf9ee4bd7c7a1, 0x67cb7704d8bcd0bd, 0x3473ed13aa3171f6, 0xa16851cd5c0f32e6, 0x00001a40078e2029,   // XRA0
                                                     0xe2bdb7f99e549392, 0xff945cdee9768821, 0x747afc549fcc9d00, 0x1bdc07931579eb93, 0x190a4f5dfe840130, 0x09464f5a378331b7,
                                                     0xa0b318073e1dd73c, 0x5750a477fb3dfa88, 0xc98628df003ed02f, 0x54bb1e8bbc987896, 0x898a230be0d9cda8, 0x00002de00ebe03f4 }; // XRA1







// Bob's generator values {XPB0, XQB0, XRB0 + XRB1*i} in GF(p751^2), expressed in Montgomery representation
const uint64_t B_gen_Huff[6 * NWORDS64_FIELD]    = { 0x2b2d3fedbf1b4d99, 0x1f251c2204733675, 0xd844df106ba2187a, 0x4eaef28a1022e99c, 0x73aaafc086112371, 0x065dd481f2daade0,
                                                     0xd073f8f53e8b5af3, 0x77a9fd9c2ebbc91f, 0x3c0f6f288a0e8c58, 0x7e7da6c52ad77b36, 0xb4e5b7e0fef99e3b, 0x0000503c80caeb63,   // XPB0
                                                     0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000,
                                                     0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000,   // XPB1
                                                     0xb865084d0d145c44, 0xf3439c9cf67cdd73, 0xa63b57d4682658c7, 0x326c81f50c5a16a8, 0xcfb629b34a2988f5, 0x9ba96a11ec87eb83,
                                                     0x71f28015d57ff391, 0xcb424447bd75704e, 0x14ade156bb947990, 0xc7c9f77de7e0789b, 0x070a38f70a3f68b8, 0x0000086033ead2dd,   // XQB0
                                                     0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000,
                                                     0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000,   // XQB1
                                                     0xa052828d4cdaa258, 0x07b86a778074f614, 0x514abb8a920d7224, 0x7ed63447c4d89504, 0xf133b05da41433d1, 0xdb3a7b0e836ceae6,
                                                     0xfe35e8eb3689477b, 0xdeba505ca0bdedec, 0x0ef6e111b98ac83a, 0x3608488a33aa31ea, 0x2335bd9ee1fa7851, 0x00000e97149c6337,   // XRB0
                                                     0x467e98aaaf0da207, 0xaa9f9fcb39f10f78, 0x4a073a14848e2a54, 0xa28c9a7ee8fb7d49, 0x5bc91d9683f524cf, 0x3dbf950cba1b7788,
                                                     0x9db4c268b16820a3, 0xd9256fb51d886c83, 0xab76763c6e4bc390, 0xce072d82215dd71f, 0xcc551fb7476543e0, 0x0000139099eecc84 }; // XRB1



// Montgomery constant Montgomery_R2 = (2^768)^2 mod p751
const uint64_t Montgomery_R2[NWORDS64_FIELD]     = { 0x233046449DAD4058, 0xDB010161A696452A, 0x5E36941472E3FD8E, 0xF40BFE2082A2E706, 0x4932CCA8904F8751 ,0x1F735F1F1EE7FC81, 
                                                     0xA24F4D80C1048E18, 0xB56C383CCDB607C5, 0x441DD47B735F9C90, 0x5673ED2C6A6AC82A, 0x06C905261132294B, 0x000041AD830F1F35 };                                                    
// Value one in Montgomery representation 
const uint64_t Montgomery_one[NWORDS64_FIELD]    = { 0x00000000000249ad, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x8310000000000000,
                                                     0x5527b1e4375c6c66, 0x697797bf3f4f24d0, 0xc89db7b2ac5c4e2e, 0x4ca4b439d2076956, 0x10f7926c7512c7e9, 0x00002d5b24bce5e2 };


// Fixed parameters for isogeny tree computation
const unsigned int strat_Alice[MAX_Alice-1] = { 
80, 48, 27, 15, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 7, 4, 2, 1, 1, 2, 1, 
1, 3, 2, 1, 1, 1, 1, 12, 7, 4, 2, 1, 1, 2, 1, 1, 3, 2, 1, 1, 1, 1, 5, 3, 2, 1, 1, 
1, 1, 2, 1, 1, 1, 21, 12, 7, 4, 2, 1, 1, 2, 1, 1, 3, 2, 1, 1, 1, 1, 5, 3, 2, 1, 
1, 1, 1, 2, 1, 1, 1, 9, 5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1, 4, 2, 1, 1, 1, 2, 1, 1, 
33, 20, 12, 7, 4, 2, 1, 1, 2, 1, 1, 3, 2, 1, 1, 1, 1, 5, 3, 2, 1, 1, 1, 1, 2, 1, 
1, 1, 8, 5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1, 4, 2, 1, 1, 2, 1, 1, 16, 8, 4, 2, 1, 1, 
1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1 };

const unsigned int strat_Bob[MAX_Bob-1] = { 
112, 63, 32, 16, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 8, 4, 2, 1, 1, 2, 1, 
1, 4, 2, 1, 1, 2, 1, 1, 16, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 8, 4, 2, 
1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 31, 16, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 
1, 1, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 15, 8, 4, 2, 1, 1, 2, 1, 1, 4, 
2, 1, 1, 2, 1, 1, 7, 4, 2, 1, 1, 2, 1, 1, 3, 2, 1, 1, 1, 1, 49, 31, 16, 8, 4, 2, 
1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 
15, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 7, 4, 2, 1, 1, 2, 1, 1, 3, 2, 1, 
1, 1, 1, 21, 12, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 5, 3, 2, 1, 1, 1, 1, 
2, 1, 1, 1, 9, 5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1, 4, 2, 1, 1, 1, 2, 1, 1 };

// Setting up macro defines and including GF(p), GF(p^2), curve, isogeny and kex functions
#define fpcopy                        fpcopy751
#define fpzero                        fpzero751
#define fpadd                         fpadd751
#define fpsub                         fpsub751
#define fpneg                         fpneg751
#define fpdiv2                        fpdiv2_751
#define fpcorrection                  fpcorrection751
#define fpmul_mont                    fpmul751_mont
#define fpsqr_mont                    fpsqr751_mont
#define fpinv_mont                    fpinv751_mont
#define fpinv_chain_mont              fpinv751_chain_mont
#define fpinv_mont_bingcd             fpinv751_mont_bingcd
#define fp2copy                       fp2copy751
#define fp2zero                       fp2zero751
#define fp2add                        fp2add751
#define fp2sub                        fp2sub751
#define mp_sub_p2                     mp_sub751_p2
#define mp_sub_p4                     mp_sub751_p4
#define sub_p4                        mp_sub_p4
#define fp2neg                        fp2neg751
#define fp2div2                       fp2div2_751
#define fp2correction                 fp2correction751
#define fp2mul_mont                   fp2mul751_mont
#define fp2sqr_mont                   fp2sqr751_mont
#define fp2inv_mont                   fp2inv751_mont
#define fp2inv_mont_bingcd            fp2inv751_mont_bingcd
#define fpequal_non_constant_time     fpequal751_non_constant_time
#define mp_add_asm                    mp_add751_asm
#define mp_subaddx2_asm               mp_subadd751x2_asm
#define mp_dblsubx2_asm               mp_dblsub751x2_asm


#include "../fpx.c"
#include "../ec_isogeny.c"
#include "../sidh.c"
#include "../sike.c"