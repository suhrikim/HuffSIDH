/********************************************************************************************
* Supersingular Isogeny Key Encapsulation Library
*
* Abstract: supersingular isogeny parameters and generation of functions for P610
*********************************************************************************************/  

#include "api.h" 
#include "P610_internal.h"
#include "../internal.h"


// Encoding of field elements, elements over Z_order, elements over GF(p^2) and elliptic curve points:
// --------------------------------------------------------------------------------------------------
// Elements over GF(p) and Z_order are encoded with the least significant octet (and digit) located at the leftmost position (i.e., little endian format). 
// Elements (a+b*i) over GF(p^2), where a and b are defined over GF(p), are encoded as {a, b}, with a in the least significant position.
// Elliptic curve points P = (x,y) are encoded as {x, y}, with x in the least significant position. 
// Internally, the number of digits used to represent all these elements is obtained by approximating the number of bits to the immediately greater multiple of 32.
// For example, a 610-bit field element is represented with Ceil(610 / 64) = 10 64-bit digits or Ceil(610 / 32) = 20 32-bit digits.

//
// Curve isogeny system "SIDHp610". Base curve: Montgomery curve By^2 = Cx^3 + Ax^2 + Cx defined over GF(p610^2), where A=6, B=1, C=1 and p610 = 2^67*3^175*5^119-1
//
         
const uint64_t p610[NWORDS64_FIELD]              = { 0xFFFFFFFFFFFFFFFF, 0x98D4AA8E88E4E877, 0x5E44CA3E846082DD, 0xFBB5A9C32A859194, 0xE9FE448E3C61D24A,
                                                     0xE33DD37E7619FF4B, 0xD3E1FB37B4236445, 0xE360D148D97F2482, 0xC9C2C1E5A0BB4E65, 0x1998BB83972C      };
                                        
// 2p
uint64_t p610x2[NWORDS64_FIELD]                  = { 0xFFFFFFFFFFFFFFFE, 0x31A9551D11C9D0EF, 0xBC89947D08C105BB, 0xF76B5386550B2328, 0xD3FC891C78C3A495,
                                                     0xC67BA6FCEC33FE97, 0xA7C3F66F6846C88B, 0xC6C1A291B2FE4905, 0x938583CB41769CCB, 0x333177072E59      }; 

// 4p()
const uint64_t p610x4[NWORDS64_FIELD]            = { 0xfffffffffffffffc, 0x6352aa3a2393a1df, 0x791328fa11820b76, 0xeed6a70caa164651, 0xa7f91238f187492b, 
                                                     0x8cf74df9d867fd2f, 0x4f87ecded08d9117, 0x8d83452365fc920b, 0x270b079682ed3997, 0x6662ee0e5cb3 };




//p+1
const uint64_t p610p1[NWORDS64_FIELD]            = { 0x0000000000000000, 0x98D4AA8E88E4E878, 0x5E44CA3E846082DD, 0xFBB5A9C32A859194, 0xE9FE448E3C61D24A,
                                                     0xE33DD37E7619FF4B, 0xD3E1FB37B4236445, 0xE360D148D97F2482, 0xC9C2C1E5A0BB4E65, 0x1998BB83972C      };   
//16p^2()
const uint64_t p610x16p[2*NWORDS64_FIELD]        = { 0x0000000000000010, 0xe56aae2ee362f100, 0x61bb4516018b284c, 0x5c9da58cc8c4862f, 0xd2a35cde71135fd5, 
                                                     0xaa508978a12980d8, 0xb54a60f7ebb7fb52, 0x4fdad3208db62200, 0xaedf83f9b1224f7e, 0x42679cc19d338454, 
                                                     0xc7974237c2954be5, 0x423e489de8eccf89, 0x7b14bc934a18c5ee, 0x888f9ed67ba463cf, 0xa793415d98fd48b3, 
                                                     0x98a7d732714e0728, 0xc4c57605b7219055, 0x84416d34277d4684, 0x923c7172103fc851, 0x28f2fbee }; 

// 2-torsion point on Montgomery curve By^2 = Cx^3 + Ax^2 + Cx defined over GF(p610^2), where A=6, B=1, C=1
const uint64_t Mont_P2[NWORDS64_FIELD]           = { 0x847804f4d2c80c6c, 0x0b2a9ba1bd7106d2, 0xd1b56cbde45e1da0, 0x0b0387f612d2a63e, 0x5efcfecea48ba44c,
                                                     0xf50a8d46222cc8d7, 0x71768406ff916c2c, 0x70522aa79f0d4f0a, 0x99bd938147393a5b, 0x1113d784336b      };



// Order of Alice's subgroup
const uint64_t Alice_order[NWORDS64_ORDER]       = { 0x74DCFF880845A503, 0x9EE786D7F6818A34, 0x8C77192637BFE85B, 0x19237ADBC4D30B80, 0x3B9F1FB1CC258 }; 
// Order of Bob's subgroup
const uint64_t Bob_order[NWORDS64_ORDER]         = { 0x03DF99092E953E01, 0x2374E42F0F1538FD, 0xC404DC08D3CFF5EC, 0xA6337F19BCCDB0DA, 0x24EE91F2603 };
// Alice's generator values {XPA0 + XPA1*i, XQA0 + xQA1*i, XRA0 + XRA1*i} in GF(p610^2), expressed in Montgomery representation
const uint64_t A_gen[6 * NWORDS64_FIELD]         = { 0x98ce7aa352ab3c36, 0x80da195d419e75c3, 0xec99654b7df77525, 0x0cb7c8647f383294, 0xf6c085c0d0afd19b,
                                                     0x559f1e82c4d56afe, 0x319be0ab282b67c1, 0x7c967e51ad51e086, 0xfa106315cbd78789, 0x5636b9f0388,   // XPA0
                                                     0x7dc9579b47a6cd88, 0x67ceea82f7d0e5ad, 0x61471489073cb181, 0x9ac240b434f1c6a2, 0xcc4054aeb6c24c78,
                                                     0xb1b25c311bf5b2c2, 0xf9ea2f470404b28d, 0x5dd5785e89223e18, 0x729d8ceded88a182, 0x322b1cd08a2,   // XPA1
                                                     0xe4711bbd08ce7f6e, 0x311e3e447fd41385, 0x03549ca800e277e3, 0xcbee283395d5d15c, 0xccb8310812b32df0,
                                                     0x67658edcc02419e6, 0xa3be85a8debc3e69, 0xba49d95e2ac2e0ea, 0xff453e9dcd1223af, 0x7db03f940b5,   // XQA0
                                                     0xed19f6e95a7d8216, 0xec450e388fb9a650, 0xba569dbb6119b7c0, 0x112e584c1ba43126, 0x12c270edd82b41d3,
                                                     0xbf05efaaa44e489f, 0x51ed5500996c815c, 0x3657ecc2c0de6ff8, 0xd025c44484c3b467, 0xf77947f87d4,   // XQA1
                                                     0x8b7cc3b02d2282ce, 0xd5cc8cc23d1884f6, 0x77a526a4bea1a1e8, 0x95e7b76c014e39d9, 0xd83f30f3dea42e1a,
                                                     0xac94c36e9451d299, 0x5f0b2190c4be7dfa, 0xe7d2dcb4c81ea222, 0xec781a4308163162, 0x7eb6230619e,   // XRA0
                                                     0xb35167781197cef2, 0x5fe95878d2edbb99, 0x6f7679d004b2bd87, 0xca4fd3b0b3282952, 0xedfe9e7a191ce2de,
                                                     0x4daa01ab218fce14, 0xc8e5a190fb49e51e, 0xea3e66b39e0cc336, 0x41f0f3a3f642047d, 0xfdebfc7b29e }; // XRA1






// Bob's generator values {XPB0, XQB0, XRB0 + XRB1*i} in GF(p610^2), expressed in Montgomery representation
const uint64_t B_gen[6 * NWORDS64_FIELD]         = { 0xb5a844ec4d06f342, 0x55d60882ae65363b, 0xfb4141c73d41cbdf, 0xe58d9ffa4e2e6a06, 0x82b669d8baf5746a,
                                                     0xca3ad8532e76bd87, 0x86b8a7fb1e752776, 0x9df8ea988aa98b89, 0x38ab2c590ed07e96, 0xcd35d3b051,   // XPB0
                                                     0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 
                                                     0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000,   // XPB1
                                                     0x5474ad18b02281cf, 0x7dd02abbb733c592, 0x3cf75d339b01d689, 0xa51c1ae8b0ed7e5a, 0xb42cee0eac7e479c,
                                                     0x2b5d26cd43585368, 0xce64ba4b82472488, 0x8be4ce5b973a86b8, 0x4ca40e439b2b10d0, 0x110e4630a825,   // XQB0
                                                     0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 
                                                     0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000,   // XQB1
                                                     0xb6f10f12568532fb, 0x9e850447e9bac94a, 0xf627fb69cba7cdb5, 0xe5b48b06b9f70105, 0x19fbd96c4adecb1f,
                                                     0xa15a728cdca26d30, 0x66d88b8d56a9c50a, 0xf62670bf8f83177b, 0x2b4ab0a5ff85e744, 0x18ca359557af,   // XRB0
                                                     0x9bca53ae65c86950, 0x12930679bde0bbf2, 0x85628bdfc2ac2018, 0x53b341cf42a8ab6b, 0x0313d069809b9e9a,
                                                     0xd56835498544fc45, 0xbe6847a8064f1b5b, 0x9cc24925e1ff9fdf, 0x39ce82a16b80087c, 0xf842d2e6bb2 }; // XRB1

// Alice's generator values {XPA0 + XPA1*i, XQA0 + xQA1*i, XRA0 + XRA1*i} in GF(p610^2), expressed in Montgomery representation
const uint64_t A_gen_Huff[6 * NWORDS64_FIELD]    = { 0xbecc7341ae6211fc, 0x749a4afbe27fd449, 0x57970fe0dd7e366b, 0x8049dd4aad227687, 0xad89dba5f9881da0,
                                                     0x0566a9dacf82994b, 0x4d2a55ee00af1b64, 0xd4af72a7c279caf3, 0xdfaec71a6b29c7d1, 0x7953aff119a,   // XPA0
                                                     0xe471c1a849201548, 0x76b6abf1c012c8f5, 0x16f092c04fd9361c, 0x8abe5a9a3eebae3d, 0x101f58c2c87a8810,
                                                     0xe89bf61226de6a2a, 0xcc2ee70395f91695, 0x8341a3a6faf0b90f, 0xb5d54675afae7328, 0xaa89294f301,
                                                     0x68c8b4748d7c46ee, 0xac3fd51572749b10, 0x510107a480f9d0cc, 0xb42e9e2ed4aa0508, 0x15a4332990cb8926,
                                                     0x0831dee1fbfc555e, 0x81bdfe4f1f1adcfa, 0xc6d32ad193ef97bb, 0xa722abd0c23ef026, 0x452fa33b6cc,
                                                     0xa5b52234a1127add, 0xf4d50a047bae85d2, 0xd4c64d14d78efc17, 0xeb735b0bc77de8ca, 0xafc4505c4faf6db0,
                                                     0xeece9b3a233a9b98, 0x15762b4e41c08c04, 0x13fbc90ee2929b4c, 0x3e520f63396e8d18, 0x1aaf14b9eb7,
                                                     0x9d52f367a93f3dcf, 0x1c63c39ad38c3e4c, 0x9754ffaa2d127a5a, 0x67c3d14977da29c0, 0x4086221b38f17b6c,
                                                     0x5b06dbb17a5442fb, 0xdb15e3f58f4e5e57, 0xc7886535e2f7a197, 0x53c5d30560ca4302, 0x8d56f2fc3c1,
                                                     0x8eb2e7000157aa79, 0xaec98faefa1faab4, 0x6d0c9302a0409534, 0x9a6fa5a875a89d0a, 0x9837391ccac42dd9,
                                                     0xc6fe0bf49969f02b, 0xbf9384cb59009af1, 0xd4cf83b76523551c, 0xe599de4df9ce4be8, 0x13ab6635439d}; // XRA1


// Bob's generator values {XPB0, XQB0, XRB0 + XRB1*i} in GF(p610^2), expressed in Montgomery representation
const uint64_t B_gen_Huff[6 * NWORDS64_FIELD]     = {0x189dca498154c976, 0xd505a6618001f5cb, 0xcef61751434d75c8, 0x8ea62c16699acf41, 0xb933080a50e86304,
                                                     0xe31aa4e9037d0281, 0xfd04ce276eb0d428, 0xf630b71ab0d9768a, 0x51876753c3a69f85, 0x7ad11623252, // XPB0
                                                     0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 
                                                     0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000,   // XPB1
                                                     0x02d3e035540f4214, 0xed17d3f13d81745c, 0x693c46ac2813631c, 0x71633bbf0d4c7d87, 0xbd77652c9751140c,
                                                     0x603779739a3dd281, 0x2da4a91f37feaa4c, 0xca75d5726d6a7c63, 0xf9826222c9bb073f, 0x95fb7e0baaa,    // XQB0
                                                     0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 
                                                     0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000,   // XQB1
                                                     0xd11576ad9868ae38, 0x7d69ba0e172363f6, 0x671ecd5de9f8d7c1, 0x7d2e21553a2c5ab4, 0x9a42e1bdb62a4b34,
                                                     0x466c34e0eb247db4, 0x56018026620f100d, 0x1a15718da6b15a54, 0x02cf6d0f8df0f334, 0x13e038cfc2fe,   // XRB0
                                                     0xf472b61d0c9b865c, 0xbdcf57273abc5c75, 0x79874a6f9a34548b, 0xe9e2dcfd68fc7732, 0x5b42722650f64bdc,
                                                     0x196595ad11f949d7, 0xd602bd64028d79b1, 0x0312131a05cb54f5, 0xc8cb72131741dda7, 0x189414e171cb }; // XRB1


// Initial Huff coefficient c=3+sqrt(8) in GF(p610^2)
const uint64_t Huff_C[NWORDS64_FIELD]            = { 0x847804f4d3040e75, 0x08f369e6a773ea9a, 0x910ecce3d003f6cc, 0xfc8f12766a37472b, 0x172d526877e3ab1d,
                                                     0xea9372cff1795fc0, 0x5980075de3291773, 0xa33543986ca257d2, 0x8a8f4b46ca25624e, 0x666ca363f9f      };   


// Montgomery constant Montgomery_R2 = (2^640)^2 mod p621
const uint64_t Montgomery_R2[NWORDS64_FIELD]     = { 0x6B44F89B409A6EA2, 0xFB890764ECC2540F, 0x874DE858CD9B4918, 0x209ACB6289A7FC14, 0x3FDE7FFBA4CAD985,
                                                     0x1512E8E411027FC9, 0x6079ADC31B8A595B, 0x1620B9F8FCF5BB00, 0xCC19FC3551ED0859, 0xA9D97E36E48      };                                                    
// Value one in Montgomery representation 
const uint64_t Montgomery_one[NWORDS64_FIELD]    = { 0x00000000000A0056, 0xFEFD5B2CEE69E7B0, 0x6E72E33A6AEC113F, 0xCF594F3807556978, 0xB7069C65803585B6,
                                                     0x664A56805A4D1890, 0x01E891E79161701B, 0x70A1328F822D7539, 0xA59A958AF118C7FD, 0x138D1A0B5558     };

// Fixed parameters for isogeny tree computation


const unsigned int strat_Alice_Huff[MAX_Alice] = { 
0, 1, 1, 1, 2, 2, 2, 3, 3, 3, 3, 4, 5, 5, 5, 5, 5, 6, 6, 7, 8,
8, 8, 8, 8, 8, 9, 9, 9, 9, 10, 11, 12, 13, 13, 13, 13, 13, 13, 13, 14, 14, 14,
14, 14, 14, 14, 15, 16, 17, 18, 18, 19, 20, 21, 21, 21, 21, 21, 21, 21, 21, 22,
22, 22, 22, 22, 22, 22, 22, 22, 22, 22, 23, 24, 25, 26, 27, 27, 27, 27, 27, 28,
29, 30, 31, 32, 33, 34, 34, 34, 34, 34, 34, 34, 34, 34, 35, 35, 35, 35, 35, 35,
35, 35, 35, 35, 35, 35, 35, 35, 35, 35, 36, 37, 38, 39, 40, 41, 41, 41, 41, 41,
41, 41, 41, 41, 41, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 51, 52, 53, 54,
55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56,
56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56  };


const unsigned int strat_Bob_Huff[MAX_Bob] = { 
0, 1, 1, 2, 2, 2, 3, 3, 4, 4, 4, 5, 5, 5, 6, 7, 7, 7, 7, 8, 8,
9, 9, 9, 9, 10, 11, 12, 12, 12, 12, 12, 13, 13, 13, 14, 15, 16, 16, 16, 16, 16,
16, 16, 17, 18, 19, 20, 20, 20, 20, 20, 20, 21, 21, 22, 22, 22, 22, 22, 23, 24,
25, 26, 27, 28, 28, 28, 28, 28, 28, 28, 28, 28, 28, 28, 29, 30, 31, 32, 33, 33,
33, 33, 33, 33, 33, 34, 34, 34, 35, 36, 37, 38, 38, 38, 38, 38, 38, 38, 38, 38,
38, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 48, 48, 48, 48, 48 };



// Setting up macro defines and including GF(p), GF(p^2), curve, isogeny and kex functions
#define fpcopy                        fpcopy610
#define fpzero                        fpzero610
#define fpadd                         fpadd610
#define fpsub                         fpsub610
#define fpneg                         fpneg610
#define fpdiv2                        fpdiv2_610
#define fpcorrection                  fpcorrection610
#define fpmul_mont                    fpmul610_mont
#define fpsqr_mont                    fpsqr610_mont
#define fpinv_mont                    fpinv610_mont
#define fpsqrt_mont                   fpsqrt610_mont
#define fpinv_chain_mont              fpinv610_chain_mont
#define fpinv_mont_bingcd             fpinv610_mont_bingcd
#define fp2copy                       fp2copy610
#define fp2zero                       fp2zero610
#define fp2add                        fp2add610
#define fp2sub                        fp2sub610
#define mp_sub_p2                     mp_sub610_p2
#define mp_sub_p4                     mp_sub610_p4
#define sub_p4                        mp_sub_p4
#define fp2neg                        fp2neg610
#define fp2div2                       fp2div2_610
#define fp2correction                 fp2correction610
#define fp2mul_mont                   fp2mul610_mont
#define fp2sqr_mont                   fp2sqr610_mont
#define fp2inv_mont                   fp2inv610_mont
#define fp2inv_mont_bingcd            fp2inv610_mont_bingcd
#define fp2sqrt_mont                  fp2sqrt610_mont
#define fpequal_non_constant_time     fpequal610_non_constant_time
#define mp_add_asm                    mp_add610_asm
#define mp_subaddx2_asm               mp_subadd610x2_asm
#define mp_dblsubx2_asm               mp_dblsub610x2_asm


#include "../fpx.c"
#include "../ec_isogeny.c"
#include "../sidh.c"
#include "../sike.c"