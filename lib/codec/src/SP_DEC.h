#ifndef __SP_DEC
#define __SP_DEC

#include "typedefs.h"

/*_________________________________________________________________________
 |                                                                         |
 |                            Function Prototypes                          |
 |_________________________________________________________________________|
*/

void   speechDecoder(Shortword pswParameters[],
                            Shortword pswDecodedSpeechFrame[]);

  void   aFlatRcDp(Longword *pL_R, Shortword *pswRc);

  void   b_con(Shortword swCodeWord, short siNumBits,
                      Shortword pswVectOut[]);

  void   fp_ex(Shortword swOrigLagIn, Shortword pswLTPState[]);

  Shortword g_corr1(Shortword *pswIn, Longword *pL_out);

  Shortword g_corr1s(Shortword pswIn[], Shortword swEngyRShft,
                            Longword *pL_out);

  void   getSfrmLpc(short int siSoftInterpolation,
                           Shortword swPrevR0, Shortword swNewR0,
                           Shortword pswPrevFrmKs[],
                           Shortword pswPrevFrmAs[],
                           Shortword pswPrevFrmPFNum[],
                           Shortword pswPrevFrmPFDenom[],
                           Shortword pswNewFrmKs[],
                           Shortword pswNewFrmAs[],
                           Shortword pswNewFrmPFNum[],
                           Shortword pswNewFrmPFDenom[],
                           struct NormSw *psnsSqrtRs,
                           Shortword *ppswSynthAs[],
                           Shortword *ppswPFNumAs[],
                           Shortword *ppswPFDenomAs[]);

  void   get_ipjj(Shortword swLagIn,
                         Shortword *pswIp, Shortword *pswJj);

  short int interpolateCheck(Shortword pswRefKs[],
                                    Shortword pswRefCoefsA[],
                                    Shortword pswOldCoefsA[],
                                    Shortword pswNewCoefsA[],
                                    Shortword swOldPer,
                                    Shortword swNewPer,
                                    Shortword swRq,
                                    struct NormSw *psnsSqrtRsOut,
                                    Shortword pswCoefOutA[]);

  void   lpcFir(Shortword pswInput[], Shortword pswCoef[],
                       Shortword pswState[], Shortword pswFiltOut[]);

  void   lpcIir(Shortword pswInput[], Shortword pswCoef[],
                       Shortword pswState[], Shortword pswFiltOut[]);

  void   lpcIrZsIir(Shortword pswCoef[], Shortword pswFiltOut[]);

  void   lpcZiIir(Shortword pswCoef[], Shortword pswState[],
                         Shortword pswFiltOut[]);

  void   lpcZsFir(Shortword pswInput[], Shortword pswCoef[],
                         Shortword pswFiltOut[]);

  void   lpcZsIir(Shortword pswInput[], Shortword pswCoef[],
                         Shortword pswFiltOut[]);

  void   lpcZsIirP(Shortword pswCommonIO[], Shortword pswCoef[]);

  Shortword r0BasedEnergyShft(Shortword swR0Index);

  short  rcToADp(Shortword swAscale, Shortword pswRc[],
                        Shortword pswA[]);

  void   rcToCorrDpL(Shortword swAshift, Shortword swAscale,
                            Shortword pswRc[], Longword pL_R[]);

  void   res_eng(Shortword pswReflecCoefIn[], Shortword swRq,
                        struct NormSw *psnsSqrtRsOut);

  void   rs_rr(Shortword pswExcitation[], struct NormSw snsSqrtRs,
                      struct NormSw *snsSqrtRsRr);

  void   rs_rrNs(Shortword pswExcitation[], struct NormSw snsSqrtRs,
                        struct NormSw *snsSqrtRsRr);

  Shortword scaleExcite(Shortword pswVect[],
                               Shortword swErrTerm, struct NormSw snsRS,
                               Shortword pswScldVect[]);

  Shortword sqroot(Longword L_SqrtIn);

  void   v_con(Shortword pswBVects[], Shortword pswOutVect[],
                      Shortword pswBitArray[], short int siNumBVctrs);

#endif
