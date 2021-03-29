#ifndef __SP_SFRM
#define __SP_SFRM

#include "typedefs.h"

/*_________________________________________________________________________
 |                                                                         |
 |                           Function Prototypes                           |
 |_________________________________________________________________________|
*/

Shortword g_corr2(Shortword *pswIn, Shortword *pswIn2,
                         Longword *pL_out);

int    closedLoopLagSearch(Shortword pswLagList[],
                                  int iNumLags,
                                  Shortword pswLtpState[],
                                  Shortword pswHCoefs[],
                                  Shortword pswPVect[],
                                  Shortword *pswLag,
                                  Shortword *pswLtpShift);

  void   decorr(int iNumVects,
                       Shortword pswGivenVect[],
                       Shortword pswVects[]);

  Shortword g_quant_vl(Shortword swUVCode,
                              Shortword pswWInput[],
                              Shortword swWIShift,
                              Shortword pswWLTPVec[],
                              Shortword pswWVSVec1[],
                              Shortword pswWVSVec2[],
                              struct NormSw snsRs00,
                              struct NormSw snsRs11,
                              struct NormSw snsRs22);

  void   gainTweak(struct NormSw *psErrorTerm);

  void   hnwFilt(Shortword pswInSample[],
                        Shortword pswOutSample[],
                        Shortword pswState[],
                        Shortword pswInCoef[],
                        int iStateOffset,
                        Shortword swZeroState,
                        int iNumSamples);

  void   sfrmAnalysis(Shortword *pswWSpeech,
                             Shortword swVoicingMode,
                             struct NormSw snsSqrtRs,
                             Shortword *pswHCoefs,
                             Shortword *pswLagList,
                             short siNumLags,
                             Shortword swPitch,
                             Shortword swHNWCoef,
                             short *psiLagCode,
                             short *psiVSCode1,
                             short *psiVSCode2,
                             short *psiGsp0Code,
                             Shortword swSP);

  Shortword v_srch(Shortword pswWInput[],
                          Shortword pswWBasisVecs[],
                          short int siNumBasis);

#endif
