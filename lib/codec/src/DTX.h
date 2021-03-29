#ifndef __DTX
#define __DTX

#include "typedefs.h"


#define PN_INIT_SEED (Longword)0x1091988L       /* initial seed for Comfort
                                                 * noise pn-generator */

#define CNINTPER    12                 /* inperpolation period of CN
                                        * parameters */

#define SPEECH      1
#define CNIFIRSTSID 2
#define CNICONT     3
#define CNIBFI      4

#define VALIDSID    11
#define INVALIDSID  22
#define GOODSPEECH  33
#define UNUSABLE    44

/*________________________________________________________________________
 |                                                                        |
 |                      Function Prototypes                               |
 |________________________________________________________________________|
*/

void   avgCNHist(Longword pL_R0History[],
                        Longword ppL_CorrHistory[OVERHANG][NP + 1],
                        Longword *pL_AvgdR0,
                        Longword pL_AvgdCorrSeq[]);

  void   avgGsHistQntz(Longword pL_GsHistory[], Longword *pL_GsAvgd);

  Shortword swComfortNoise(Shortword swVadFlag,
                            Longword L_UnqntzdR0, Longword *pL_UnqntzdCorr);

  Shortword getPnBits(int iBits, Longword *L_PnSeed);

  Shortword gsQuant(Longword L_GsIn, Shortword swVoicingMode);

  void   updateCNHist(Longword L_UnqntzdR0,
                             Longword *pL_UnqntzdCorr,
                             Longword pL_R0Hist[],
                             Longword ppL_CorrHist[OVERHANG][NP + 1]);

  void   lpcCorrQntz(Longword pL_CorrelSeq[],
                            Shortword pswFinalRc[],
                            int piVQCodewds[]);

  Longword linInterpSid(Longword L_New, Longword L_Old, Shortword swDtxState);

  Shortword linInterpSidShort(Shortword swNew,
                                     Shortword swOld,
                                     Shortword swDtxState);

  void   rxInterpR0Lpc(Shortword *pswOldKs, Shortword *pswNewKs,
                              Shortword swRxDTXState,
                              Shortword swDecoMode, Shortword swFrameType);

#endif
