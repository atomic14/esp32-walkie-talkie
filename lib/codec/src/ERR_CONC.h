#ifndef __ERR_CONC
#define __ERR_CONC

#include "typedefs.h"

/*_________________________________________________________________________
 |                                                                         |
 |                            Function Prototypes                          |
 |_________________________________________________________________________|
*/

void   para_conceal_speech_decoder(Shortword pswErrorFlag[],
                       Shortword pswSpeechPara[], Shortword *pswMutePermit);

  Shortword level_calc(Shortword swInd, Longword *pl_en);

  void   level_estimator(Shortword update, Shortword *pswLevelMean,
                                Shortword *pswLevelMax,
                                Shortword pswDecodedSpeechFrame[]);

  void   signal_conceal_sub(Shortword pswPPFExcit[],
                     Shortword ppswSynthAs[], Shortword pswSynthFiltState[],
                       Shortword pswLtpStateOut[], Shortword pswPPreState[],
                                Shortword swLevelMean, Shortword swLevelMax,
                            Shortword swErrorFlag1, Shortword swMuteFlagOld,
                            Shortword *pswMuteFlag, Shortword swMutePermit);

#endif
