#ifndef __HOMING
#define __HOMING

#include "typedefs.h"

#define EHF_MASK 0x0008                /* Encoder Homing Frame pattern */

/*_________________________________________________________________________
 |                                                                         |
 |                           Function Prototypes                           |
 |_________________________________________________________________________|
*/

int    decoderHomingFrameTest(Shortword pswSpeechPara[], int iLastPara);

void   decoderReset(void);

int    encoderHomingFrameTest(Shortword pswSpeech[]);

void   encoderReset(void);

void   resetDec(void);

void   resetEnc(void);

void   dtxResetTx(void);

void   dtxResetRx(void);

#endif
