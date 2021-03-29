#ifndef __HOSTGSM
#define __HOSTGSM

#include <stdio.h>
#include "typedefs.h"

/*_________________________________________________________________________
 |                                                                         |
 |                           Function Prototypes                           |
 |_________________________________________________________________________|
*/

void   fillBitAlloc(int iVoicing, int iR0, int *piVqIndeces,
                           int iSoftInterp, int *piLags,
                           int *piCodeWrdsA, int *piCodeWrdsB,
                           int *piGsp0s, Shortword swVadFlag,
                           Shortword swSP, Shortword *pswBAlloc);

int    hostEncoderInterface(FILE *pfileInSpeech, int iNumToRead,
                                   Shortword pswSamplesRead[]);

  int    readDecfile(FILE *infile, Shortword pswSpeechPara[]);

  void   speechDecoderHostInterface(Shortword pswDecodedSpeechFrame[],
                                           FILE *fpfileSpeechOut);

  int    writeEncfile(Shortword pswOutBit[], FILE *fpfileEnc);

#endif
