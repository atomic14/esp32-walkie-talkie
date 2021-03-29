#ifndef __GSM_HR
#define __GSM_HR

#include <stdio.h>

/*_________________________________________________________________________
 |                                                                         |
 |                           Function Prototypes                           |
 |_________________________________________________________________________|
*/

int    encode(FILE *pfileSpeechIn, FILE *pfileEnc);
int    decode(FILE *pfileDec, FILE *pfileSpeechOut);

#endif
