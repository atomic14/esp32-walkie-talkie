/***************************************************************************
 *
 *   File Name:  host.c
 *
 *   Purpose:  Contains functions for file I/O and formatting, no signal
 *      processing.
 *
 *      The functions in this file are listed below.  All are level 2
 *      fuctions, where level 0 is main(), except for fillBitAlloc() which
 *      is level 3.  The two "Interface" routines perform truncation of the
 *      three least significant bits of the 16 bit linear input.  The others
 *      are simply file I/O functions and data reformatters.
 *
 *      fillBitAlloc()
 *      hostEncoderInterface()
 *      readDecfile()
 *      speechDecoderHostInterface()
 *      writeEncfile()
 *
 **************************************************************************/

/*_________________________________________________________________________
 |                                                                         |
 |                            Include Files                                |
 |_________________________________________________________________________|
*/

#include <stdio.h>
#include "typedefs.h"

/***************************************************************************
 *
 *   FUNCTION NAME: fillBitAlloc
 *
 *   PURPOSE:
 *
 *     Arrange speech parameters for encoder output
 *
 *   INPUTS:
 *
 *     The speechcoders codewords:
 *     iR0 - Frame energy
 *     piVqIndeces[0:2] - LPC vector quantizer codewords
 *     iSoftInterp - Soft interpolation bit 1 or 0
 *     iVoicing - voicing mode 0,1,2,3
 *     piLags[0:3] - Frame and delta lag codewords
 *     piCodeWrdsA[0:3] - VSELP codevector 1
 *     piCodeWrdsB[0:3] - VSELP codevector 2 (n/a for voiced modes)
 *     piGsp0s[0:3] - GSP0 codewords
 *     swVadFlag - voice activity detection flag
 *     swSP - Speech flag
 *
 *   OUTPUTS:
 *
 *     pswBAlloc[0:20] - an array into which the coded data is moved
 *
 *   RETURN VALUE:
 *
 *     none
 *
 *   REFERENCES: Sub-clause 2.1 and 4.1.12 of GSM Recomendation 06.20
 *
 **************************************************************************/

void   fillBitAlloc(int iVoicing, int iR0, int *piVqIndeces,
                           int iSoftInterp, int *piLags,
                           int *piCodeWrdsA, int *piCodeWrdsB,
                           int *piGsp0s, Shortword swVadFlag,
                           Shortword swSP, Shortword *pswBAlloc)
{

/*_________________________________________________________________________
 |                                                                         |
 |                            Automatic Variables                          |
 |_________________________________________________________________________|
*/

  int    i;
  Shortword *pswNxt;

/*_________________________________________________________________________
 |                                                                         |
 |                            Executable Code                              |
 |_________________________________________________________________________|
*/

  pswNxt = pswBAlloc;
  *pswNxt++ = iR0;
  for (i = 0; i < 3; i++)
    *pswNxt++ = *piVqIndeces++;
  *pswNxt++ = iSoftInterp;
  *pswNxt++ = iVoicing;

  /* check voicing mode */
  if (iVoicing)
  {
    /* voiced mode */
    for (i = 0; i < N_SUB; i++)
    {
      *pswNxt++ = *piLags++;
      *pswNxt++ = *piCodeWrdsA++;
      *pswNxt++ = *piGsp0s++;
    }
  }
  else
  {                                    /* unvoiced frame */
    for (i = 0; i < N_SUB; i++)
    {
      *pswNxt++ = *piCodeWrdsA++;
      *pswNxt++ = *piCodeWrdsB++;
      *pswNxt++ = *piGsp0s++;
    }
  }
  *pswNxt++ = swVadFlag;
  *pswNxt++ = swSP;
}

/***************************************************************************
 *
 *   FUNCTION NAME: hostEncoderInterface
 *
 *   PURPOSE:
 *
 *     Read in speech data from a file.  Zero the least significant 3 bits.
 *
 *
 *   INPUTS:
 *
 *     pfileInSpeech
 *                     FILE pointer to the binary input file
 *
 *     iNumToRead
 *                     Number of samples to read from the file, typically
 *                     160 (20 ms of speech).
 *
 *
 *   OUTPUTS:
 *
 *     pswSamplesRead[]
 *                     The speech samples read in from the file.
 *
 *
 *   RETURN VALUE:
 *
 *     iNumRead
 *                     The number of samples actually read.
 *
 *   IMPLEMENTATION:
 *
 *     The input speech file should be in "native" format. This means that
 *     the file is to be read (by this program) and written (by another
 *     program) as short ints (not chars).
 *
 *     If not enough samples are available in the file, the number actually
 *     read is returned.  If the read fails to fill the requested iNumToRead
 *     samples, then the rest are zeroed.
 *
 *     In all cases the least significant 3 bits of all speech samples are
 *     zeroed.
 *
 *   REFERENCES: Sub-clause 4.1 of GSM Recomendation 06.20
 *
 *   KEYWORDS: read, read speech, get speech data
 *
 **************************************************************************/

int    hostEncoderInterface(FILE *pfileInSpeech, int iNumToRead,
                                   Shortword pswSamplesRead[])
{

/*_________________________________________________________________________
 |                                                                         |
 |                            Automatic Variables                          |
 |_________________________________________________________________________|
*/
  int    iNumRead,
         i;

/*_________________________________________________________________________
 |                                                                         |
 |                              Executable Code                            |
 |_________________________________________________________________________|
*/

  iNumRead = fread((char *) pswSamplesRead, sizeof (Shortword),
                   iNumToRead, pfileInSpeech);

  /* Delete the 3 LSB's - 13 bit speech input */
  /*------------------------------------------*/

  for (i = 0; i < iNumRead; i++)
  {
    pswSamplesRead[i] &= 0xfff8;
  }


  /* Fill out the iNumToRead buffer with zeroes */
  /*--------------------------------------------*/

  if (iNumRead != iNumToRead)
  {
    for (i = iNumRead; i < iNumToRead; i++)
    {
      pswSamplesRead[i] = 0;
    }
  }
  return (iNumRead);
}

/***************************************************************************
 *
 *   FUNCTION NAME:  readDecfile
 *
 *   PURPOSE:
 *      Reads decoder parameter input file
 *
 *   INPUT:
 *      infile           decoder parameter input file.
 *
 *   OUTPUT:
 *      pswSpeechPara    array of received 16-bit values
 *
 *   RETURN:
 *      0                successfully read a complete frame of data
 *
 *   REFERENCES: Sub-clause 4.2 of GSM Recomendation 06.20
 *
 *   KEYWORDS: pswSpeechPara
 *
 **************************************************************************/

int    readDecfile(FILE *infile, Shortword pswSpeechPara[])
{
  if (fread((char *) pswSpeechPara, sizeof (Shortword), 22, infile) == 0)
    return (1);
  else
    return (0);
}

/***************************************************************************
 *
 *   FUNCTION NAME: speechDecoderHostInterface
 *
 *   PURPOSE:
 *     The purpose of this function is to truncate the linear pcm and write
 *     it into the output file
 *
 *   INPUTS:
 *
 *     F_LEN
 *
 *                     160 = linear pcm output frame size
 *
 *     pswDecodedSpeechFrame[0:F_LEN-1]
 *
 *                     16 bit linear pcm
 *
 *   OUTPUTS:
 *
 *     fpfileSpeechOut
 *
 *                     13 bit linear pcm stored to file given by this pointer
 *
 *   RETURN VALUE:
 *
 *     none
 *
 *   IMPLEMENTATION:
 *
 *   REFERENCES: Sub-clause 4.2 of GSM Recomendation 06.20
 *
 *   KEYWORDS: synthesis, speechdecoder, decoding, truncation
 *
 **************************************************************************/

void   speechDecoderHostInterface(Shortword pswDecodedSpeechFrame[],
                                         FILE *fpfileSpeechOut)
{

/*_________________________________________________________________________
 |                                                                         |
 |                              Local Constants                            |
 |_________________________________________________________________________|
*/

#define  PCM_MASK     0xfff8           /* 16 to 13 bit linear PCM mask */

/*_________________________________________________________________________
 |                                                                         |
 |                            Automatic Variables                          |
 |_________________________________________________________________________|
*/

  short int i;

/*_________________________________________________________________________
 |                                                                         |
 |                              Executable Code                            |
 |_________________________________________________________________________|
*/

  /* truncate the 16 bit linear pcm to 13 bits */
  /* ----------------------------------------- */

  for (i = 0; i < F_LEN; i++)
  {
    pswDecodedSpeechFrame[i] = pswDecodedSpeechFrame[i] & PCM_MASK;
  }

  /* F_LEN samples of linear pcm to output file */
  /* ------------------------------------------ */

  fwrite((char *) pswDecodedSpeechFrame, sizeof (Shortword),
         F_LEN, fpfileSpeechOut);
}

/***************************************************************************
 *
 *   FUNCTION NAME:  writeEncfile
 *
 *   PURPOSE:
 *      Writes encoded parameters to ouput file
 *
 *   INPUT:
 *      pswSpeechPara        array of encoded parameter words.
 *
 *   OUTPUT:
 *      fpfileEnc        16-bit encoded values.
 *
 *   RETURN:
 *      i                number of bytes written
 *
 *   REFERENCES: Sub-clause 4.1 of GSM Recomendation 06.20
 *
 *   KEYWORDS: pswSpeechPara, fpfileEnc
 *
 **************************************************************************
*/

int    writeEncfile(Shortword pswSpeechPara[], FILE *fpfileEnc)
{
  int    i;

  i = fwrite((char *) pswSpeechPara, sizeof (Shortword), 20, fpfileEnc);

  return (i);
}
