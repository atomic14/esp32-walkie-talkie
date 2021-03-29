/**************************************************************************
 *
 *   File Name:  gsm_hr.c
 *
 *   Purpose:
 *
 *     This file contains the main routine for the GSM Half Rate speech
 *     codec.  The code for the speech coder, including the VAD, DTX,
 *     Comfort Noise, bad frame handling (error concealment), and codec
 *     homing functions is in the files listed below:
 *
 *     gsm_hr.c      globdefs.c      homing.c      host.c
 *     mathdp31.c    mathhalf.c      sp_dec.c      sp_enc.c
 *     sp_frm.c      sp_rom.c        sp_sfrm.c     vad.c
 *     dtx.c         err_conc.c
 *
 *     This particular file only has the level 0, main(), and all level
 *     1, main()'s daughter functions, in it.
 *
 *     Details on how to run gsm_hr are in the readme.txt file.
 *
 *     Below is a listing of all the functions appearing in the file.
 *     The ordering is hierarchical.
 *
 *     main() "gsm_hr", main program.
 *       encode() - encodes a speech file, outputs an encoded parameter file
 *       decode() - decodes a parameter file, outputs a decoded speech file
 *
 **************************************************************************/

/*_________________________________________________________________________
 |                                                                         |
 |                              Include Files                              |
 |_________________________________________________________________________|
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#include "gsm_hr.h"
#include "homing.h"
#include "host.h"
#include "sp_dec.h"
#include "sp_enc.h"
#include "typedefs.h"

/****************************************************************************
 *
 *   PROGRAM NAME:  gsm_hr
 *
 *   PURPOSE:
 *      gsm_hr - main program for the GSM Half-Rate speech coder.
 *      Allows running one of the following operational modes:
 *         0. Encoder only with speech input / encoded parameter data output.
 *         1. Decoder only with speech parameter data input / speech output.
 *
 *   INPUT:
 *      inputFileName     0. speech file in 8 kHz, pcm format.
 *                        1. speech parameter file in decoder input format.
 *      <number of frames> - number of frames to be processed (optional).
 *      <dtx|nodtx>        - switch to enable/disable the DTX functions.
 *
 *   OUTPUT:
 *      outputFileName    0. encoded paramter output file.
 *                        1. synthesized speech file in 8 kHz, pcm format.
 *
 *   RETURN:
 *      none
 *
 *   REFERENCES: Sub-clause 4.0 of GSM Recomendation 06.20
 *
 *   KEYWORDS: main, gsm_hr, speech codec, top level
 *
 ***************************************************************************/

int    main(int argc, char *argv[])
{

/*_________________________________________________________________________
 |                                                                         |
 |                            Local Constants                              |
 |_________________________________________________________________________|
*/

#define  DEFAULT_NUMFRAMES  32766

/*_________________________________________________________________________
 |                                                                         |
 |                            Automatic Variables                          |
 |_________________________________________________________________________|
*/

  int    iDoneFrm,
         iMaxNumFrms,
         option,
         i;

  FILE  *pfileInFile,
        *pfileOutFile;

/*_________________________________________________________________________
 |                                                                         |
 |                              Executable Code                            |
 |_________________________________________________________________________|
*/

  /* check command line arguments */
  /* ---------------------------- */

  iMaxNumFrms = DEFAULT_NUMFRAMES;
  giDTXon = 1;

  /* check command line arguments */
  /* ---------------------------- */

  switch (argc)
  {
    case 4:

      break;

    case 5:
    case 6:

      for (i = 4; i < argc; i++)
      {
        if (!strcmp(argv[i], "dtx"))
          giDTXon = 1;
        else if (!strcmp(argv[i], "nodtx"))
          giDTXon = 0;
        else if ((iMaxNumFrms = atoi(argv[i])) <= 0)
        {
          printf("invalid number of frames or wrong DTX switch, %s", argv[i]);
          exit(1);
        }
      }
      break;

    default:

      printf("\n\nUsage:\n\n");
      printf("gsm_hr option inputFileName outputFileName ");
      printf("<number of frames> <dtx|nodtx>\n");
      printf(" or\n");
      printf("gsm_hr option inputFileName outputFileName ");
      printf("<dtx|nodtx> <number of frames>\n");
      printf("\nOptions: ");
      printf("enc    ");
      puts("inputFileName: speech --> outputFileName: encoder output");
      printf("         ");
      printf("dec    ");
      puts("inputFileName: decoder input --> outputFileName: speech");
      printf("\n");
      exit(1);
  }

  /* check encoding and decoding options */
  /* ----------------------------------- */

  if (!strcmp(argv[1], "enc"))
    option = 0;
  else if (!strcmp(argv[1], "dec"))
    option = 1;
  else
  {
    printf("error in option selection\n");
    printf(" Your entry : %s \n", argv[1]);
    printf("\n\nUsage:\n\n");
    printf("gsm_hr option inputFileName outputFileName ");
    printf("<number of frames> <dtx|nodtx>\n");
    printf(" or\n");
    printf("gsm_hr option inputFileName outputFileName ");
    printf("<dtx|nodtx> <number of frames>\n");
    printf("\nOptions: ");
    printf("enc    ");
    puts("inputFileName: speech --> outputFileName: encoder output");
    printf("         ");
    printf("dec    ");
    puts("inputFileName: decoder input --> outputFileName: speech");
    printf("\n");
    exit(1);
  }


  /* open the input and output files */
  /* ------------------------------- */

#ifdef VAX
  pfileInFile = fopen(argv[2], "rb", "mrs=2", "rfm=fix", "ctx=stm");
  pfileOutFile = fopen(argv[3], "wb", "mrs=2", "rfm=fix", "ctx=stm");
#else
  pfileInFile = fopen(argv[2], "rb");
  pfileOutFile = fopen(argv[3], "wb");
#endif

  if (pfileInFile == NULL)
  {
    printf("error: can not open file %s\n", argv[2]);
    exit(1);
  }
  if (pfileOutFile == NULL)
  {
    printf("error: can not open file %s\n", argv[3]);
    exit(1);
  }

  puts("\n\n");
  puts("   ****************************************");
  puts("   *                                      *");
  puts("   *      GSM Half-Rate Speech Coder      *");
  puts("   *                                      *");
  puts("   *                                      *");
  printf("   *        %s            *\n", VERSION);
  printf("   *        %s            *\n", DATE);
  puts("   *                                      *");
  puts("   ****************************************");
  puts("\n");


  printf("option:       ");

  switch (option)
  {
    case 0:
      puts("enc (speech encoder)");
      break;
    case 1:
      puts("dec (speech decoder)");
      break;
    default:
      puts("invalid option");
      exit(1);
  }

  if (giDTXon)
    printf("DTX mode:     enabled\n");
  else
    printf("DTX mode:     disabled\n");

  printf("input file:   %s\n", argv[2]);
  printf("output file:  %s\n\n", argv[3]);


  switch (option)
  {
    case 0:                            /* encode */

      /* start the encoder, VAD, and transmit DTX in the home state */
      /* ---------------------------------------------------------- */

      resetEnc();

      /* encode:  analyze 8 kHz speech and output encoded parameter file */
      /* --------------------------------------------------------------- */

      for (giFrmCnt = 1, iDoneFrm = 0;
           !iDoneFrm && giFrmCnt <= iMaxNumFrms;
           giFrmCnt++)
      {

#ifndef SILENT
        printf("encoding frame %d \r", giFrmCnt);
#endif

        if (encode(pfileInFile, pfileOutFile))
          iDoneFrm = 1;
      }

      if (iDoneFrm)
        giFrmCnt--;

      printf(" %d speech frames encoded\n", giFrmCnt - 1);
      break;

    case 1:                            /* decode */

      /* start the decoder and receive DTX in the home state */
      /* --------------------------------------------------- */

      resetDec();

      /* decode:  synthesize speech */
      /* -------------------------- */

      for (giFrmCnt = 1, iDoneFrm = 0;
           !iDoneFrm && giFrmCnt <= iMaxNumFrms;
           giFrmCnt++)
      {

#ifndef SILENT
        printf("decoding frame %d \r", giFrmCnt);
#endif

        if (decode(pfileInFile, pfileOutFile))
          iDoneFrm = 1;
      }

      if (iDoneFrm)
        giFrmCnt--;

      printf(" %d speech frames decoded\n", giFrmCnt - 1);
      break;
  }

  fclose(pfileInFile);
  fclose(pfileOutFile);

  return (0);
}

/**************************************************************************
 *
 *   FUNCTION NAME:  decode
 *
 *   PURPOSE:
 *      Reads in one frame of speech parameters and outputs one frame of
 *      synthesized speech.  Resets the decoder to the home state if the
 *      Decoder Homing Frame pattern is detected.
 *
 *   INPUT:
 *      pfileDec        speech parameter input file.
 *
 *   OUTPUT:
 *      pfileSpeechOut  synthesized speech file
 *
 *   RETURN:
 *      0               successfully synthesized a complete frame of speech
 *      1               failed to synthesize a complete frame of speech
 *
 *   REFERENCES: Sub-clause 4.2 of GSM Recomendation 06.20
 *
 *   KEYWORDS:
 *      pfileDec, pfileSpeechOut
 **************************************************************************/

int    decode(FILE *pfileDec, FILE *pfileSpeechOut)
{

/*_________________________________________________________________________
 |                                                                         |
 |                            Local Constants                              |
 |_________________________________________________________________________|
*/

/* These constants define the number of consecutive */
/* parameters that decoderHomingFrameTest checks.   */
/* -------------------------------------------------*/

#define  WHOLE_FRAME  18
#define  TO_FIRST_SUBFRAME  9

/*_________________________________________________________________________
 |                                                                         |
 |                            Static Variables                             |
 |_________________________________________________________________________|
*/
  static Shortword pswSpeechPara[22],
         pswDecodedSpeechFrame[F_LEN];

  static int reset_flag_decoder_old = 1;

/*_________________________________________________________________________
 |                                                                         |
 |                            Automatic Variables                          |
 |_________________________________________________________________________|
*/
  int    i,
         reset_flag_decoder;

/*_________________________________________________________________________
 |                                                                         |
 |                              Executable Code                            |
 |_________________________________________________________________________|
*/

  if (readDecfile(pfileDec, pswSpeechPara))
    return (1);

  if(reset_flag_decoder_old)
    reset_flag_decoder=decoderHomingFrameTest(pswSpeechPara,
					      TO_FIRST_SUBFRAME);
  else
    reset_flag_decoder=0;

  if (reset_flag_decoder && reset_flag_decoder_old)
  {
    /* force the output to be the encoder homing frame pattern */

    for (i = 0; i < F_LEN; i++)
      pswDecodedSpeechFrame[i] = EHF_MASK;
  }
  else
  {
    speechDecoder(pswSpeechPara, pswDecodedSpeechFrame);
  }

  speechDecoderHostInterface(pswDecodedSpeechFrame, pfileSpeechOut);

  if(!reset_flag_decoder_old)
    reset_flag_decoder=decoderHomingFrameTest(pswSpeechPara, WHOLE_FRAME);

  if (reset_flag_decoder)
    resetDec();                        /* bring the decoder and receive DTX
                                        * to the home state */

  reset_flag_decoder_old = reset_flag_decoder;

  return (0);
}

/**************************************************************************
 *
 *   FUNCTION NAME:  encode
 *
 *   PURPOSE:
 *      Reads in one frame of speech samples and outputs one frame of
 *      speech parameters.  Resets the encoder to the home state if the
 *      Encoder Homing Frame pattern is detected.
 *
 *   INPUT:
 *      pfileSpeechIn   speech file
 *
 *   OUTPUT:
 *      pfileEnc        speech, encoded paramater data
 *
 *   RETURN:
 *      0                successfully wrote a complete frame of data
 *      1                failed to write a complete frame of data
 *
 *   REFERENCES: Sub-clause 4.1 of GSM Recomendation 06.20
 *
 *   KEYWORDS:
 *      pfileSpeechIn, pfileEnc
 **************************************************************************/

int    encode(FILE *pfileSpeechIn, FILE *pfileEnc)
{

/*_________________________________________________________________________
 |                                                                         |
 |                            Static Variables                             |
 |_________________________________________________________________________|
*/
  static Shortword pswSpeechPara[20];
  static Shortword pswSpeechBuff[F_LEN];

/*_________________________________________________________________________
 |                                                                         |
 |                            Automatic Variables                          |
 |_________________________________________________________________________|
*/
  int    iNumRead,
         reset_flag;

/*_________________________________________________________________________
 |                                                                         |
 |                              Executable Code                            |
 |_________________________________________________________________________|
*/

  iNumRead = hostEncoderInterface(pfileSpeechIn, F_LEN, &pswSpeechBuff[0]);

  if (iNumRead < F_LEN)
    return (1);

  reset_flag = encoderHomingFrameTest(&pswSpeechBuff[0]);

  speechEncoder(&pswSpeechBuff[0], pswSpeechPara);

  if (writeEncfile(pswSpeechPara, pfileEnc) != 20)
    return (1);

  if (reset_flag)
    resetEnc();                        /* Bring the encoder, VAD, and DTX to
                                        * the home state */

  return (0);
}
