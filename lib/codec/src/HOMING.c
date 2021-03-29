/**************************************************************************
 *
 *   File Name:  homing.c
 *
 *   Purpose:
 *      This file contains the following functions:
 *
 *      decoderHomingFrameTest() - checks if a frame of input speech
 *                                 parameters matches the Decoder Homing
 *                                 Frame pattern.
 *
 *      decoderReset() - called by resetDec() to reset all of the state
 *                       variables for the decoder
 *
 *      encoderHomingFrameTest() - checks if a frame of input samples
 *                                 matches the Encoder Homing Frame pattern.
 *
 *      encoderReset() - called by resetEnc() to reset all the state
 *                       variables for the encoder.
 *
 *      resetDec() - calls functions to reset the state variables for the
 *                   decoder, and for the receive DTX and Comfort Noise.
 *
 *      resetEnc() - calls functions to reset the state variables for the
 *                   encoder and VAD, and for the transmit DTX and
 *                   Comfort Noise.
 *
 *      dtxResetTx() - called by resetEnc() to reset all of the transmit
 *                     DTX and Comfort Noise state variables
 *
 *      dtxResetRx() - called by resetDec() to reset all of the receive
 *                     DTX and Comfort Noise state variables
 *
 **************************************************************************/

/*_________________________________________________________________________
 |                                                                         |
 |                              Include Files                              |
 |_________________________________________________________________________|
*/

#include "typedefs.h"
#include "vad.h"
#include "dtx.h"
#include "homing.h"

/*_________________________________________________________________________
 |                                                                         |
 |                              Local Defines                              |
 |_________________________________________________________________________|
*/

#define EHF_MASK 0x0008                /* Encoder Homing Frame pattern */
#define LMAX 142                       /* largest lag (integer sense) */
#define CG_INT_MACS 6                  /* Number of multiply-accumulates in
                                        * one interpolation           */
#define NUM_CLOSED 3                   /* maximum number of lags searched */
#define LPCSTARTINDEX 25               /* where the LPC analysis window
                                        * starts                        */
#define INBUFFSZ LPCSTARTINDEX + A_LEN /* input buffer size */


#define LTP_LEN         147            /* maximum ltp lag */
#define LSMAX           (LMAX + CG_INT_MACS/2)
#define HNW_BUFF_LEN    LSMAX


/***************************************************************************
 *
 *   FUNCTION NAME:  decoderHomingFrameTest
 *
 *   PURPOSE:
 *      Checks if a frame of input speech parameters matches the Decoder
 *      Homing Frame pattern, which is:
 *
 *      parameter  decimal value  hexidecimal value
 *      ---------  -------------  -----------------
 *      R0         0              0x0000
 *      LPC1       881            0x0371
 *      LPC2       350            0x015E
 *      LPC3       195            0x00c3
 *      INT_LPC    1              0x0001
 *      MODE       0              0x0000
 *      CODE1_1    71             0x0047
 *      CODE2_1    74             0x004a
 *      GSP0_1     0              0x0000
 *      CODE1_2    9              0x0009
 *      CODE2_2    38             0x0026
 *      GSP0_2     7              0x0007
 *      CODE1_3    0              0x0000
 *      CODE2_3    0              0x0000
 *      GSP0_3     0              0x0000
 *      CODE1_4    0              0x0000
 *      CODE2_4    0              0x0000
 *      GSP0_4     0              0x0000
 *
 *   INPUT:
 *      pswSpeechPara[] - one frame of speech parameters
 *                        in decoder input format
 *
 *      iLastPara - the number of consecutive parameters in 
 *                  pswSpeechPara[] to match.
 *
 *   OUTPUT:
 *      None
 *
 *   RETURN:
 *      0       input frame does not match the Decoder Homing Frame pattern.
 *      1       input frame matches the Decoder Homing Frame pattern.
 *
 *   REFERENCES: Sub-clause 10 of GSM Recomendation 06.02
 *
 *   KEYWORDS:
 *      pswSpeechPara
 **************************************************************************/

int    decoderHomingFrameTest(Shortword pswSpeechPara[], int iLastPara)
{
  /* the n[] array contains the number of bits in each speech parameter */
  static int n[] = {5, 11, 9, 8, 1, 2, 7, 7, 5, 7, 7, 5, 7, 7, 5, 7, 7, 5};

  static Shortword dhf_mask[] =
  {
    0x0000,                            /* R0      */
    0x0371,                            /* LPC1    */
    0x015E,                            /* LPC2    */
    0x00c3,                            /* LPC3    */
    0x0001,                            /* INT_LPC */
    0x0000,                            /* MODE    */
    0x0047,                            /* CODE1_1 */
    0x004a,                            /* CODE2_1 */
    0x0000,                            /* GSP0_1  */
    0x0009,                            /* CODE1_2 */
    0x0026,                            /* CODE2_2 */
    0x0007,                            /* GSP0_2  */
    0x0000,                            /* CODE1_3 */
    0x0000,                            /* CODE2_3 */
    0x0000,                            /* GSP0_3  */
    0x0000,                            /* CODE1_4 */
    0x0000,                            /* CODE2_4 */
    0x0000                             /* GSP0_4  */
  };

  int    i;
  int    j;

  for (i = 0; i < iLastPara; i++)
  {
    j = ((pswSpeechPara[i] & ~(~0 << n[i])) ^ dhf_mask[i]);
    if (j)
      break;
  }

  return !j;
}


/***************************************************************************
 *
 *   FUNCTION NAME:  decoderReset
 *
 *   PURPOSE:
 *      resets all of the state variables for the decoder
 *
 *   INPUT:
 *      None
 *
 *   OUTPUT:
 *      None
 *
 *   RETURN:
 *      None
 *
 *   REFERENCES: Sub-clause 10 of GSM Recomendation 06.02
 *
 *   KEYWORDS:
 **************************************************************************/

void   decoderReset(void)
{
/*_________________________________________________________________________
 |                                                                         |
 |   External declarations for decoder variables which need to be reset    |
 |_________________________________________________________________________|
*/

  /* variables defined in sp_dec.c */
  /* ----------------------------- */

  extern Shortword gswPostFiltAgcGain,
         gpswPostFiltStateNum[NP],
         gpswPostFiltStateDenom[NP],
         swPostEmphasisState,
         pswSynthFiltState[NP],
         pswOldFrmKsDec[NP],
         pswOldFrmAsDec[NP],
         pswOldFrmPFNum[NP],
         pswOldFrmPFDenom[NP],
         swOldR0Dec,
         pswLtpStateBaseDec[LTP_LEN + S_LEN],
         pswPPreState[LTP_LEN + S_LEN];

  extern Shortword swMuteFlagOld;      /* error concealment */


  /* variables defined in err_conc.c *//* error concealment */
  /* ------------------------------- *//* error concealment */

  extern Longword plSubfrEnergyMem[4]; /* error concealment */
  extern Shortword swLevelMem[4],
         lastR0,                       /* error concealment */
         pswLastGood[18],              /* error concealment */
         swState,
         swLastFlag;                   /* error concealment */

/*_________________________________________________________________________
 |                                                                         |
 |   Automatic Variables                                                   |
 |_________________________________________________________________________|
*/

  int    i;

/*_________________________________________________________________________
 |                                                                         |
 |   Executable code                                                       |
 |_________________________________________________________________________|
*/

  /* reset all the decoder state variables */
  /* ------------------------------------- */

  swOldR0Dec = 0;

  gswPostFiltAgcGain = 0;

  swPostEmphasisState = 0;

  for (i = 0; i < NP; i++)
  {
    gpswPostFiltStateNum[i] = 0;
    gpswPostFiltStateDenom[i] = 0;
    pswSynthFiltState[i] = 0;
    pswOldFrmKsDec[i] = 0;
    pswOldFrmAsDec[i] = 0;
    pswOldFrmPFNum[i] = 0;
    pswOldFrmPFDenom[i] = 0;
  }

  for (i = 0; i < (LTP_LEN + S_LEN); i++)
  {
    pswLtpStateBaseDec[i] = 0;
    pswPPreState[i] = 0;
  }


  /* reset all the error concealment state variables */
  /* ----------------------------------------------- */

  swMuteFlagOld = 0;

  lastR0 = 0;
  swState = 7;
  swLastFlag = 0;
  for (i = 0; i < 3; i++)
  {
    plSubfrEnergyMem[i] = 80;
    swLevelMem[i] = -72;
  }
  for (i = 0; i < 18; i++)
  {
    pswLastGood[i] = 0;
  }


}

/***************************************************************************
 *
 *   FUNCTION NAME:  encoderHomingFrameTest
 *
 *   PURPOSE:
 *      Checks if a frame of input samples matches the Encoder Homing Frame
 *      pattern, which is 0x0008 for all 160 samples in the frame.
 *
 *   INPUT:
 *      pswSpeech[]    one frame of speech samples
 *
 *   OUTPUT:
 *      None
 *
 *   RETURN:
 *      0       input frame does not match the Encoder Homing Frame pattern.
 *      1       input frame matches the Encoder Homing Frame pattern.
 *
 *   REFERENCES: Sub-clause 10 of GSM Recomendation 06.02
 *
 *   KEYWORDS:
 *      pswSpeech
 **************************************************************************/

int    encoderHomingFrameTest(Shortword pswSpeech[])
{
  int    i;
  Shortword j;

  for (i = 0; i < F_LEN; i++)
  {
    j = pswSpeech[i] ^ EHF_MASK;
    if (j)
      break;
  }

  return !j;
}

/***************************************************************************
 *
 *   FUNCTION NAME:  encoderReset
 *
 *   PURPOSE:
 *      resets all of the state variables for the encoder
 *
 *   INPUT:
 *      None
 *
 *   OUTPUT:
 *      None
 *
 *   RETURN:
 *      None
 *
 *   REFERENCES: Sub-clause 10 of GSM Recomendation 06.02
 *
 *   KEYWORDS:
 **************************************************************************/

void   encoderReset(void)
{

/*_________________________________________________________________________
 |                                                                         |
 |   External declarations for encoder variables which need to be reset    |
 |_________________________________________________________________________|
*/

  /* variables defined in sp_enc.c */
  /* ----------------------------- */

  extern Shortword swOldR0;
  extern Shortword swOldR0Index;

  extern struct NormSw psnsWSfrmEngSpace[];

  extern Shortword pswHPFXState[];
  extern Shortword pswHPFYState[];
  extern Shortword pswOldFrmKs[];
  extern Shortword pswOldFrmAs[];
  extern Shortword pswOldFrmSNWCoefs[];
  extern Shortword pswWgtSpeechSpace[];

  extern Shortword pswSpeech[];        /* input speech */

  extern Shortword swPtch;


  /* variables defined in sp_frm.c */
  /* ----------------------------- */

  extern Shortword pswAnalysisState[NP];

  extern Shortword pswWStateNum[NP],
         pswWStateDenom[NP];


  /* variables defined in sp_sfrm.c */
  /* ------------------------------ */

  extern Shortword pswLtpStateBase[LTP_LEN + S_LEN];
  extern Shortword pswHState[NP];
  extern Shortword pswHNWState[HNW_BUFF_LEN];

/*_________________________________________________________________________
 |                                                                         |
 |   Automatic Variables                                                   |
 |_________________________________________________________________________|
*/

  int    i;

/*_________________________________________________________________________
 |                                                                         |
 |   Executable code                                                       |
 |_________________________________________________________________________|
*/

  /* reset all the encoder state variables */
  /* ------------------------------------- */

  swOldR0Index = 0;
  swOldR0 = 0;

  for (i = 0; i < 2 * N_SUB; i++)
  {
    psnsWSfrmEngSpace[i].man = 0;
    psnsWSfrmEngSpace[i].sh = 0;
  }

  for (i = 0; i < 4; i++)
    pswHPFXState[i] = 0;

  for (i = 0; i < 8; i++)
    pswHPFYState[i] = 0;

  for (i = 0; i < NP; i++)
  {
    pswOldFrmKs[i] = 0;
    pswOldFrmAs[i] = 0;
    pswOldFrmSNWCoefs[i] = 0;
    pswAnalysisState[i] = 0;
    pswWStateNum[i] = 0;
    pswWStateDenom[i] = 0;
    pswHState[i] = 0;
  }

  for (i = 0; i < (F_LEN + LMAX + CG_INT_MACS / 2); i++)
    pswWgtSpeechSpace[i] = 0;

  for (i = 0; i < INBUFFSZ; i++)
    pswSpeech[i] = 0;

  for (i = 0; i < (LTP_LEN + S_LEN); i++)
    pswLtpStateBase[i] = 0;

  for (i = 0; i < HNW_BUFF_LEN; i++)
    pswHNWState[i] = 0;

  swPtch = 1;
}

/***************************************************************************
 *
 *   FUNCTION NAME:  resetDec
 *
 *   PURPOSE:
 *      resets all of the state variables for the decoder, and for the
 *      receive DTX and Comfort Noise.
 *
 *   INPUT:
 *      None
 *
 *   OUTPUT:
 *      None
 *
 *   RETURN:
 *      None
 *
 *   REFERENCES: Sub-clause 10 of GSM Recomendation 06.02
 *
 *   KEYWORDS:
 **************************************************************************/

void   resetDec(void)
{
  decoderReset();                      /* reset all the state variables in
                                        * the speech decoder */
  dtxResetRx();                        /* reset all the receive DTX and CN
                                        * state variables */
}

/***************************************************************************
 *
 *   FUNCTION NAME:  resetEnc
 *
 *   PURPOSE:
 *      resets all of the state variables for the encoder and VAD, and for
 *      the transmit DTX and Comfort Noise.
 *
 *   INPUT:
 *      None
 *
 *   OUTPUT:
 *      None
 *
 *   RETURN:
 *      None
 *
 *   REFERENCES: Sub-clause 10 of GSM Recomendation 06.02
 *
 *   KEYWORDS:
 **************************************************************************/

void   resetEnc(void)
{

  encoderReset();                      /* reset all the state variables in
                                        * the speech encoder */
  vad_reset();                         /* reset all the VAD state variables */

  dtxResetTx();                        /* reset all the transmit DTX and CN
                                        * state variables */
}

/***************************************************************************
 *
 *   FUNCTION NAME:  dtxResetTx
 *
 *   PURPOSE:
 *      reset all the transmit DTX and CN state variables
 *
 *   INPUT:
 *      None
 *
 *   OUTPUT:
 *      None
 *
 *   RETURN:
 *      None
 *
 *   REFERENCES: Sub-clause 10 of GSM Recomendation 06.02
 *
 *   KEYWORDS:
 **************************************************************************/

void   dtxResetTx(void)
{

/*_________________________________________________________________________
 |                                                                         |
 |   External declarations for encoder variables which need to be reset    |
 |_________________________________________________________________________|
*/

  /* variables defined in sp_enc.c */
  /* ----------------------------- */

  extern Shortword swTxGsHistPtr;


  /* variables defined in dtx.c */
  /* -------------------------- */

  extern Shortword swVadFrmCnt;        /* Indicates the number of sequential
                                        * frames where VAD == 0 */
  extern short int siUpdPointer;
  extern Shortword swNElapsed;

/*_________________________________________________________________________
 |                                                                         |
 |   Automatic Variables                                                   |
 |_________________________________________________________________________|
*/


/*_________________________________________________________________________
 |                                                                         |
 |   Executable code                                                       |
 |_________________________________________________________________________|
*/

  /* reset all the transmit DTX and CN state variables */
  /* ------------------------------------------------- */

  swTxGsHistPtr = 0;

  swVadFrmCnt = 0;

  siUpdPointer = 0;

  swNElapsed = 50;

}

/***************************************************************************
 *
 *   FUNCTION NAME:  dtxResetRx
 *
 *   PURPOSE:
 *      reset all the receive DTX and CN state variables
 *
 *   INPUT:
 *      None
 *
 *   OUTPUT:
 *      None
 *
 *   RETURN:
 *      None
 *
 *   REFERENCES: Sub-clause 10 of GSM Recomendation 06.02
 *
 *   KEYWORDS:
 **************************************************************************/

void   dtxResetRx(void)
{

/*_________________________________________________________________________
 |                                                                         |
 |   External declarations for encoder variables which need to be reset    |
 |_________________________________________________________________________|
*/

  /* variables defined in sp_dec.c */
  /* ----------------------------- */

  extern Shortword swRxDTXState;
  extern Shortword swDecoMode;
  extern Shortword swDtxMuting;
  extern Shortword swDtxBfiCnt;

  extern Shortword swOldR0IndexDec;

  extern Shortword swRxGsHistPtr;
  extern Longword pL_RxGsHist[(OVERHANG - 1) * N_SUB];


/*_________________________________________________________________________
 |                                                                         |
 |   Automatic Variables                                                   |
 |_________________________________________________________________________|
*/

  int    i;


/*_________________________________________________________________________
 |                                                                         |
 |   Executable code                                                       |
 |_________________________________________________________________________|
*/

  /* reset all the receive DTX and CN state variables */
  /* ------------------------------------------------ */

  swRxDTXState = CNINTPER - 1;
  swDecoMode = SPEECH;
  swDtxMuting = 0;
  swDtxBfiCnt = 0;

  swOldR0IndexDec = 0;

  swRxGsHistPtr = 0;

  for (i = 0; i < (OVERHANG - 1) * N_SUB; i++)
    pL_RxGsHist[i] = 0;

}
