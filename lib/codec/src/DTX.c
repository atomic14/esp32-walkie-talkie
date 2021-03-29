/***************************************************************************
 *
 *   File Name: dtx.c
 *
 *   Purpose:   DTX and comfort noise functions of the GSM half rate
 *              system
 *
 *   Reference: Recommendation GSM 06.41 (DTX)
 *              Recommendation GSM 06.22 (Comfort Noise)
 *
 *     Below is a listing of all the functions appearing in the file.
 *     The functions are arranged according to their purpose.  Under
 *     each heading, the ordering is hierarchical.
 *
 *     Evaluation of comfort noise parameters
 *       swComfortNoise()
 *         updateCNHist()
 *         avgGsHistQntz()
 *         gsQuant()
 *         avgCNHist()
 *         lpcCorrQntz()
 *         getPnBits()
 *
 *     Interpolation of comfort noise parameters
 *        rxInterpR0Lpc()
 *         linInterpSid()
 *
 **************************************************************************/

/*________________________________________________________________________
 |                                                                        |
 |                         Include Files                                  |
 |________________________________________________________________________|
*/

#include "typedefs.h"
#include "mathhalf.h"
#include "mathdp31.h"
#include "dtx.h"
#include "sp_dec.h"
#include "sp_rom.h"
#include "sp_frm.h"

/*________________________________________________________________________
 |                                                                        |
 |                            Defines                                     |
 |________________________________________________________________________|
*/

#define PN_XOR_REG (Longword)0x00000005L
#define PN_XOR_ADD (Longword)0x40000000L

#define OH_SHIFT 3                     /* shift corresponding to OVERHANG */

#define NP_AFLAT 4
#define LPC_VQ_SEG 3

#define ASHIFT 4
#define ASCALE 0x0800


/*________________________________________________________________________
 |                                                                        |
 |                        Global Variables                                |
 |________________________________________________________________________|
*/

Shortword swVadFrmCnt = 0;             /* Indicates the number of sequential
                                        * frames where VAD == 0 */

short int siUpdPointer = 0;
Shortword swNElapsed = 50;


Longword pL_GsHist[N_SUB * (OVERHANG - 1)];


/*________________________________________________________________________
 |                                                                        |
 |                     Other External Variables                           |
 |________________________________________________________________________|
*/

extern int iLimit;

extern Shortword swR0Dec,
       swOldR0Dec,
       swR0NewCN;

extern Shortword swCNR0,
       pswCNLpc[],
       pswCNGsp0Code[],
       pswCNVSCode1[],
       pswCNVSCode2[];

/*________________________________________________________________________
 |                                                                        |
 |                         DTX Rom Tables                                 |
 |________________________________________________________________________|
*/

/* interpolation curve for comfort noise (i*1/12) i=1..12 */
Shortword psrCNNewFactor[12] = {0x0aaa, 0x1554, 0x1ffe, 0x2aa8, 0x3552,
  0x3ffc, 0x4aa6, 0x5550, 0x5ffa, 0x6aa4,
0x754e, 0x7fff};


/* Values of GS for voicing state 0, all values shifted down by 2
   shifts */
LongwordRom ppLr_gsTable[4][32] =
{
  {
    0x000011ab, 0x000038d2, 0x0000773e, 0x000144ef,
    0x00035675, 0x000648c5, 0x000c3d65, 0x0017ae17,
    0x002a3dbb, 0x005238e7, 0x00695c1a, 0x00a60d45,
    0x00e4cc68, 0x01c3ba6a, 0x019e3c96, 0x02d1fbac,
    0x030453ec, 0x0549a998, 0x05190298, 0x08258920,
    0x08daff30, 0x0c3150e0, 0x0e45d850, 0x14c111a0,
    0x0ff7e1c0, 0x18a06860, 0x13810400, 0x1abc9ee0,
    0x28500940, 0x41f22800, 0x22fc5040, 0x2cd90180
  },

  {
    0x00003ede, 0x00021fc9, 0x0013f0c3, 0x003a7be2,
    0x007a6663, 0x00fe3773, 0x012fabf4, 0x02275cd0,
    0x01c0ef14, 0x02c0b1d8, 0x0350fc70, 0x05505078,
    0x04175f30, 0x052c1098, 0x08ed3310, 0x0a63b470,
    0x05417870, 0x08995ee0, 0x07bbe018, 0x0a19fa10,
    0x0b5818c0, 0x0fd96ea0, 0x0e5cad10, 0x13b40d40,
    0x12d45840, 0x14577320, 0x2b2e5e00, 0x333e9640,
    0x194c35c0, 0x1c30f8c0, 0x2d16db00, 0x2cc970ff
  },
  {
    0x002f18e7, 0x00a47be0, 0x01222efe, 0x01c42df8,
    0x024be794, 0x03424c40, 0x036950fc, 0x04973108,
    0x038405b4, 0x05d8c8f0, 0x05063e08, 0x070cdea0,
    0x05812be8, 0x06da5fc8, 0x088fcd60, 0x0a013cb0,
    0x0909a460, 0x09e6cf40, 0x0ee581d0, 0x0ec99f20,
    0x0b4e7470, 0x0c730e80, 0x0ff39d20, 0x105d0d80,
    0x158b0b00, 0x172babe0, 0x14576460, 0x181a6720,
    0x26126e80, 0x1f590180, 0x1fdaad60, 0x2e0e8000
  },
  {
    0x00c7f603, 0x01260cda, 0x01b3926a, 0x026d82bc,
    0x0228fba0, 0x036ec5b0, 0x034bf4cc, 0x043a55d0,
    0x044f9c20, 0x05c66f50, 0x0515f890, 0x06065300,
    0x0665dc00, 0x0802b630, 0x0737a1c0, 0x087294e0,
    0x09253fc0, 0x0a619760, 0x097bd060, 0x0a6d4e50,
    0x0d19e520, 0x0e15c420, 0x0c4e4eb0, 0x0e8880e0,
    0x11cdf480, 0x12c85800, 0x10f4c0a0, 0x13e51b00,
    0x189dbaa0, 0x18a6bb60, 0x22e31500, 0x21615240
  }
};

/*************************************************************************
 *
 *   FUNCTION NAME: swComfortNoise
 *
 *   PURPOSE:
 *
 *   This routine perform the following tasks:
 *     - generation of the speech flag (swSP)
 *     - averaging and encoding of the comfort noise parameters
 *     - randomization of the codebook indices
 *
 *
 *   INPUTS:
 *
 *   swVadFrmCnt (global) - swVadFlag=0 frame counter.
 *   If swVadFlag=1 then this counter is 0, the first frame with
 *   swVadFlag=0 will set this counter to 1, with each additional
 *   swVadFlag=0 frame the counter is incremented.
 *
 *   swVadFlag - voise activity flag. swVadFlag=0 frame with
 *   no voice activity, swVadFlag=0 frame with voice activity
 *
 *   L_UnqntzdR0 - unquantized R(0), 32 bit value, output of
 *   FLAT.
 *
 *   pL_UnqntzdCorr[NP+1] - unquantized correlation sequence,
 *   also an output of FLAT.
 *
 *
 *   OUTPUTS:
 *
 *   swCNR0 - global variable, the output quantized R0 index
 *
 *   pswCNLpc[3]  - global variable, the output quantized LPC to the
 *   transmitted in the SID frame
 *
 *   pswCNGsp0Code[N_SUB] - global variable, the output quantized GSP0 indices
 *
 *   pswCNVSCode1[N_SUB] - global variable, the output quantized codevector 1
 *   indices.
 *
 *   pswCNVSCode2[N_SUB] - global variable, the output quantized codevector 2
 *   indices.
 *
 *
 *   RETURN VALUE:
 *
 *   swSP - speech flag, swSP=1 speech frames are generated, swSP=0
 *   SID frames are generated.
 *
 *************************************************************************/

Shortword swComfortNoise(Shortword swVadFlag,
                             Longword L_UnqntzdR0, Longword *pL_UnqntzdCorr)
{

/*________________________________________________________________________
 |                                                                        |
 |                        Static Variables                                |
 |________________________________________________________________________|
*/

  /* history of unquantized parameters */
  static Longword pL_R0Hist[OVERHANG];
  static Longword ppL_CorrHist[OVERHANG][NP + 1];

  /* quantized reference parameters */
  static Shortword swQntRefR0,
         swRefGsIndex;
  static int piRefVqCodewds[3];

  /* handling of short speech bursts */
  static Shortword swShortBurst;

  /* state value of random generator */
  static Longword L_TxPNSeed;

/*_________________________________________________________________________
 |                                                                         |
 |                            Automatic Variables                          |
 |_________________________________________________________________________|
*/

  Shortword swSP;
  Shortword pswFinalRc[NP];

  /* unquantized reference parameters */
  Longword L_RefR0;
  Longword pL_RefCorr[NP + 1];
  Longword L_RefGs;

  int    i;


/*_________________________________________________________________________
 |                                                                         |
 |                              Executable Code                            |
 |_________________________________________________________________________|
*/

  swSP = 1;

  /* VadFrmCnt will indicate the number of sequential frames where */
  /* swVadFlag == 0                                                */
  /* ------------------------------------------------------------- */

  if (swVadFlag)
    swVadFrmCnt = 0;                   /* Voice acitvity present */
  else
    swVadFrmCnt = add(swVadFrmCnt, 1); /* no voice activity */


  /* swNElapsed will indicate the number of frames that have elapsed */
  /* since the last SID frame with updated comfort noise parameters  */
  /* was generated                                                   */
  /* --------------------------------------------------------------- */

  swNElapsed = add(swNElapsed, 1);


  /* If no voice activity was detected.  */
  /* ----------------------------------- */

  if (swVadFrmCnt)
  {

    /* Short speech burst ? */
    /* -------------------- */

    if (swVadFrmCnt == 1)
    {
      if (sub(swNElapsed, 24) < 0)
        swShortBurst = 1;              /* short speech burst detected */
      else
        swShortBurst = 0;              /* long speech burst detected */
    }


    /* Update history, with this frames data */
    /* ------------------------------------- */

    updateCNHist(L_UnqntzdR0, pL_UnqntzdCorr,
                 pL_R0Hist, ppL_CorrHist);


    /* first SID frame */
    /* --------------- */

    if (((swShortBurst == 0) && (swVadFrmCnt == OVERHANG)) ||
        ((swShortBurst == 1) && (swVadFrmCnt == 1)))
    {

      /* init. random generator */
      /* ---------------------- */
      L_TxPNSeed = PN_INIT_SEED;


      /* average GS */
      /* ---------- */
      avgGsHistQntz(pL_GsHist, &L_RefGs);


      /* GS quantization */
      /* --------------- */
      swRefGsIndex = gsQuant(L_RefGs, 0);

    }


    /* No Overhang in case of short speech bursts,                */
    /* generate SID frames with repeated comfort noise parameters */
    /* ---------------------------------------------------------- */

    if ((swShortBurst == 1) && (swVadFrmCnt < OVERHANG))
    {

      /* generate a SID frame with repeated parameters */
      /* --------------------------------------------- */

      swSP = 0;


      /* repeat data: r0, LPC, GS */
      /* ------------------------ */

      swCNR0 = swQntRefR0;

      for (i = 0; i < 3; i++)
        pswCNLpc[i] = piRefVqCodewds[i];

      for (i = 0; i < N_SUB; i++)
        pswCNGsp0Code[i] = swRefGsIndex;

    }


    /* generate SID frames with updated comfort noise parameters */
    /* --------------------------------------------------------- */

    if (swVadFrmCnt >= OVERHANG)
    {

      /* A SID frame with updated parameters */
      /* ----------------------------------- */

      swSP = 0;
      swNElapsed = 0;


      /* average R0 and correlation values */
      /* --------------------------------- */

      avgCNHist(pL_R0Hist, ppL_CorrHist, &L_RefR0,
                pL_RefCorr);


      /* now quantize the averaged R(0) */
      /* ------------------------------ */

      swQntRefR0 = r0Quant(L_RefR0);


      /* Quantize the averaged correlation */
      /* --------------------------------- */

      lpcCorrQntz(pL_RefCorr,
                  pswFinalRc,
                  piRefVqCodewds);


      /* update frame data: r0, LPC */
      /* -------------------------- */

      swCNR0 = swQntRefR0;
      for (i = 0; i < 3; i++)
        pswCNLpc[i] = piRefVqCodewds[i];


      /* update subframe data (unvoiced mode): GSP0 */
      /* ------------------------------------------ */

      for (i = 0; i < N_SUB; i++)
        pswCNGsp0Code[i] = swRefGsIndex;

    }


    /* random codevectors */
    /* ------------------ */

    if (swSP == 0)
    {
      for (i = 0; i < N_SUB; i++)
      {
        pswCNVSCode1[i] = getPnBits(7, &L_TxPNSeed);
        pswCNVSCode2[i] = getPnBits(7, &L_TxPNSeed);
      }
    }


  }

  return (swSP);
}


/*************************************************************************
 *
 *   FUNCTION NAME:  updateCNHist
 *
 *   PURPOSE:
 *
 *     Add current frame's unquantized R(0) and LPC information to the
 *     comfort noise history, so that it will be available for
 *     averaging.
 *
 *   INPUTS:
 *
 *     Unquantized values from the coder:
 *
 *
 *     L_UnqntzdR0 - unquantized frame energy R(0), an output of FLAT
 *
 *     pL_UnqntzdCorr[NP+1] - unquantized correlation coefficient
 *     array.  Also an output of FLAT.
 *
 *     siUpdPointer (global) - A modulo counter which counts up from
 *     0 to OVERHANG-1.
 *
 *   OUTPUTS:
 *
 *     pL_R0History[OVERHANG] - history of the OVERHANG frames worth of
 *     R(0).
 *
 *     ppL_CorrHistory[OVERHANG][NP+1] - - history of the OVERHANG
 *     frames worth of pL_UnqntzdCorr[].
 *
 *   RETURN VALUE:
 *
 *     none
 *
 *************************************************************************/

void   updateCNHist(Longword L_UnqntzdR0,
                           Longword *pL_UnqntzdCorr,
                           Longword pL_R0History[],
                           Longword ppL_CorrHistory[OVERHANG][NP + 1])
{

/*_________________________________________________________________________
 |                                                                         |
 |                            Automatic Variables                          |
 |_________________________________________________________________________|
*/

  int    i;


/*_________________________________________________________________________
 |                                                                         |
 |                              Executable Code                            |
 |_________________________________________________________________________|
*/

  /* update */
  pL_R0History[siUpdPointer] = L_UnqntzdR0;

  for (i = 0; i < NP + 1; i++)
    ppL_CorrHistory[siUpdPointer][i] = pL_UnqntzdCorr[i];

  siUpdPointer = (siUpdPointer + 1) % OVERHANG;
}


/*************************************************************************
 *
 *   FUNCTION NAME: avgGsHistQntz
 *
 *   PURPOSE:
 *
 *     Average gs history, where history is of length OVERHANG-1
 *     frames.  The last frame's (i.e. this frame) gs values are not
 *     available since quantization would have occured only after the
 *     VAD decision is made.
 *
 *   INPUTS:
 *
 *     pL_GsHistory[(OVERHANG-1)*N_SUB] - the GS of the past
 *     OVERHANG-1 frames. The GS values are stored shifted down by 2
 *     shifts to avoid overflow (the largest GS is greater than 2.0).
 *
 *
 *   OUTPUTS:
 *
 *     *pL_GsAvgd - the average of pL_GsHistory[], also shifted down
 *     by two shifts.
 *
 *   RETURN VALUE:
 *
 *     none.
 *
 *
 *************************************************************************/

void   avgGsHistQntz(Longword pL_GsHistory[], Longword *pL_GsAvgd)
{

/*_________________________________________________________________________
 |                                                                         |
 |                            Automatic Variables                          |
 |_________________________________________________________________________|
*/

  int    i;
  Longword L_avg;

/*_________________________________________________________________________
 |                                                                         |
 |                              Executable Code                            |
 |_________________________________________________________________________|
*/

  L_avg = L_shift_r(pL_GsHistory[0], -(OH_SHIFT + 2));

  for (i = 1; i < N_SUB * (OVERHANG - 1); i++)
    L_avg = L_add(L_shift_r(pL_GsHistory[i], -(OH_SHIFT + 2)), L_avg);

  /* avg number x/32 not x/28 */

  *pL_GsAvgd = L_add(L_avg, L_mpy_ls(L_avg, 0x1249));   /* L_avg *= 32/28 */

}


/*************************************************************************
 *
 *   FUNCTION NAME: gsQuant
 *
 *   PURPOSE:
 *
 *     Quantize a value of gs in any of the voicing modes.  Input GS
 *     is a 32 bit number.  The GSP0 index is returned.
 *
 *   INPUTS:
 *
 *     L_GsIn - 32 bit GS value,  shifted down by 2 shifts.
 *
 *     swVoicingMode - voicing level
 *
 *     ppLr_gsTable[4][32] - Rom GS Table. (global), all GS values
 *     have been shifted down by 2 from their true value.
 *
 *   OUTPUTS:
 *
 *     none
 *
 *   RETURN VALUE:
 *
 *
 *     GSP0 Index closest to the input value of GS.
 *
 *
 *************************************************************************/

Shortword gsQuant(Longword L_GsIn, Shortword swVoicingMode)
{

/*_________________________________________________________________________
 |                                                                         |
 |                            Automatic Variables                          |
 |_________________________________________________________________________|
*/

  Shortword swGsIndex,
         swBestGs;
  Longword L_diff,
         L_min = LW_MAX;


/*_________________________________________________________________________
 |                                                                         |
 |                              Executable Code                            |
 |_________________________________________________________________________|
*/

  for (swGsIndex = 0; swGsIndex < 32; swGsIndex++)
  {
    L_diff = L_abs(L_sub(L_GsIn, ppLr_gsTable[swVoicingMode][swGsIndex]));

    if (L_sub(L_diff, L_min) < 0)
    {
      /* new minimum */
      /* ----------- */

      swBestGs = swGsIndex;
      L_min = L_diff;

    }
  }

  return (swBestGs);

}


/*************************************************************************
 *
 *   FUNCTION NAME: avgCNHist
 *
 *   PURPOSE:
 *
 *     Average the unquantized R0 and LPC data stored at the encoder
 *     to arrive at an average R0 and LPC frame for use in a SID
 *     frame.
 *
 *   INPUTS:
 *
 *   pL_R0History[OVERHANG] - contains unquantized R(0) data from the
 *   most recent OVERHANG frame (including this one).
 *
 *   ppL_CorrHistory[OVERHANG][NP+1] - Unquantized correlation
 *   coefficients from the most recent OVERHANG frame (including this
 *   one).  The data stored here is an output of FLAT.
 *
 *   OUTPUTS:
 *
 *   *pL_AvgdR0 - the average of pL_R0History[]
 *
 *   pL_AvgdCorrSeq[NP+1] - the average of ppL_CorrHistory[][].
 *
 *
 *   RETURN VALUE:
 *
 *     none
 *
 *************************************************************************/

void   avgCNHist(Longword pL_R0History[],
                        Longword ppL_CorrHistory[OVERHANG][NP + 1],
                        Longword *pL_AvgdR0,
                        Longword pL_AvgdCorrSeq[])
{

/*_________________________________________________________________________
 |                                                                         |
 |                            Automatic Variables                          |
 |_________________________________________________________________________|
*/

  int    i,
         j;
  Longword L_avg;

/*_________________________________________________________________________
 |                                                                         |
 |                              Executable Code                            |
 |_________________________________________________________________________|
*/

  /* R0 Averaging */
  /* ------------ */

  for (L_avg = 0, i = 0; i < OVERHANG; i++)
    L_avg = L_add(L_shr(pL_R0History[i], OH_SHIFT), L_avg);

  *pL_AvgdR0 = L_avg;


  /* LPC: average the last OVERHANG frames */
  /* ------------------------------------- */

  for (j = 0; j < NP + 1; j++)
  {
    for (L_avg = 0, i = 0; i < OVERHANG; i++)
    {
      L_avg = L_add(L_shift_r(ppL_CorrHistory[i][j], -OH_SHIFT), L_avg);
    }

    pL_AvgdCorrSeq[j] = L_avg;
  }

}


/***************************************************************************
 *
 *    FUNCTION NAME: lpcCorrQntz
 *
 *    PURPOSE:  Quantize a correlation sequence
 *
 *
 *    INPUT:
 *
 *         pL_CorrelSeq[NP+1]
 *                     Correlation sequence to quantize.
 *
 *    OUTPUTS:
 *
 *        pswFinalRc[0:NP-1]
 *                     A quantized set of NP reflection coefficients.
 *
 *        piVQCodewds[0:2]
 *                     An array containing the indices of the 3 reflection
 *                     coefficient vectors selected from the three segment
 *                     Rc-VQ.
 *
 *    RETURN:
 *        None.
 *
 *    KEYWORDS: AFLAT,aflat,flat,vectorquantization, reflectioncoefficients
 *
 *************************************************************************/

void   lpcCorrQntz(Longword pL_CorrelSeq[],
                          Shortword pswFinalRc[],
                          int piVQCodewds[])
{

/*_________________________________________________________________________
 |                                                                         |
 |                            Automatic Variables                          |
 |_________________________________________________________________________|
*/

  Shortword pswPOldSpace[NP_AFLAT],
         pswPNewSpace[NP_AFLAT],
         pswVOldSpace[2 * NP_AFLAT - 1],
         pswVNewSpace[2 * NP_AFLAT - 1],
        *ppswPAddrs[2],
        *ppswVAddrs[2],
        *pswVBar,
         pswPBar[NP_AFLAT],
         pswVBarSpace[2 * NP_AFLAT - 1],
         pswFlatsRc[NP],               /* Unquantized Rc's computed by FLAT */
         pswRc[NP + 1];                /* Temp list for the converted RC's */
  Longword *pL_VBarFull,
         pL_PBarFull[NP],
         pL_VBarFullSpace[2 * NP - 1];

  int    i,
         iVec,
         iSeg,
         iCnt;                         /* Loop counter */
  struct QuantList quantList,          /* A list of vectors */
         bestPql[4];                   /* The four best vectors from
                                        * the PreQ */
  struct QuantList bestQl[LPC_VQ_SEG + 1];      /* Best vectors for each of
                                                 * the three segments */

/*_________________________________________________________________________
 |                                                                         |
 |                              Executable Code                            |
 |_________________________________________________________________________|
*/

  /* Setup pointers temporary space */
  /*--------------------------------*/

  pswVBar = pswVBarSpace + NP_AFLAT - 1;
  pL_VBarFull = pL_VBarFullSpace + NP - 1;
  ppswPAddrs[0] = pswPOldSpace;
  ppswPAddrs[1] = pswPNewSpace;
  ppswVAddrs[0] = pswVOldSpace + NP_AFLAT - 1;
  ppswVAddrs[1] = pswVNewSpace + NP_AFLAT - 1;


  /* Set up pL_PBarFull and pL_VBarFull initial conditions, using the   */
  /* autocorrelation sequence derived from the optimal reflection       */
  /* coefficients computed by FLAT. The initial conditions are shifted  */
  /* right by RSHIFT bits. These initial conditions, stored as          */
  /* Longwords, are used to initialize PBar and VBar arrays for the     */
  /* next VQ segment.                                                   */
  /*--------------------------------------------------------------------*/

  initPBarFullVBarFullL(pL_CorrelSeq, pL_PBarFull, pL_VBarFull);

  /* Set up initial PBar and VBar initial conditions, using pL_PBarFull */
  /* and pL_VBarFull arrays initialized above. These are the initial    */
  /* PBar and VBar conditions to be used by the AFLAT recursion at the  */
  /* 1-st Rc-VQ segment.                                                */
  /*--------------------------------------------------------------------*/

  initPBarVBarL(pL_PBarFull, pswPBar, pswVBar);

  for (iSeg = 1; iSeg <= LPC_VQ_SEG; iSeg++)
  {
    /* initialize candidate list */
    /*---------------------------*/

    quantList.iNum = psrPreQSz[iSeg - 1];
    quantList.iRCIndex = 0;

    /* do aflat for all vectors in the list */
    /*--------------------------------------*/

    setupPreQ(iSeg, quantList.iRCIndex);        /* set up vector ptrs */

    for (iCnt = 0; iCnt < quantList.iNum; iCnt++)
    {
      /* get a vector */
      /*--------------*/

      getNextVec(pswRc);

      /* clear the limiter flag */
      /*------------------------*/

      iLimit = 0;

      /* find the error values for each vector */
      /*---------------------------------------*/

      quantList.pswPredErr[iCnt] =
              aflatRecursion(&pswRc[psvqIndex[iSeg - 1].l],
                             pswPBar, pswVBar,
                             ppswPAddrs, ppswVAddrs,
                             psvqIndex[iSeg - 1].len);

      /* check the limiter flag */
      /*------------------------*/

      if (iLimit)
        quantList.pswPredErr[iCnt] = 0x7fff;    /* set error to bad value */

    }                                  /* done list loop */

    /* find 4 best prequantizer levels */
    /*---------------------------------*/

    findBestInQuantList(quantList, 4, bestPql);

    for (iVec = 0; iVec < 4; iVec++)
    {

      /* initialize quantizer list */
      /*---------------------------*/

      quantList.iNum = psrQuantSz[iSeg - 1];
      quantList.iRCIndex = bestPql[iVec].iRCIndex * psrQuantSz[iSeg - 1];

      setupQuant(iSeg, quantList.iRCIndex);     /* set up vector ptrs */

      /* do aflat recursion on each element of list */
      /*--------------------------------------------*/

      for (iCnt = 0; iCnt < quantList.iNum; iCnt++)
      {
        /* get a vector */
        /*--------------*/

        getNextVec(pswRc);

        /* clear the limiter flag */
        /*------------------------*/

        iLimit = 0;

        /* find the error values for each vector */
        /*---------------------------------------*/

        quantList.pswPredErr[iCnt] =
                aflatRecursion(&pswRc[psvqIndex[iSeg - 1].l],
                               pswPBar, pswVBar,
                               ppswPAddrs, ppswVAddrs,
                               psvqIndex[iSeg - 1].len);

        /* check the limiter flag */
        /*------------------------*/

        if (iLimit)
          quantList.pswPredErr[iCnt] = 0x7fff;  /* set error to the worst
                                                 * value */

      }                                /* done list loop */

      /* find best quantizer vector for this segment, and save it */
      /*----------------------------------------------------------*/

      findBestInQuantList(quantList, 1, bestQl);
      if (iVec == 0)
        bestQl[iSeg] = bestQl[0];
      else if (sub(bestQl[iSeg].pswPredErr[0], bestQl[0].pswPredErr[0]) > 0)
        bestQl[iSeg] = bestQl[0];

    }

    /* find the quantized reflection coefficients */
    /*--------------------------------------------*/

    setupQuant(iSeg, bestQl[iSeg].iRCIndex);    /* set up vector ptrs */
    getNextVec((Shortword *) (pswFinalRc - 1));


    /* Update pBarFull and vBarFull for the next Rc-VQ segment, and */
    /* update the pswPBar and pswVBar for the next Rc-VQ segment    */
    /*--------------------------------------------------------------*/

    if (iSeg < LPC_VQ_SEG)
      aflatNewBarRecursionL(&pswFinalRc[psvqIndex[iSeg - 1].l - 1], iSeg,
                            pL_PBarFull, pL_VBarFull, pswPBar, pswVBar);

  }

  /* find the quantizer index (the values to be output in the symbol file) */
  /*-----------------------------------------------------------------*/

  for (iSeg = 1; iSeg <= LPC_VQ_SEG; iSeg++)
    piVQCodewds[iSeg - 1] = bestQl[iSeg].iRCIndex;

}


/*************************************************************************
 *
 *   FUNCTION NAME: getPnBits
 *
 *   PURPOSE:
 *
 *     Generate iBits pseudo-random bits using *pL_PNSeed as the
 *     pn-generators seed.
 *
 *   INPUTS:
 *
 *     iBits - integer indicating how many random bits to return.
 *     range [0,15], 0 yields 1 bit output
 *
 *     *pL_PNSeed - 32 bit seed (changed by function)
 *
 *   OUTPUTS:
 *
 *     *pL_PNSeed - 32 bit seed, modified.
 *
 *   RETURN VALUE:
 *
 *    random bits in iBits LSB's.
 *
 *
 *   IMPLEMENTATION:
 *
 *    implementation of x**31 + x**3 + 1 == PN_XOR_REG | PN_XOR_ADD a
 *    PN sequence generator using Longwords generating a 2**31 -1
 *    length pn-sequence.
 *
 *************************************************************************/

Shortword getPnBits(int iBits, Longword *pL_PNSeed)
{

/*_________________________________________________________________________
 |                                                                         |
 |                            Automatic Variables                          |
 |_________________________________________________________________________|
*/

  Shortword swPnBits = 0;
  Longword L_Taps,
         L_FeedBack;
  int    i;

/*_________________________________________________________________________
 |                                                                         |
 |                              Executable Code                            |
 |_________________________________________________________________________|
*/

  for (i = 0; i < iBits; i++)
  {
    /* update the state */
    /* ---------------- */

    L_Taps = *pL_PNSeed & PN_XOR_REG;
    L_FeedBack = L_Taps;               /* Xor tap bits to yield
                                        * feedback bit */
    L_Taps = L_shr(L_Taps, 1);

    while (L_Taps)
    {
      L_FeedBack = L_FeedBack ^ L_Taps;
      L_Taps = L_shr(L_Taps, 1);
    }

    /* LSB of L_FeedBack is next MSB of PN register */

    *pL_PNSeed = L_shr(*pL_PNSeed, 1);
    if (L_FeedBack & 1)
      *pL_PNSeed = *pL_PNSeed | PN_XOR_ADD;

    /* State update complete.  Get the output bit from the state, add/or it
     * into output */

    swPnBits = shl(swPnBits, 1);
    swPnBits = swPnBits | (extract_l(*pL_PNSeed) & 0x0001);

  }
  return (swPnBits);
}


/*************************************************************************
 *
 *   FUNCTION NAME: rxInterpR0Lpc
 *
 *   PURPOSE:
 *
 *     Perform part of the comfort noise algorithm at the decoder.
 *     LPC and R0 are derived in this routine
 *
 *   INPUTS:
 *
 *     pswOldKs - Last frame's reflection coeffs.
 *
 *     pswNewKs - This frame's decoded/received reflection coeffs.
 *     This will serve a new endpoint in interpolation.
 *
 *     swRxDTXState - primary DTX state variable (at the receiver).  A
 *     modulo 12 counter, which is 0 at SID frame.
 *
 *     swDecoMode - actual mode the decoder: speech decoding mode
 *     or comfort noise insertion mode (SPEECH = speech decoding;
 *     CNIFIRSTSID = comfort noise, 1st SID received; CNICONT = comfort
 *     noise, SID frame received, but not 1st SID; CNIBFI = comfort
 *     noise, bad frame received)
 *
 *     swFrameType - type of the received frame (VALIDSID, INVALIDSID
 *     GOODSPEECH or UNUSABLE)
 *
 *     swOldR0Dec - global variable, the decoded R0 value from the last
 *     frame .  This will be modified.
 *
 *     swR0NewCN - global variable the decoded R0 value from the frame
 *     just received. Valid information if current frame is a SID frame.
 *
 *
 *   OUTPUTS:
 *
 *     pswNewKs - This frames LPC coeffs. modified to reflect
 *     interpolated correlation sequence pL_CorrSeq[].
 *
 *     swR0Dec - global variable, interpolated R0 value
 *
 *     swR0OldCN - global variable, R0 interpolation point to
 *     interpolate from.
 *
 *     swR0NewCN - global variable, R0 interpolation point to
 *     interpolate to.
 *
 *     pL_OldCorrSeq[NP+1] - global variable, starting point for
 *     interpolation of LPC information.
 *
 *     pL_NewCorrSeq[NP+1] - global variable, end point for
 *     interpolation of LPC information.
 *
 *     pL_CorrSeq[NP+1] - global variable, interpolated value of LPC
 *     information to be used in this frame.
 *
 *
 *   RETURN VALUE:
 *
 *     None.
 *
 *   KEYWORDS: interpolation, comfort noise, SID, DTX
 *
 *************************************************************************/

void   rxInterpR0Lpc(Shortword *pswOldKs, Shortword *pswNewKs,
                            Shortword swRxDTXState,
                            Shortword swDecoMode, Shortword swFrameType)
{

/*________________________________________________________________________
 |                                                                        |
 |                        Static Variables                                |
 |________________________________________________________________________|
*/

  static Shortword swR0OldCN;
  static Longword pL_OldCorrSeq[NP + 1],
         pL_NewCorrSeq[NP + 1],
         pL_CorrSeq[NP + 1];


/*_________________________________________________________________________
 |                                                                         |
 |                            Automatic Variables                          |
 |_________________________________________________________________________|
*/

  int    i;


/*_________________________________________________________________________
 |                                                                         |
 |                              Executable Code                            |
 |_________________________________________________________________________|
*/

  if (swDecoMode == CNIFIRSTSID)
  {
    /* first SID frame arrived */
    /* ----------------------- */

    /* use tx'd R0 frame as both endpoints of interp curve. */
    /* i.e. no interpolation for the first frames           */
    /* ---------------------------------------------------- */


    swR0OldCN = swOldR0Dec;            /* last non-SID, received R0 */
    swR0Dec = linInterpSidShort(swR0NewCN, swR0OldCN, swRxDTXState);


    /* generate the LPC end points for interpolation */
    /* --------------------------------------------- */

    rcToCorrDpL(ASHIFT, ASCALE, pswOldKs, pL_OldCorrSeq);
    rcToCorrDpL(ASHIFT, ASCALE, pswNewKs, pL_NewCorrSeq);

    /* linearly interpolate between the two sets of correlation coefs */
    /* -------------------------------------------------------------- */

    for (i = 0; i < NP + 1; i++)
    {
      pL_CorrSeq[i] = linInterpSid(pL_NewCorrSeq[i], pL_OldCorrSeq[i],
                                   swRxDTXState);
    }

    /* Generate this frames K's (overwrite input) */
    /* ------------------------------------------ */

    aFlatRcDp(pL_CorrSeq, pswNewKs);

  }
  else if ((swDecoMode == CNICONT) && (swFrameType == VALIDSID))
  {
    /* new (not the first) SID frame arrived */
    /* ------------------------------------- */

    swR0OldCN = swOldR0Dec;            /* move current state of R0 to old */
    swR0Dec = linInterpSidShort(swR0NewCN, swR0OldCN, swRxDTXState);


    /* LPC: generate new endpoints for interpolation */
    /* --------------------------------------------- */

    for (i = 0; i < NP + 1; i++)
    {
      pL_OldCorrSeq[i] = pL_CorrSeq[i];
    }

    rcToCorrDpL(ASHIFT, ASCALE, pswNewKs, pL_NewCorrSeq);


    /* linearly interpolate between the two sets of correlation coefs */
    /* -------------------------------------------------------------- */

    for (i = 0; i < NP + 1; i++)
    {
      pL_CorrSeq[i] = linInterpSid(pL_NewCorrSeq[i], pL_OldCorrSeq[i],
                                   swRxDTXState);
    }


    /* Use interpolated LPC for this frame, overwrite the input K's */
    /* ------------------------------------------------------------ */

    aFlatRcDp(pL_CorrSeq, pswNewKs);

  }
  else
  {
    /* in between SID frames / invalid SID frames */
    /* ------------------------------------------ */

    swR0Dec = linInterpSidShort(swR0NewCN, swR0OldCN, swRxDTXState);


    /* linearly interpolate between the two sets of correlation coefs */
    /* -------------------------------------------------------------- */

    for (i = 0; i < NP + 1; i++)
    {
      pL_CorrSeq[i] = linInterpSid(pL_NewCorrSeq[i], pL_OldCorrSeq[i],
                                   swRxDTXState);
    }


    /* Use interpolated LPC for this frame, overwrite the input K's */
    /* ------------------------------------------------------------ */

    aFlatRcDp(pL_CorrSeq, pswNewKs);

  }
}


/*************************************************************************
 *
 *   FUNCTION NAME: linInterpSid
 *
 *   PURPOSE:
 *
 *     Linearly interpolate between two input numbers based on what the
 *     current DtxState is.
 *
 *   INPUTS:
 *
 *     L_New - longword more current value
 *
 *     L_Old - longword oldest value
 *
 *     swDtxState - state is 0 at the transmitted SID Frame.
 *
 *
 *   OUTPUTS:
 *
 *     none
 *
 *   RETURN VALUE:
 *
 *     A value between old and new inputs with dtxState+1/12 of the new
 *     (dtxState+1)-12/12 of the old
 *
 *
 *************************************************************************/

Longword linInterpSid(Longword L_New, Longword L_Old, Shortword swDtxState)
{

/*_________________________________________________________________________
 |                                                                         |
 |                            Automatic Variables                          |
 |_________________________________________________________________________|
*/

  Shortword swOldFactor;


/*_________________________________________________________________________
 |                                                                         |
 |                              Executable Code                            |
 |_________________________________________________________________________|
*/

  /* old factor = (1.0 - newFactor) */
  /* ------------------------------ */

  swOldFactor = sub(0x7fff, psrCNNewFactor[swDtxState]);
  swOldFactor = add(0x1, swOldFactor);


  /* contributions from new and old */
  /* ------------------------------ */

  L_New = L_mpy_ls(L_New, psrCNNewFactor[swDtxState]);
  L_Old = L_mpy_ls(L_Old, swOldFactor);

  return (L_add(L_New, L_Old));

}


/*************************************************************************
 *
 *   FUNCTION NAME: linInterpSidShort
 *
 *   PURPOSE:
 *
 *     Linearly interpolate between two input numbers based on what
 *     the current DtxState is.
 *
 *   INPUTS:
 *
 *     swNew - 16 bit,  more current value
 *
 *     swOld - 16 bit, oldest value
 *
 *     swDtxState - state is 0 at the transmitted SID Frame.
 *
 *
 *   OUTPUTS:
 *
 *     none
 *
 *   RETURN VALUE:
 *
 *     A value between old and new inputs with dtxState+1/12 of the new
 *     (dtxState+1)-12/12 of the old
 *
 *************************************************************************/

Shortword linInterpSidShort(Shortword swNew, Shortword swOld,
                                   Shortword swDtxState)
{

/*_________________________________________________________________________
 |                                                                         |
 |                            Automatic Variables                          |
 |_________________________________________________________________________|
*/

  Shortword swOldFactor;
  Longword L_New,
         L_Old;


/*_________________________________________________________________________
 |                                                                         |
 |                              Executable Code                            |
 |_________________________________________________________________________|
*/

  /* old factor = (1.0 - newFactor) */
  /* ------------------------------ */

  swOldFactor = sub(0x7fff, psrCNNewFactor[swDtxState]);
  swOldFactor = add(0x1, swOldFactor);


  /* contributions from new and old */
  /* ------------------------------ */

  L_New = L_mult(swNew, psrCNNewFactor[swDtxState]);
  L_Old = L_mult(swOld, swOldFactor);


  return (round(L_add(L_New, L_Old)));

}
