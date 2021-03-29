/***************************************************************************
 *
 *   File Name:  sp_frm.c
 *
 *   Purpose:  Contains all functions for frame-based processing in the
 *      speech encoder.  The frame-based processing yields the following:
 *      energy in the speech signal, LPC filter coefficients, perceptually-
 *      weighted filter coefficients (for H(z) and C(z)), perceptually-
 *      weighted speech, voicing level, and constrained adaptive-codebook
 *      (long-term predictor) choices.
 *
 *     Below is a listing of all the functions appearing in the file.
 *     The functions are arranged according to their purpose.  Under
 *     each heading, the ordering is hierarchical.
 *
 *     High pass filtering:
 *       filt4_2nd()
 *         iir_d()
 *
 *     AFLAT, vector quantization of LPC coefficients:
 *       aflat()
 *         aflatNewBarRecursionL()
 *         aflatRecursion()
 *         findBestInQuantList()
 *         getNextVec()
 *         initPBarVBarFullL()
 *         initPBarVBarL()
 *         setupPreQ()
 *         setupQuant()
 *
 *     FLAT:  derivation of the unquantized LPC coefficients:
 *       flat()
 *         cov32()
 *         r0Quant()
 *
 *
 *     Generation of LPC filters for each subframe:
 *       getSfrmLpcTx()
 *         compResidEnergy()
 *
 *     Perceptual weighting:
 *       weightSpeechFrame()
 *
 *     Generation of the noise weighting filter:
 *       getNWCoefs()
 *
 *     Open loop lag search:
 *       openLoopLagSearch()
 *         bestDelta()
 *           maxCCOverGWithSign()
 *         getCCThreshold()
 *           fnExp2()
 *           fnLog2()
 *         pitchLags()
 *           CGInterp()
 *           CGInterpValid()
 *           findPeak()
 *           fnBest_CG()
 *           quantLag()
 *
 **************************************************************************/

/*_________________________________________________________________________
 |                                                                         |
 |                            Include Files                                |
 |_________________________________________________________________________|
*/

#include "mathhalf.h"
#include "mathdp31.h"
#include "sp_rom.h"
#include "sp_dec.h"
#include "sp_frm.h"
#include "sp_sfrm.h"
#include "vad.h"
#include "dtx.h"

/*_________________________________________________________________________
 |                                                                         |
 |                            Local Constants                              |
 |_________________________________________________________________________|
*/

#define ASCALE  0x0800
#define ASHIFT 4
#define CG_INT_MACS     6
#define CG_TERMS        (LSMAX - LSMIN + 1)
#define CVSHIFT 2                      /* Number of right shifts to be
                                        * applied to the normalized Phi
                                        * array in cov32, also used in flat
                                        * to shift down normalized F, B, C
                                        * matrices.                        */
#define C_FRAME_LEN     (N_SUB * CG_TERMS)
#define DELTA_LEVELS    16
#define G_FRAME_LEN     (LSMAX + (N_SUB-1) * S_LEN - LSMIN  + 1)
#define HIGH 1
#define INV_OS_FCTR     0x1555         /* 1.0/6.0 */
#define LAG_TABLE_LEN   (1 << L_BITS)
#define LMAX            142
#define LMAX_FR         (LMAX * OS_FCTR)
#define LMIN            21
#define LMIN_FR         (LMIN * OS_FCTR)
#define LOW 0
#define LPC_VQ_SEG 3
#define LSMAX           (LMAX + CG_INT_MACS/2)
#define LSMIN           (LMIN - CG_INT_MACS/2)
#define LSP_MASK  0xffff
#define L_BITS          8
#define L_ROUND (Longword)0x8000       /* Preload accumulator value for
                                        * rounding  */
#define NP_AFLAT     4
#define NUM_CLOSED      3
#define NUM_TRAJ_MAX    2
#define ONE_EIGHTH      0x1000
#define ONE_HALF        0x4000
#define ONE_QUARTER     0x2000
#define PEAK_VICINITY   3
#define PGAIN_CLAMP    0x0021          /* 0.001 */
#define PGAIN_SCALE    0x6000          /* 0.75 */
#define PW_FRAC         0x3333         /* 0.4 */
#define R0BITS 5
#define RSHIFT  2
#define S_SH    6                      /* Shift offset for computing frame
                                        * energy */
#define UV_SCALE0       -0x2976
#define UV_SCALE1       -0x46d3
#define UV_SCALE2       -0x6676
#define W_F_BUFF_LEN  (F_LEN + LSMAX)
#define high(x) (shr(x,8) & 0x00ff)
#define low(x) x & 0x00ff              /* This macro will return the low
                                        * byte of a word */
#define odd(x) (x & 0x0001)            /* This macro will determine if an
                                        * integer is odd */

/*_________________________________________________________________________
 |                                                                         |
 |                         State Variables (globals)                       |
 |_________________________________________________________________________|
*/

Shortword pswAnalysisState[NP];

Shortword pswWStateNum[NP],
       pswWStateDenom[NP];

/*_________________________________________________________________________
 |                                                                         |
 |                         Other External Variables                        |
 |_________________________________________________________________________|
*/

static ShortwordRom *psrTable;         /* points to correct table of
                                        * vectors */
int iLimit;                            /* accessible to all in this file
                                        * and to lpcCorrQntz() in dtx.c */
static int iLow;                       /* the low element in this segment */
static int iThree;                     /* boolean, is this a three element
                                        * vector */
static int iWordHalfPtr;               /* points to the next byte */
static int iWordPtr;                   /* points to the next word to be
                                        * read */

extern Shortword pswCNVSCode1[],       /* comfort noise parameters */
                 pswCNVSCode2[],
                 pswCNGsp0Code[],
                 pswCNLpc[],
                 swCNR0;

/***************************************************************************
 *
 *    FUNCTION NAME: aflat
 *
 *    PURPOSE:  Given a vector of high-pass filtered input speech samples
 *              (A_LEN samples), function aflat computes the NP unquantized
 *              reflection coefficients using the FLAT algorithm, searches
 *              the three segment Rc-VQ based on the AFLAT recursion, and
 *              outputs a quantized set of NP reflection coefficients, along
 *              with the three indices specifying the selected vectors
 *              from the Rc-VQ. The index of the quantized frame energy R0
 *              is also output.
 *
 *
 *    INPUT:
 *
 *        pswSpeechToLpc[0:A_LEN-1]
 *                     A vector of high-pass filtered input speech, from
 *                     which the unquantized reflection coefficients and
 *                     the index of the quantized frame energy are
 *                     computed.
 *
 *    OUTPUTS:
 *
 *        piR0Index[0:0]
 *                     An index into a 5 bit table of quantized frame
 *                     energies.
 *
 *        pswFinalRc[0:NP-1]
 *                     A quantized set of NP reflection coefficients.
 *
 *        piVQCodewds[0:2]
 *                     An array containing the indices of the 3 reflection
 *                     coefficient vectors selected from the three segment
 *                     Rc-VQ.
 *
 *        swPtch
 *                     Flag to indicate a periodic signal component
 *
 *        pswVadFlag
 *                     Voice activity decision flag
 *                      = 1: voice activity
 *                      = 0: no voice activity
 *
 *        pswSP
 *                     Speech flag
 *                      = 1: encoder generates speech frames
 *                      = 0: encoder generate SID frames
 *
 *
 *    RETURN:
 *        None.
 *
 *    REFERENCE:  Sub-clauses 4.1.3, 4.1.4, and 4.1.4.1
 *        of GSM Recommendation 06.20
 *
 *    KEYWORDS: AFLAT,aflat,flat,vectorquantization, reflectioncoefficients
 *
 *************************************************************************/

  void   aflat(Shortword pswSpeechToLPC[],
                      int piR0Index[],
                      Shortword pswFinalRc[],
                      int piVQCodewds[],
                      Shortword swPtch,
                      Shortword *pswVadFlag,
                      Shortword *pswSP)
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
  Longword pL_CorrelSeq[NP + 1],
        *pL_VBarFull,
         pL_PBarFull[NP],
         pL_VBarFullSpace[2 * NP - 1];

  int    i,
         iVec,
         iSeg,
         iCnt;                         /* Loop counter */
  struct QuantList quantList,          /* A list of vectors */
         bestPql[4];                   /* The four best vectors from the
                                        * PreQ */
  struct QuantList bestQl[LPC_VQ_SEG + 1];      /* Best vectors for each of
                                                 * the three segments */
  Shortword swVadScalAuto;
  Shortword pswVadRc[4];
  Longword pL_VadAcf[9];

  Longword L_R0;                       /* Normalized R0 (use swRShifts to 
                                        * unnormalize). This is done prior
                                        * to r0quant(). After this, its is
                                        * a unnormalized number */

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

  /* Given the input speech, compute the optimal reflection coefficients */
  /* using the FLAT algorithm.                                           */
  /*---------------------------------------------------------------------*/

  L_R0 = flat(pswSpeechToLPC, pswFlatsRc, piR0Index, pL_VadAcf, 
              &swVadScalAuto);

  /* Get unquantized reflection coefficients for VAD */      /* DTX mode */
  /* algorithm                                       */      /* DTX mode */
  /* ----------------------------------------------- */      /* DTX mode */

  for (i = 0; i < 4; i++)                                    /* DTX mode */
    pswVadRc[i] = pswFlatsRc[i];                             /* DTX mode */


  /* convert reflection coefficients to correlation */       /* DTX mode */
  /* sequence                                       */       /* DTX mode */
  /* ---------------------------------------------- */       /* DTX mode */

  rcToCorrDpL(ASHIFT, ASCALE, pswFlatsRc, pL_CorrelSeq);     /* DTX mode */


  /* Make the voice activity detection. Only swVadFlag is */ /* DTX mode */
  /*  modified.                                           */ /* DTX mode */
  /* ---------------------------------------------------- */ /* DTX mode */

  vad_algorithm(pL_VadAcf, swVadScalAuto, pswVadRc, swPtch,  /* DTX mode */
                pswVadFlag);


  /* if DTX mode off, then always voice activity */          /* DTX mode */
  /* ------------------------------------------- */          /* DTX mode */
  if (!giDTXon) *pswVadFlag = 1;                             /* DTX mode */


  /* determination of comfort noise parameters */            /* DTX mode */
  /* ----------------------------------------- */            /* DTX mode */

  *pswSP = swComfortNoise(*pswVadFlag,                       /* DTX mode */
                          L_R0,                              /* DTX mode */
                          pL_CorrelSeq);                     /* DTX mode */

  if (*pswSP == 0)                                           /* DTX mode */
  {   /* SID frame generation */                             /* DTX mode */

    /* use unquantized reflection coefficients in the */     /* DTX mode */
    /* encoder, when SID frames are generated         */     /* DTX mode */
    /* ---------------------------------------------- */     /* DTX mode */

    for (i = 0; i < NP; i++)                                 /* DTX mode */
      pswFinalRc[i] = pswFlatsRc[i];                         /* DTX mode */

  }                                                          /* DTX mode */
  else                                                       /* DTX mode */
  { /* speech frame generation */

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
        {
          quantList.pswPredErr[iCnt] = 0x7fff;    /* set error to bad value */
        }

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
          {
            quantList.pswPredErr[iCnt] = 0x7fff;  /* set error to the worst
                                                   * value */
          }

        }                                /* done list loop */

        /* find best quantizer vector for this segment, and save it */
        /*----------------------------------------------------------*/

        findBestInQuantList(quantList, 1, bestQl);
        if (iVec == 0)
        {
          bestQl[iSeg] = bestQl[0];
        }
        else
        {
          if (sub(bestQl[iSeg].pswPredErr[0],
                  bestQl[0].pswPredErr[0]) > 0)
          {
            bestQl[iSeg] = bestQl[0];
          }
        }
      }

      /* find the quantized reflection coefficients */
      /*--------------------------------------------*/

      setupQuant(iSeg, bestQl[iSeg].iRCIndex);    /* set up vector ptrs */
      getNextVec((Shortword *) (pswFinalRc - 1));


      /* Update pBarFull and vBarFull for the next Rc-VQ segment, and */
      /* update the pswPBar and pswVBar for the next Rc-VQ segment    */
      /*--------------------------------------------------------------*/

      if (iSeg < LPC_VQ_SEG)
      {

        aflatNewBarRecursionL(&pswFinalRc[psvqIndex[iSeg - 1].l - 1], iSeg,
                              pL_PBarFull, pL_VBarFull, pswPBar, pswVBar);

      }

    }

    /* find the quantizer index (the values */
    /* to be output in the symbol file)     */
    /*--------------------------------------*/

    for (iSeg = 1; iSeg <= LPC_VQ_SEG; iSeg++)
    {
      piVQCodewds[iSeg - 1] = bestQl[iSeg].iRCIndex;
    }

  }
  
}

/***************************************************************************
 *
 *    FUNCTION NAME: aflatNewBarRecursionL
 *
 *    PURPOSE:  Given the Longword initial condition arrays, pL_PBarFull and
 *              pL_VBarFull, a reflection coefficient vector selected from
 *              the Rc-VQ at the current stage, and index of the current
 *              Rc-VQ stage, the AFLAT recursion is evaluated to obtain the
 *              updated initial conditions for the AFLAT recursion at the
 *              next Rc-VQ stage. At each lattice stage the pL_PBarFull and
 *              pL_VBarFull arrays are shifted to be RSHIFT down from full
 *              scale. Two sets of initial conditions are output:
 *
 *              1) pswPBar and pswVBar Shortword arrays are used at the
 *                 next Rc-VQ segment as the AFLAT initial conditions
 *                 for the Rc prequantizer and the Rc quantizer searches.
 *              2) pL_PBarFull and pL_VBarFull arrays are output and serve
 *                 as the initial conditions for the function call to
 *                 aflatNewBarRecursionL at the next lattice stage.
 *
 *
 *              This is an implementation of equations 4.24 through
 *              4.27.
 *    INPUTS:
 *
 *        pswQntRc[0:NP_AFLAT-1]
 *                     An input reflection coefficient vector selected from
 *                     the Rc-VQ quantizer at the current stage.
 *
 *        iSegment
 *                    An input describing the current Vector quantizer
 *                    quantizer segment (1, 2, or 3).
 *
 *        RSHIFT      The number of shifts down from full scale the
 *                     pL_PBarFull and pL_VBarFull arrays are to be shifted
 *                     at each lattice stage. RSHIFT is a global constant.
 *
 *        pL_PBar[0:NP-1]
 *                     A Longword input array containing the P initial
 *                     conditions for the full 10-th order LPC filter.
 *                     The address of the 0-th element of  pL_PBarFull
 *                     is passed in when function aflatNewBarRecursionL
 *                     is called.
 *
 *        pL_VBar[-NP+1:NP-1]
 *                     A Longword input array containing the V initial
 *                     conditions for the full 10-th order LPC filter.
 *                     The address of the 0-th element of  pL_VBarFull
 *                     is passed in when function aflatNewBarRecursionL
 *                     is called.
 *
 *    OUTPUTS:
 *
 *        pL_PBar[0:NP-1]
 *                     A Longword output array containing the updated P
 *                     initial conditions for the full 10-th order LPC
 *                     filter.
 *
 *        pL_VBar[-NP+1:NP-1]
 *                     A Longword output array containing the updated V
 *                     initial conditions for the full 10-th order LPC
 *                     filter.
 *
 *        pswPBar[0:NP_AFLAT-1]
 *                     An output Shortword array containing the P initial
 *                     conditions for the P-V AFLAT recursion for the next
 *                     Rc-VQ segment. The address of the 0-th element of
 *                     pswVBar is passed in.
 *
 *        pswVBar[-NP_AFLAT+1:NP_AFLAT-1]
 *                     The output Shortword array containing the V initial
 *                     conditions for the P-V AFLAT recursion, for the next
 *                     Rc-VQ segment. The address of the 0-th element of
 *                     pswVBar is passed in.
 *
 *    RETURN:
 *        None.
 *
 *    REFERENCE:  Sub-clause 4.1.4.1 GSM Recommendation 06.20
 *
 *************************************************************************/

void   aflatNewBarRecursionL(Shortword pswQntRc[], int iSegment,
                                    Longword pL_PBar[], Longword pL_VBar[],
                                 Shortword pswPBar[], Shortword pswVBar[])
{

/*_________________________________________________________________________
 |                                                                         |
 |                            Automatic Variables                          |
 |_________________________________________________________________________|
*/
  Longword *pL_VOld,
        *pL_VNew,
        *pL_POld,
        *pL_PNew,
        *ppL_PAddrs[2],
        *ppL_VAddrs[2],
         pL_VOldSpace[2 * NP - 1],
         pL_VNewSpace[2 * NP - 1],
         pL_POldSpace[NP],
         pL_PNewSpace[NP],
         L_temp,
         L_sum;
  Shortword swQntRcSq,
         swNShift;
  short int i,
         j,
         bound;

/*_________________________________________________________________________
 |                                                                         |
 |                              Executable Code                            |
 |_________________________________________________________________________|
*/
  /* Copy the addresses of the input PBar and VBar arrays into  */
  /* pL_POld and pL_VOld respectively.                          */
  /*------------------------------------------------------------*/

  pL_POld = pL_PBar;
  pL_VOld = pL_VBar;

  /* Point to PNew and VNew temporary arrays */
  /*-----------------------------------------*/

  pL_PNew = pL_PNewSpace;
  pL_VNew = pL_VNewSpace + NP - 1;

  /* Load the addresses of the temporary buffers into the address arrays. */
  /* The address arrays are used to swap PNew and POld (VNew and VOLd)    */
  /* buffers to avoid copying of the buffer contents at the end of a      */
  /* lattice filter stage.                                                */
  /*----------------------------------------------------------------------*/

  ppL_PAddrs[0] = pL_POldSpace;
  ppL_PAddrs[1] = pL_PNewSpace;
  ppL_VAddrs[0] = pL_VOldSpace + NP - 1;
  ppL_VAddrs[1] = pL_VNewSpace + NP - 1;


  /* Update AFLAT recursion initial conditions for searching the Rc vector */
  /* quantizer at the next VQ segment.                                     */
  /*-------------------------------------------------------------------*/

  for (j = 0; j < psvqIndex[iSegment - 1].len; j++)
  {
    bound = NP - psvqIndex[iSegment - 1].l - j - 1;

    /* Compute rc squared, used by the recursion at the j-th lattice stage. */
    /*---------------------------------------------------------------------*/

    swQntRcSq = mult_r(pswQntRc[j], pswQntRc[j]);

    /* Calculate PNew(i) */
    /*-------------------*/

    L_temp = L_mpy_ls(pL_VOld[0], pswQntRc[j]);
    L_sum = L_add(L_temp, pL_POld[0]);
    L_temp = L_mpy_ls(pL_POld[0], swQntRcSq);
    L_sum = L_add(L_temp, L_sum);
    L_temp = L_mpy_ls(pL_VOld[0], pswQntRc[j]);
    L_temp = L_add(L_sum, L_temp);

    /* Compute the number of bits to shift left by to achieve  */
    /* the nominal value of PNew[0] which is right shifted by  */
    /* RSHIFT bits relative to full scale.                     */
    /*---------------------------------------------------------*/

    swNShift = sub(norm_s(extract_h(L_temp)), RSHIFT);

    /* Rescale PNew[0] by shifting left by swNShift bits */
    /*---------------------------------------------------*/

    pL_PNew[0] = L_shl(L_temp, swNShift);

    for (i = 1; i <= bound; i++)
    {
      L_temp = L_mpy_ls(pL_VOld[i], pswQntRc[j]);
      L_sum = L_add(L_temp, pL_POld[i]);
      L_temp = L_mpy_ls(pL_POld[i], swQntRcSq);
      L_sum = L_add(L_temp, L_sum);
      L_temp = L_mpy_ls(pL_VOld[-i], pswQntRc[j]);
      L_temp = L_add(L_sum, L_temp);
      pL_PNew[i] = L_shl(L_temp, swNShift);
    }

    /* Calculate VNew(i) */
    /*-------------------*/

    for (i = -bound; i < 0; i++)
    {
      L_temp = L_mpy_ls(pL_VOld[-i - 1], swQntRcSq);
      L_sum = L_add(L_temp, pL_VOld[i + 1]);
      L_temp = L_mpy_ls(pL_POld[-i - 1], pswQntRc[j]);
      L_temp = L_shl(L_temp, 1);
      L_temp = L_add(L_temp, L_sum);
      pL_VNew[i] = L_shl(L_temp, swNShift);
    }
    for (i = 0; i <= bound; i++)
    {
      L_temp = L_mpy_ls(pL_VOld[-i - 1], swQntRcSq);
      L_sum = L_add(L_temp, pL_VOld[i + 1]);
      L_temp = L_mpy_ls(pL_POld[i + 1], pswQntRc[j]);
      L_temp = L_shl(L_temp, 1);
      L_temp = L_add(L_temp, L_sum);
      pL_VNew[i] = L_shl(L_temp, swNShift);
    }

    if (j < psvqIndex[iSegment - 1].len - 2)
    {

      /* Swap POld and PNew buffers, using modulo addressing */
      /*-----------------------------------------------------*/

      pL_POld = ppL_PAddrs[(j + 1) % 2];
      pL_PNew = ppL_PAddrs[j % 2];

      /* Swap VOld and VNew buffers, using modulo addressing */
      /*-----------------------------------------------------*/

      pL_VOld = ppL_VAddrs[(j + 1) % 2];
      pL_VNew = ppL_VAddrs[j % 2];
    }
    else
    {
      if (j == psvqIndex[iSegment - 1].len - 2)
      {

        /* Then recursion to be done for one more lattice stage */
        /*------------------------------------------------------*/

        /* Copy address of PNew into POld */
        /*--------------------------------*/
        pL_POld = ppL_PAddrs[(j + 1) % 2];

        /* Copy address of the input pL_PBar array into pswPNew; this will */
        /* cause the PNew array to overwrite the input pL_PBar array, thus */
        /* updating it at the final lattice stage of the current segment   */
        /*-----------------------------------------------------------------*/

        pL_PNew = pL_PBar;

        /* Copy address of VNew into VOld */
        /*--------------------------------*/

        pL_VOld = ppL_VAddrs[(j + 1) % 2];

        /* Copy address of the input pL_VBar array into pswVNew; this will */
        /* cause the VNew array to overwrite the input pL_VBar array, thus */
        /* updating it at the final lattice stage of the current segment   */
        /*-----------------------------------------------------------------*/

        pL_VNew = pL_VBar;

      }
    }
  }

  /* Update the pswPBar and pswVBar initial conditions for the AFLAT      */
  /* Rc-VQ search at the next segment.                                    */
  /*----------------------------------------------------------------------*/

  bound = psvqIndex[iSegment].len - 1;

  for (i = 0; i <= bound; i++)
  {
    pswPBar[i] = round(pL_PBar[i]);
    pswVBar[i] = round(pL_VBar[i]);
  }
  for (i = -bound; i < 0; i++)
  {
    pswVBar[i] = round(pL_VBar[i]);
  }

  return;
}

/***************************************************************************
 *
 *    FUNCTION NAME: aflatRecursion
 *
 *    PURPOSE:  Given the Shortword initial condition arrays, pswPBar and
 *              pswVBar, a reflection coefficient vector from the quantizer
 *              (or a prequantizer), and the order of the current Rc-VQ
 *              segment, function aflatRecursion computes and returns the
 *              residual error energy by evaluating the AFLAT recursion.
 *
 *              This is an implementation of equations 4.18 to 4.23.
 *    INPUTS:
 *
 *        pswQntRc[0:NP_AFLAT-1]
 *                     An input reflection coefficient vector from the
 *                     Rc-prequantizer or the Rc-VQ codebook.
 *
 *        pswPBar[0:NP_AFLAT-1]
 *                     The input Shortword array containing the P initial
 *                     conditions for the P-V AFLAT recursion at the current
 *                     Rc-VQ segment. The address of the 0-th element of
 *                     pswVBar is passed in.
 *
 *        pswVBar[-NP_AFLAT+1:NP_AFLAT-1]
 *                     The input Shortword array containing the V initial
 *                     conditions for the P-V AFLAT recursion, at the current
 *                     Rc-VQ segment. The address of the 0-th element of
 *                     pswVBar is passed in.
 *
 *        *ppswPAddrs[0:1]
 *                     An input array containing the address of temporary
 *                     space P1 in its 0-th element, and the address of
 *                     temporary space P2 in its 1-st element. Each of
 *                     these addresses is alternately assigned onto
 *                     pswPNew and pswPOld pointers using modulo
 *                     arithmetic, so as to avoid copying the contents of
 *                     pswPNew array into the pswPOld array at the end of
 *                     each lattice stage of the AFLAT recursion.
 *                     Temporary space P1 and P2 is allocated outside
 *                     aflatRecursion by the calling function aflat.
 *
 *        *ppswVAddrs[0:1]
 *                     An input array containing the address of temporary
 *                     space V1 in its 0-th element, and the address of
 *                     temporary space V2 in its 1-st element. Each of
 *                     these addresses is alternately assigned onto
 *                     pswVNew and pswVOld pointers using modulo
 *                     arithmetic, so as to avoid copying the contents of
 *                     pswVNew array into the pswVOld array at the end of
 *                     each lattice stage of the AFLAT recursion.
 *                     Temporary space V1 and V2 is allocated outside
 *                     aflatRecursion by the calling function aflat.
 *
 *        swSegmentOrder
 *                     This input short word describes the number of
 *                     stages needed to compute the vector
 *                     quantization of the given segment.
 *
 *    OUTPUTS:
 *        None.
 *
 *    RETURN:
 *        swRe         The Shortword value of residual energy for the
 *                     Rc vector, given the pswPBar and pswVBar initial
 *                     conditions.
 *
 *    REFERENCE:  Sub-clause 4.1.4.1 GSM Recommendation 06.20
 *
 *************************************************************************/

Shortword aflatRecursion(Shortword pswQntRc[],
                                Shortword pswPBar[],
                                Shortword pswVBar[],
                                Shortword *ppswPAddrs[],
                                Shortword *ppswVAddrs[],
                                Shortword swSegmentOrder)
{

/*_________________________________________________________________________
 |                                                                         |
 |                            Automatic Variables                          |
 |_________________________________________________________________________|
*/

  Shortword *pswPOld,
        *pswPNew,
        *pswVOld,
        *pswVNew,
         pswQntRcSqd[NP_AFLAT],
         swRe;
  Longword L_sum;
  short int i,
         j,
         bound;                        /* loop control variables */

/*_________________________________________________________________________
 |                                                                         |
 |                              Executable Code                            |
 |_________________________________________________________________________|
*/

  /* Point to PBar and VBar, the initial condition arrays for the AFLAT  */
  /* recursion.                                                          */
  /*---------------------------------------------------------------------*/

  pswPOld = pswPBar;
  pswVOld = pswVBar;

  /* Point to PNew and VNew, the arrays into which updated values of  P  */
  /* and V functions will be written.                                    */
  /*---------------------------------------------------------------------*/

  pswPNew = ppswPAddrs[1];
  pswVNew = ppswVAddrs[1];

  /* Compute the residual error energy due to the selected Rc vector */
  /* using the AFLAT recursion.                                      */
  /*-----------------------------------------------------------------*/

  /* Compute rc squared, used by the recursion */
  /*-------------------------------------------*/

  for (j = 0; j < swSegmentOrder; j++)
  {
    pswQntRcSqd[j] = mult_r(pswQntRc[j], pswQntRc[j]);
  }

  /* Compute the residual error energy due to the selected Rc vector */
  /* using the AFLAT recursion.                                      */
  /*-----------------------------------------------------------------*/

  for (j = 0; j < swSegmentOrder - 1; j++)
  {
    bound = swSegmentOrder - j - 2;

    /* Compute Psubj(i), for i = 0, bound  */
    /*-------------------------------------*/

    for (i = 0; i <= bound; i++)
    {
      L_sum = L_mac(L_ROUND, pswVOld[i], pswQntRc[j]);
      L_sum = L_mac(L_sum, pswVOld[-i], pswQntRc[j]);
      L_sum = L_mac(L_sum, pswPOld[i], pswQntRcSqd[j]);
      L_sum = L_msu(L_sum, pswPOld[i], SW_MIN);
      pswPNew[i] = extract_h(L_sum);
    }

    /* Check if potential for limiting exists. */
    /*-----------------------------------------*/

    if (sub(pswPNew[0], 0x4000) >= 0)
      iLimit = 1;

    /* Compute the new Vsubj(i) */
    /*--------------------------*/

    for (i = -bound; i < 0; i++)
    {
      L_sum = L_msu(L_ROUND, pswVOld[i + 1], SW_MIN);
      L_sum = L_mac(L_sum, pswQntRcSqd[j], pswVOld[-i - 1]);
      L_sum = L_mac(L_sum, pswQntRc[j], pswPOld[-i - 1]);
      L_sum = L_mac(L_sum, pswQntRc[j], pswPOld[-i - 1]);
      pswVNew[i] = extract_h(L_sum);
    }

    for (i = 0; i <= bound; i++)
    {
      L_sum = L_msu(L_ROUND, pswVOld[i + 1], SW_MIN);
      L_sum = L_mac(L_sum, pswQntRcSqd[j], pswVOld[-i - 1]);
      L_sum = L_mac(L_sum, pswQntRc[j], pswPOld[i + 1]);
      L_sum = L_mac(L_sum, pswQntRc[j], pswPOld[i + 1]);
      pswVNew[i] = extract_h(L_sum);
    }

    if (j < swSegmentOrder - 2)
    {

      /* Swap POld and PNew buffers, using modulo addressing */
      /*-----------------------------------------------------*/

      pswPOld = ppswPAddrs[(j + 1) % 2];
      pswPNew = ppswPAddrs[j % 2];

      /* Swap VOld and VNew buffers, using modulo addressing */
      /*-----------------------------------------------------*/

      pswVOld = ppswVAddrs[(j + 1) % 2];
      pswVNew = ppswVAddrs[j % 2];

    }
  }

  /* Computing Psubj(0) for the last lattice stage */
  /*-----------------------------------------------*/

  j = swSegmentOrder - 1;

  L_sum = L_mac(L_ROUND, pswVNew[0], pswQntRc[j]);
  L_sum = L_mac(L_sum, pswVNew[0], pswQntRc[j]);
  L_sum = L_mac(L_sum, pswPNew[0], pswQntRcSqd[j]);
  L_sum = L_msu(L_sum, pswPNew[0], SW_MIN);
  swRe = extract_h(L_sum);

  /* Return the residual energy corresponding to the reflection   */
  /* coefficient vector being evaluated.                          */
  /*--------------------------------------------------------------*/

  return (swRe);                       /* residual error is returned */

}

/***************************************************************************
 *
 *   FUNCTION NAME: bestDelta
 *
 *   PURPOSE:
 *
 *     This function finds the delta-codeable lag which maximizes CC/G.
 *
 *   INPUTS:
 *
 *     pswLagList[0:siNumLags-1]
 *
 *                     List of delta-codeable lags over which search is done.
 *
 *     pswCSfrm[0:127]
 *
 *                     C(k) sequence, k integer.
 *
 *     pswGSfrm[0:127]
 *
 *                     G(k) sequence, k integer.
 *
 *     siNumLags
 *
 *                     Number of lags in contention.
 *
 *     siSfrmIndex
 *
 *                     The index of the subframe to which the delta-code
 *                     applies.
 *
 *
 *   OUTPUTS:
 *
 *     pswLTraj[0:3]
 *
 *                     The winning lag is put into this array at
 *                     pswLTraj[siSfrmIndex]
 *
 *     pswCCTraj[0:3]
 *
 *                     The corresponding winning C**2 is put into this
 *                     array at pswCCTraj[siSfrmIndex]
 *
 *     pswGTraj[0:3]
 *
 *                     The corresponding winning G is put into this arrray
 *                     at pswGTraj[siSfrmIndex]
 *
 *   RETURN VALUE:
 *
 *     none
 *
 *   DESCRIPTION:
 *
 *   REFERENCE:  Sub-clause 4.1.8.3 of GSM Recommendation 06.20
 *
 *   KEYWORDS:
 *
 *************************************************************************/

void   bestDelta(Shortword pswLagList[],
                        Shortword pswCSfrm[],
                        Shortword pswGSfrm[],
                        short int siNumLags,
                        short int siSfrmIndex,
                        Shortword pswLTraj[],
                        Shortword pswCCTraj[],
                        Shortword pswGTraj[])
{

/*_________________________________________________________________________
 |                                                                         |
 |                            Automatic Variables                          |
 |_________________________________________________________________________|
*/

  Shortword pswCBuf[DELTA_LEVELS + CG_INT_MACS + 2],
         pswGBuf[DELTA_LEVELS + CG_INT_MACS + 2],
         pswCInterp[DELTA_LEVELS + 2],
         pswGInterp[DELTA_LEVELS + 2],
        *psw1,
        *psw2,
         swCmaxSqr,
         swGmax,
         swPeak;
  short int siIPLo,
         siRemLo,
         siIPHi,
         siRemHi,
         siLoLag,
         siHiLag,
         siI;

/*_________________________________________________________________________
 |                                                                         |
 |                              Executable Code                            |
 |_________________________________________________________________________|
*/

  /* get bounds for integer C's and G's needed for interpolation */
  /* get integer and fractional portions of boundary lags        */
  /* ----------------------------------------------------------- */

  get_ipjj(pswLagList[0], &siIPLo, &siRemLo);

  get_ipjj(pswLagList[siNumLags - 1], &siIPHi, &siRemHi);

  /* get lag for first and last C and G required */
  /* ------------------------------------------- */

  siLoLag = sub(siIPLo, CG_INT_MACS / 2 - 1);

  if (siRemHi != 0)
  {
    siHiLag = add(siIPHi, CG_INT_MACS / 2);
  }
  else
  {
    siHiLag = add(siIPHi, CG_INT_MACS / 2 - 1);
  }

  /* transfer needed integer C's and G's to temp buffers */
  /* --------------------------------------------------- */

  psw1 = pswCBuf;
  psw2 = pswGBuf;

  if (siRemLo == 0)
  {

    /* first lag in list is integer: don't care about first entries */
    /* (they will be paired with zero tap in interpolating filter)  */
    /* ------------------------------------------------------------ */

    psw1[0] = 0;
    psw2[0] = 0;
    psw1 = &psw1[1];
    psw2 = &psw2[1];
  }

  for (siI = siLoLag; siI <= siHiLag; siI++)
  {
    psw1[siI - siLoLag] = pswCSfrm[siI - LSMIN];
    psw2[siI - siLoLag] = pswGSfrm[siI - LSMIN];
  }

  if (siRemLo == 0)
  {
    /* make siLoLag correspond to first entry in temp buffers */
    /* ------------------------------------------------------ */
    siLoLag = sub(siLoLag, 1);
  }

  /* interpolate to get C's and G's which correspond to lags in list */
  /* --------------------------------------------------------------- */

  CGInterp(pswLagList, siNumLags, pswCBuf, pswGBuf, siLoLag,
           pswCInterp, pswGInterp);

  /* find max C*C*sgn(C)/G */
  /* --------------------- */

  swPeak = maxCCOverGWithSign(pswCInterp, pswGInterp, &swCmaxSqr, &swGmax,
                              siNumLags);

  /* store best lag and corresponding C*C and G */
  /* ------------------------------------------ */

  pswLTraj[siSfrmIndex] = pswLagList[swPeak];
  pswCCTraj[siSfrmIndex] = swCmaxSqr;
  pswGTraj[siSfrmIndex] = swGmax;

}


/***************************************************************************
 *
 *   FUNCTION NAME: CGInterp
 *
 *   PURPOSE:
 *
 *     Given a list of fractional lags, a C(k) array, and a G(k) array
 *     (k integer), this function generates arrays of C's and G's
 *     corresponding to the list of fractional lags by interpolating the
 *     integer C(k) and G(k) arrays.
 *
 *   INPUTS:
 *
 *     pswLIn[0:siNum-1]
 *
 *                     List of valid lags
 *
 *     siNum
 *
 *                     Length of output lists
 *
 *     pswCIn[0:variable]
 *
 *                     C(k) sequence, k integer.  The zero index corresponds
 *                     to k = siLoIntLag.
 *
 *     pswGIn[0:variable]
 *
 *                     G(k) sequence, k integer.  The zero index corresponds
 *                     to k = siLoIntLag.
 *
 *     siLoIntLag
 *
 *                     Integer lag corresponding to the first entry in the
 *                     C(k) and G(k) input arrays.
 *
 *     ppsrCGIntFilt[0:5][0:5]
 *
 *                     The FIR interpolation filter for C's and G's.
 *
 *   OUTPUTS:
 *
 *     pswCOut[0:siNum-1]
 *
 *                     List of interpolated C's corresponding to pswLIn.
 *
 *     pswGOut[0:siNum-1]
 *
 *                     List of interpolated G's corresponding to pswLIn
 *
 *   RETURN VALUE: none
 *
 *   DESCRIPTION:
 *
 *
 *   REFERENCE:  Sub-clause 4.1.8.2, 4.1.8.3 of GSM Recommendation 06.20
 *
 *   KEYWORDS: lag, interpolateCG
 *
 *************************************************************************/

void   CGInterp(Shortword pswLIn[],
                       short siNum,
                       Shortword pswCIn[],
                       Shortword pswGIn[],
                       short siLoIntLag,
                       Shortword pswCOut[],
                       Shortword pswGOut[])
{

/*_________________________________________________________________________
 |                                                                         |
 |                            Automatic Variables                          |
 |_________________________________________________________________________|
*/
  Shortword i,
         swBig,
         swLoIntLag;
  Shortword swLagInt,
         swTempRem,
         swLagRem;
  Longword L_Temp,
         L_Temp1,
         L_Temp2;

/*_________________________________________________________________________
 |                                                                         |
 |                              Executable Code                            |
 |_________________________________________________________________________|
*/

  swLoIntLag = add(siLoIntLag, (CG_INT_MACS / 2) - 1);

  for (swBig = 0; swBig < siNum; swBig++)
  {

    /* Separate integer and fractional portions of lag */
    /*-------------------------------------------------*/
    L_Temp = L_mult(pswLIn[swBig], INV_OS_FCTR);
    swLagInt = extract_h(L_Temp);

    /* swLagRem = (OS_FCTR - pswLIn[iBig] % OS_FCTR)) */
    /*---------------------------------------------------*/
    swTempRem = extract_l(L_Temp);
    swTempRem = shr(swTempRem, 1);
    swLagRem = swTempRem & SW_MAX;
    swLagRem = mult_r(swLagRem, OS_FCTR);
    swLagRem = sub(OS_FCTR, swLagRem);

    /* Get interpolated C and G values */
    /*--------------------------*/

    L_Temp1 = L_mac(32768, pswCIn[swLagInt - swLoIntLag],
                    ppsrCGIntFilt[0][swLagRem]);
    L_Temp2 = L_mac(32768, pswGIn[swLagInt - swLoIntLag],
                    ppsrCGIntFilt[0][swLagRem]);

    for (i = 1; i <= CG_INT_MACS - 1; i++)
    {
      L_Temp1 = L_mac(L_Temp1, pswCIn[i + swLagInt - swLoIntLag],
                      ppsrCGIntFilt[i][swLagRem]);
      L_Temp2 = L_mac(L_Temp2, pswGIn[i + swLagInt - swLoIntLag],
                      ppsrCGIntFilt[i][swLagRem]);

    }
    pswCOut[swBig] = extract_h(L_Temp1);
    pswGOut[swBig] = extract_h(L_Temp2);
  }
}

/***************************************************************************
 *
 *   FUNCTION NAME: CGInterpValid
 *
 *   PURPOSE:
 *
 *     The purpose of this function is to retrieve the valid (codeable) lags
 *     within one (exclusive) integer sample of the given integer lag, and
 *     interpolate the corresponding C's and G's from the integer arrays
 *
 *   INPUTS:
 *
 *     swFullResLag
 *
 *                     integer lag * OS_FCTR
 *
 *     pswCIn[0:127]
 *
 *                     integer C sequence
 *
 *     pswGIn[0:127]
 *
 *                     integer G sequence
 *
 *     psrLagTbl[0:255]
 *
 *                     reference table of valid (codeable) lags
 *
 *
 *   OUTPUTS:
 *
 *     pswLOut[0:*psiNum-1]
 *
 *                     list of valid lags within 1 of swFullResLag
 *
 *     pswCOut[0:*psiNum-1]
 *
 *                     list of interpolated C's corresponding to pswLOut
 *
 *     pswGOut[0:*psiNum-1]
 *
 *                     list of interpolated G's corresponding to pswLOut
 *
 *   RETURN VALUE:
 *
 *     siNum
 *
 *                     length of output lists
 *
 *   DESCRIPTION:
 *
 *   REFERENCE:  Sub-clause 4.1.8.2, 4.1.9 of GSM Recommendation 06.20
 *
 *   KEYWORDS: CGInterpValid, cginterpvalid, CG_INT_VALID
 *
 *************************************************************************/

short  CGInterpValid(Shortword swFullResLag,
                            Shortword pswCIn[],
                            Shortword pswGIn[],
                            Shortword pswLOut[],
                            Shortword pswCOut[],
                            Shortword pswGOut[])
{

/*_________________________________________________________________________
 |                                                                         |
 |                            Automatic Variables                          |
 |_________________________________________________________________________|
*/

  short int siLowerBound,
         siUpperBound,
         siNum,
         siI;
  Shortword swLag;

/*_________________________________________________________________________
 |                                                                         |
 |                              Executable Code                            |
 |_________________________________________________________________________|
*/

  /* Get lower and upper bounds for valid lags     */
  /* within 1 (exclusive) integer lag of input lag */
  /* --------------------------------------------- */

  swLag = sub(swFullResLag, OS_FCTR);
  swLag = quantLag(swLag, &siLowerBound);
  if (sub(swLag, swFullResLag) != 0)
  {
    siLowerBound = add(siLowerBound, 1);
  }

  swLag = add(swFullResLag, OS_FCTR);
  swLag = quantLag(swLag, &siUpperBound);
  if (sub(swLag, swFullResLag) != 0)
  {
    siUpperBound = sub(siUpperBound, 1);
  }

  /* Get list of full resolution lags whose */
  /* C's and G's will be interpolated       */
  /* -------------------------------------- */

  siNum = sub(siUpperBound, siLowerBound);
  siNum = add(siNum, 1);

  for (siI = 0; siI < siNum; siI++)
  {
    pswLOut[siI] = psrLagTbl[siI + siLowerBound];
  }

  /* Interpolate C's and G's */
  /* ----------------------- */

  CGInterp(pswLOut, siNum, pswCIn, pswGIn, LSMIN, pswCOut,
           pswGOut);

  /* Return the length of the output lists */
  /* ------------------------------------- */

  return (siNum);
}

/***************************************************************************
 *
 *   FUNCTION NAME: compResidEnergy
 *
 *   PURPOSE:
 *
 *     Computes and compares the residual energy from interpolated and
 *     non-interpolated coefficients. From the difference determines
 *     the soft interpolation decision.
 *
 *   INPUTS:
 *
 *     pswSpeech[0:159] ( [0:F_LEN-1] )
 *
 *                     Input speech frame (after high-pass filtering).
 *
 *     ppswInterpCoef[0:3][0:9] ( [0:N_SUB-1][0:NP-1] )
 *
 *                     Set of interpolated LPC direct-form coefficients for
 *                     each subframe.
 *
 *     pswPreviousCoef[0:9} ( [0:NP-1] )
 *
 *                     Set of LPC direct-form coefficients corresponding to
 *                     the previous frame
 *
 *     pswCurrentCoef[0:9} ( [0:NP-1] )
 *
 *                     Set of LPC direct-form coefficients corresponding to
 *                     the current frame
 *
 *     psnsSqrtRs[0:3] ( [0:N_SUB-1] )
 *
 *                     Array of residual energy estimates for each subframe
 *                     based on interpolated coefficients.  Used for scaling.
 *
 *   RETURN:
 *
 *     Returned value indicates the coefficients to use for each subframe:
 *     One indicates interpolated coefficients are to be used, zero indicates
 *     un-interpolated coefficients are to be used.
 *
 *   DESCRIPTION:
 *
 *
 *   REFERENCE:  Sub-clause 4.1.6 of GSM Recommendation 06.20
 *
 *   Keywords: openlooplagsearch, openloop, lag, pitch
 *
 **************************************************************************/



short  compResidEnergy(Shortword pswSpeech[],
                              Shortword ppswInterpCoef[][NP],
                              Shortword pswPreviousCoef[],
                              Shortword pswCurrentCoef[],
                              struct NormSw psnsSqrtRs[])
{

/*_________________________________________________________________________
 |                                                                         |
 |                            Automatic Variables                          |
 |_________________________________________________________________________|
*/

  short  i,
         j,
         siOverflowPossible,
         siInterpDecision;
  Shortword swMinShift,
         swShiftFactor,
         swSample,
        *pswCoef;
  Shortword pswTempState[NP];
  Shortword pswResidual[S_LEN];
  Longword L_ResidualEng;

/*_________________________________________________________________________
 |                                                                         |
 |                            Executable Code                              |
 |_________________________________________________________________________|
*/

  /* Find minimum shift count of the square-root of residual energy */
  /* estimates over the four subframes.  According to this minimum, */
  /* find a shift count for the residual signal which will be used  */
  /* to avoid overflow when the actual residual energies are        */
  /* calculated over the frame                                      */
  /*----------------------------------------------------------------*/

  swMinShift = SW_MAX;
  for (i = 0; i < N_SUB; i++)
  {

    if (sub(psnsSqrtRs[i].sh, swMinShift) < 0 && psnsSqrtRs[i].man > 0)
      swMinShift = psnsSqrtRs[i].sh;
  }

  if (sub(swMinShift, 1) >= 0)
  {

    siOverflowPossible = 0;
  }

  else if (swMinShift == 0)
  {
    siOverflowPossible = 1;
    swShiftFactor = ONE_HALF;
  }

  else if (sub(swMinShift, -1) == 0)
  {
    siOverflowPossible = 1;
    swShiftFactor = ONE_QUARTER;
  }

  else
  {
    siOverflowPossible = 1;
    swShiftFactor = ONE_EIGHTH;
  }

  /* Copy analysis filter state into temporary buffer */
  /*--------------------------------------------------*/

  for (i = 0; i < NP; i++)
    pswTempState[i] = pswAnalysisState[i];

  /* Send the speech frame, one subframe at a time, through the analysis */
  /* filter which is based on interpolated coefficients.  After each     */
  /* subframe, accumulate the energy in the residual signal, scaling to  */
  /* avoid overflow if necessary.                                        */
  /*---------------------------------------------------------------------*/

  L_ResidualEng = 0;

  for (i = 0; i < N_SUB; i++)
  {

    lpcFir(&pswSpeech[i * S_LEN], ppswInterpCoef[i], pswTempState,
           pswResidual);

    if (siOverflowPossible)
    {

      for (j = 0; j < S_LEN; j++)
      {

        swSample = mult_r(swShiftFactor, pswResidual[j]);
        L_ResidualEng = L_mac(L_ResidualEng, swSample, swSample);
      }
    }

    else
    {

      for (j = 0; j < S_LEN; j++)
      {

        L_ResidualEng = L_mac(L_ResidualEng, pswResidual[j], pswResidual[j]);
      }
    }
  }

  /* Send the speech frame, one subframe at a time, through the analysis */
  /* filter which is based on un-interpolated coefficients.  After each  */
  /* subframe, subtract the energy in the residual signal from the       */
  /* accumulated residual energy due to the interpolated coefficient     */
  /* analysis filter, again scaling to avoid overflow if necessary.      */
  /* Note that the analysis filter state is updated during these         */
  /* filtering operations.                                               */
  /*---------------------------------------------------------------------*/

  for (i = 0; i < N_SUB; i++)
  {

    switch (i)
    {

      case 0:

        pswCoef = pswPreviousCoef;
        break;

      case 1:
      case 2:
      case 3:

        pswCoef = pswCurrentCoef;
        break;
    }

    lpcFir(&pswSpeech[i * S_LEN], pswCoef, pswAnalysisState,
           pswResidual);

    if (siOverflowPossible)
    {

      for (j = 0; j < S_LEN; j++)
      {

        swSample = mult_r(swShiftFactor, pswResidual[j]);
        L_ResidualEng = L_msu(L_ResidualEng, swSample, swSample);
      }
    }

    else
    {

      for (j = 0; j < S_LEN; j++)
      {

        L_ResidualEng = L_msu(L_ResidualEng, pswResidual[j], pswResidual[j]);
      }
    }
  }

  /* Make soft-interpolation decision based on the difference in residual */
  /* energies                                                             */
  /*----------------------------------------------------------------------*/

  if (L_ResidualEng < 0)
    siInterpDecision = 1;

  else
    siInterpDecision = 0;

  return siInterpDecision;
}

/***************************************************************************
 *
 *    FUNCTION NAME: cov32
 *
 *    PURPOSE: Calculates B, F, and C correlation matrices from which
 *             the reflection coefficients are computed using the FLAT
 *             algorithm. The Spectral Smoothing Technique (SST) is applied
 *             to the correlations. End point correction is employed
 *             in computing the correlations to minimize computation.
 *
 *    INPUT:
 *
 *       pswIn[0:169]
 *                     A sampled speech vector used to compute
 *                     correlations need for generating the optimal
 *                     reflection coefficients via the FLAT algorithm.
 *
 *       CVSHIFT       The number of right shifts by which the normalized
 *                     correlations are to be shifted down prior to being
 *                     rounded into the Shortword output correlation arrays
 *                     B, F, and C.
 *
 *       pL_rFlatSstCoefs[NP]
 *
 *                     A table stored in Rom containing the spectral
 *                     smoothing function coefficients.
 *
 *    OUTPUTS:
 *
 *       pppL_B[0:NP-1][0:NP-1][0:1]
 *                     An output correlation array containing the backward
 *                     correlations of the input signal. It is a square
 *                     matrix symmetric about the diagonal. Only the upper
 *                     right hand triangular region of this matrix is
 *                     initialized, but two dimensional indexing is retained
 *                     to enhance clarity. The third array dimension is used
 *                     by function flat to swap the current and the past
 *                     values of B array, eliminating the need to copy
 *                     the updated B values onto the old B values at the
 *                     end of a given lattice stage. The third dimension
 *                     is similarily employed in arrays F and C.
 *
 *       pppL_F[0:NP-1][0:NP-1][0:1]
 *                     An output correlation array containing the forward
 *                     correlations of the input signal. It is a square
 *                     matrix symmetric about the diagonal. Only the upper
 *                     right hand triangular region of this matrix is
 *                     initialized.
 *
 *       pppL_C[0:NP-1][0:NP-1][0:1]
 *                     An output correlation array containing the cross
 *                     correlations of the input signal. It is a square
 *                     matrix which is not symmetric. All its elements
 *                     are initialized, for the third dimension index = 0.
 *
 *       pL_R0         Average normalized signal power over F_LEN
 *                     samples, given by 0.5*(Phi(0,0)+Phi(NP,NP)), where
 *                     Phi(0,0) and Phi(NP,NP) are normalized signal
 *                     autocorrelations.  The average unnormalized signal
 *                     power over the frame is given by adjusting L_R0 by
 *                     the shift count which is returned. pL_R0 along
 *                     with the returned shift count are the inputs to
 *                     the frame energy quantizer.
 *
 *        Longword pL_VadAcf[4]
 *                     An array with the autocorrelation coefficients to be
 *                     used by the VAD.
 *
 *        Shortword *pswVadScalAuto
 *                     Input scaling factor used by the VAD.
 *
 *    RETURN:
 *
 *       swNormPwr     The shift count to be applied to pL_R0 for
 *                     reconstructing the average unnormalized
 *                     signal power over the frame.
 *                     Negative shift count means that a left shift was
 *                     applied to the correlations to achieve a normalized
 *                     value of pL_R0.
 *
 *   DESCRIPTION:
 *
 *
 *      The input energy of the signal is assumed unknown.  It maximum
 *      can be F_LEN*0.5. The 0.5 factor accounts for scaling down of the
 *      input signal in the high-pass filter.  Therefore the signal is
 *      shifted down by 3 shifts producing an energy reduction of 2^(2*3)=64.
 *      The resulting energy is then normalized.  Based on the shift count,
 *      the correlations F, B, and C are computed using as few shifts as
 *      possible, so a precise result is attained.
 *      This is an implementation of equations: 2.1 through 2.11.
 *
 *   REFERENCE:  Sub-clause 4.1.3 of GSM Recommendation 06.20
 *
 *   keywords: energy, autocorrelation, correlation, cross-correlation
 *   keywords: spectral smoothing, SST, LPC, FLAT, flat
 *
 *************************************************************************/

Shortword cov32(Shortword pswIn[],
                       Longword pppL_B[NP][NP][2],
                       Longword pppL_F[NP][NP][2],
                       Longword pppL_C[NP][NP][2],
                       Longword *pL_R0,
                       Longword pL_VadAcf[],
                       Shortword *pswVadScalAuto)
{

/*_________________________________________________________________________
 |                                                                         |
 |                            Automatic Variables                          |
 |_________________________________________________________________________|
*/

  Longword L_max,
         L_Pwr0,
         L_Pwr,
         L_temp,
         pL_Phi[NP + 1];
  Shortword swTemp,
         swNorm,
         swNormSig,
         swNormPwr,
         pswInScale[A_LEN],
         swPhiNorm;
  short int i,
         k,
         kk,
         n;

/*_________________________________________________________________________
 |                                                                         |
 |                              Executable Code                            |
 |_________________________________________________________________________|
*/

  /* Calculate energy in the frame vector (160 samples) for each   */
  /* of NP frame placements. The energy is reduced by 64. This is  */
  /* accomplished by shifting the input right by 3 bits. An offset */
  /* of 0x117f0b is placed into the accumulator to account for     */
  /* the worst case power gain due to the 3 LSB's of the input     */
  /* signal which were right shifted. The worst case is that the   */
  /* 3 LSB's were all set to 1 for each of the samples. Scaling of */
  /* the input by a half is assumed here.                          */
  /*---------------------------------------------------------------*/

  L_max = 0;
  for (L_Pwr = 0x117f0b, i = 0; i < F_LEN; i++)
  {
    swTemp = shr(pswIn[i], 3);
    L_Pwr = L_mac(L_Pwr, swTemp, swTemp);
  }
  L_max |= L_Pwr;

  /* L_max tracks the maximum power over NP window placements */
  /*----------------------------------------------------------*/

  for (i = 1; i <= NP; i++)
  {

    /* Subtract the power due to 1-st sample from previous window
     * placement. */
    /*-----------------------------------------------------------*/

    swTemp = shr(pswIn[i - 1], 3);
    L_Pwr = L_msu(L_Pwr, swTemp, swTemp);

    /* Add the power due to new sample at the current window placement. */
    /*------------------------------------------------------------------*/

    swTemp = shr(pswIn[F_LEN + i - 1], 3);
    L_Pwr = L_mac(L_Pwr, swTemp, swTemp);

    L_max |= L_Pwr;

  }

  /* Compute the shift count needed to achieve normalized value */
  /* of the correlations.                                       */
  /*------------------------------------------------------------*/

  swTemp = norm_l(L_max);
  swNorm = sub(6, swTemp);

  if (swNorm >= 0)
  {

    /* The input signal needs to be shifted down, to avoid limiting */
    /* so compute the shift count to be applied to the input.       */
    /*--------------------------------------------------------------*/

    swTemp = add(swNorm, 1);
    swNormSig = shr(swTemp, 1);
    swNormSig = add(swNormSig, 0x0001);

  }
  else
  {
    /* No scaling down of the input is necessary */
    /*-------------------------------------------*/

    swNormSig = 0;

  }

  /* Convert the scaling down, if any, which was done to the time signal */
  /* to the power domain, and save.                                      */
  /*---------------------------------------------------------------------*/

  swNormPwr = shl(swNormSig, 1);

  /* Buffer the input signal, scaling it down if needed. */
  /*-----------------------------------------------------*/

  for (i = 0; i < A_LEN; i++)
  {
    pswInScale[i] = shr(pswIn[i], swNormSig);
  }

  /* Compute from buffered (scaled) input signal the correlations     */
  /* needed for the computing the reflection coefficients.            */
  /*------------------------------------------------------------------*/

  /* Compute correlation Phi(0,0) */
  /*------------------------------*/

  L_Pwr = L_mult(pswInScale[NP], pswInScale[NP]);
  for (n = 1; n < F_LEN; n++)
  {
    L_Pwr = L_mac(L_Pwr, pswInScale[NP + n], pswInScale[NP + n]);
  }
  pL_Phi[0] = L_Pwr;

  /* Get ACF[0] and input scaling factor for VAD algorithm */
  *pswVadScalAuto = swNormSig;
  pL_VadAcf[0] = L_Pwr;

  /* Compute the remaining correlations along the diagonal which */
  /* starts at Phi(0,0). End-point correction is employed to     */
  /* limit computation.                                          */
  /*-------------------------------------------------------------*/

  for (i = 1; i <= NP; i++)
  {

    /* Compute the power in the last sample from the previous         */
    /* window placement, and subtract it from correlation accumulated */
    /* at the previous window placement.                              */
    /*----------------------------------------------------------------*/

    L_Pwr = L_msu(L_Pwr, pswInScale[NP + F_LEN - i],
                  pswInScale[NP + F_LEN - i]);

    /* Compute the power in the new sample for the current window       */
    /* placement, and add it to L_Pwr to obtain the value of Phi(i,i). */
    /*------------------------------------------------------------------*/

    L_Pwr = L_mac(L_Pwr, pswInScale[NP - i], pswInScale[NP - i]);

    pL_Phi[i] = L_Pwr;

  }

  /* Compute the shift count necessary to normalize the Phi array  */
  /*---------------------------------------------------------------*/

  L_max = 0;
  for (i = 0; i <= NP; i++)
  {
    L_max |= pL_Phi[i];
  }
  swPhiNorm = norm_l(L_max);

  /* Adjust the shift count to be returned to account for any scaling */
  /* down which might have been done to the input signal prior to     */
  /* computing the correlations.                                      */
  /*------------------------------------------------------------------*/

  swNormPwr = sub(swNormPwr, swPhiNorm);

  /* Compute the average power over the frame; i.e.,                   */
  /* 0.5*(Phi(0,0)+Phi(NP,NP)), given a normalized pL_Phi array.       */
  /*-------------------------------------------------------------------*/

  swTemp = sub(swPhiNorm, 1);
  L_Pwr0 = L_shl(pL_Phi[0], swTemp);
  L_Pwr = L_shl(pL_Phi[NP], swTemp);
  *pL_R0 = L_add(L_Pwr, L_Pwr0);       /* Copy power to output pointer */

  /* Check if the average power is normalized; if not, shift left by 1 bit */
  /*-----------------------------------------------------------------------*/

  if (!(*pL_R0 & 0x40000000))
  {
    *pL_R0 = L_shl(*pL_R0, 1);         /* normalize the average power    */
    swNormPwr = sub(swNormPwr, 1);     /* adjust the shift count         */
  }

  /* Reduce the shift count needed to normalize the correlations   */
  /* by CVSHIFT bits.                                              */
  /*---------------------------------------------------------------*/

  swNorm = sub(swPhiNorm, CVSHIFT);

  /* Initialize the F, B, and C output correlation arrays, using the */
  /* Phi correlations computed along the diagonal of symmetry.       */
  /*-----------------------------------------------------------------*/

  L_temp = L_shl(pL_Phi[0], swNorm);   /* Normalize the result     */

  pppL_F[0][0][0] = L_temp;            /* Write to output array    */

  for (i = 1; i <= NP - 1; i++)
  {

    L_temp = L_shl(pL_Phi[i], swNorm); /* Normalize the result     */


    pppL_F[i][i][0] = L_temp;          /* Write to output array    */
    pppL_B[i - 1][i - 1][0] = L_temp;  /* Write to output array    */
    pppL_C[i][i - 1][0] = L_temp;      /* Write to output array    */

  }

  L_temp = L_shl(pL_Phi[NP], swNorm);  /* Normalize the result     */

  pppL_B[NP - 1][NP - 1][0] = L_temp;  /* Write to output array    */

  for (k = 1; k <= NP - 1; k++)
  {

    /* Compute correlation Phi(0,k) */
    /*------------------------------*/

    L_Pwr = L_mult(pswInScale[NP], pswInScale[NP - k]);
    for (n = 1; n < F_LEN; n++)
    {
      L_Pwr = L_mac(L_Pwr, pswInScale[NP + n], pswInScale[NP + n - k]);
    }
    /* convert covariance values to ACF and store for VAD algorithm */
    if (k < 9)
    {
      pL_VadAcf[k] = L_Pwr;
      for (kk = 0; kk < k; kk++)
      {
        pL_VadAcf[k] = L_msu(pL_VadAcf[k], pswInScale[NP + kk],
                             pswInScale[NP + kk - k]);
      }
    }

    L_temp = L_shl(L_Pwr, swNorm);     /* Normalize the result */
    L_temp = L_mpy_ll(L_temp, pL_rFlatSstCoefs[k - 1]); /* Apply SST */

    pppL_F[0][k][0] = L_temp;          /* Write to output array    */
    pppL_C[0][k - 1][0] = L_temp;      /* Write to output array    */


    /* Compute the remaining correlations along the diagonal which */
    /* starts at Phi(0,k). End-point correction is employed to     */
    /* limit computation.                                          */
    /*-------------------------------------------------------------*/

    for (kk = k + 1, i = 1; kk <= NP - 1; kk++, i++)
    {

      /* Compute the power in the last sample from the previous         */
      /* window placement, and subtract it from correlation accumulated */
      /* at the previous window placement.                              */
      /*----------------------------------------------------------------*/

      L_Pwr = L_msu(L_Pwr, pswInScale[NP + F_LEN - i],
                    pswInScale[NP + F_LEN - kk]);

      /* Compute the power in the new sample for the current window       */
      /* placement, and add it to L_Pwr to obtain the value of Phi(i,kk). */
      /*------------------------------------------------------------------*/

      L_Pwr = L_mac(L_Pwr, pswInScale[NP - i], pswInScale[NP - kk]);

      L_temp = L_shl(L_Pwr, swNorm);   /* Normalize */
      L_temp = L_mpy_ll(L_temp, pL_rFlatSstCoefs[k - 1]);     /* Apply SST */

      pppL_F[i][kk][0] = L_temp;       /* Write to output array */
      pppL_B[i - 1][kk - 1][0] = L_temp;        /* Write to output array */
      pppL_C[i][kk - 1][0] = L_temp;   /* Write to output array    */
      pppL_C[kk][i - 1][0] = L_temp;   /* Write to output array    */

    }

    /* Compute the power in the last sample from the previous         */
    /* window placement, and subtract it from correlation accumulated */
    /* at the previous window placement.                              */
    /*----------------------------------------------------------------*/

    L_Pwr = L_msu(L_Pwr, pswInScale[F_LEN + k], pswInScale[F_LEN]);

    /* Compute the power in the new sample for the current window       */
    /* placement, and add it to L_Pwr to obtain the value of Phi(i,kk). */
    /*------------------------------------------------------------------*/

    L_Pwr = L_mac(L_Pwr, pswInScale[k], pswInScale[0]);

    L_temp = L_shl(L_Pwr, swNorm);     /* Normalize the result */
    L_temp = L_mpy_ll(L_temp, pL_rFlatSstCoefs[k - 1]); /* Apply SST */

    pppL_B[NP - k - 1][NP - 1][0] = L_temp;     /* Write to output array */
    pppL_C[NP - k][NP - 1][0] = L_temp;/* Write to output array */

  }

  /* Compute correlation Phi(0,NP) */
  /*-------------------------------*/

  L_Pwr = L_mult(pswInScale[NP], pswInScale[0]);
  for (n = 1; n < F_LEN; n++)
  {
    L_Pwr = L_mac(L_Pwr, pswInScale[NP + n], pswInScale[n]);
  }

  L_temp = L_shl(L_Pwr, swNorm);       /* Normalize the result */
  L_temp = L_mpy_ll(L_temp, pL_rFlatSstCoefs[NP - 1]);  /* Apply SST */

  pppL_C[0][NP - 1][0] = L_temp;       /* Write to output array */

  return (swNormPwr);

}

/***************************************************************************
 *
 *    FUNCTION NAME: filt4_2nd
 *
 *    PURPOSE:  Implements a fourth order filter by cascading two second
 *              order sections.
 *
 *    INPUTS:
 *
 *      pswCoef[0:9]   An array of two sets of filter coefficients.
 *
 *      pswIn[0:159]   An array of input samples to be filtered, filtered
 *                     output samples written to the same array.
 *
 *      pswXstate[0:3] An array containing x-state memory for two 2nd order
 *                     filter sections.
 *
 *      pswYstate[0:7] An array containing y-state memory for two 2nd order
 *                     filter sections.
 *
 *      npts           Number of samples to filter (must be even).
 *
 *      shifts         number of shifts to be made on output y(n).
 *
 *    OUTPUTS:
 *
 *       pswIn[0:159]  Output array containing the filtered input samples.
 *
 *    RETURN:
 *
 *       none.
 *
 *    DESCRIPTION:
 *
 *    data structure:
 *
 *    Coeff array order:  (b2,b1,b0,a2,a1)Section 1;(b2,b1,b0,a2,a1)Section 2
 *    Xstate array order: (x(n-2),x(n-1))Section 1; (x(n-2),x(n-1))Section 2
 *    Ystate array order: y(n-2)MSB,y(n-2)LSB,y(n-1)MSB,y(n-1)LSB Section 1
 *                        y(n-2)MSB,y(n-2)LSB,y(n-1)MSB,y(n-1)LSB Section 2
 *
 *    REFERENCE:  Sub-clause 4.1.1 GSM Recommendation 06.20
 *
 *    KEYWORDS: highpass filter, hp, HP, filter
 *
 *************************************************************************/

void   filt4_2nd(Shortword pswCoeff[], Shortword pswIn[],
                        Shortword pswXstate[], Shortword pswYstate[],
                        int npts, int shifts)
{

/*_________________________________________________________________________
 |                                                                         |
 |                              Executable Code                            |
 |_________________________________________________________________________|
*/

  /* Do first second order section */
  /*-------------------------------*/

  iir_d(&pswCoeff[0],pswIn,&pswXstate[0],&pswYstate[0],npts,shifts,1,0);


  /* Do second second order section */
  /*--------------------------------*/

  iir_d(&pswCoeff[5],pswIn,&pswXstate[2],&pswYstate[4],npts,shifts,0,1);

}

/***************************************************************************
 *
 *   FUNCTION NAME: findBestInQuantList
 *
 *   PURPOSE:
 *     Given a list of quantizer vectors and their associated prediction
 *     errors, search the list for the iNumVectOut vectors and output them
 *     as a new list.
 *
 *   INPUTS: psqlInList, iNumVectOut
 *
 *   OUTPUTS: psqlBestOutList
 *
 *   RETURN VALUE: none
 *
 *   DESCRIPTION:
 *
 *     The AFLAT recursion yields prediction errors.  This routine finds
 *     the lowest candidate is the AFLAT recursion outputs.
 *
 *
 *   KEYWORDS: best quantlist find
 *
 *   REFERENCE:  Sub-clause 4.1.4.1 GSM Recommendation 06.20
 *
 *************************************************************************/

void findBestInQuantList(struct QuantList psqlInList,
                                int iNumVectOut,
                                struct QuantList psqlBestOutList[])
{

/*_________________________________________________________________________
 |                                                                         |
 |                            Automatic Variables                          |
 |_________________________________________________________________________|
*/
  int    quantIndex,
         bstIndex,
         i;

/*_________________________________________________________________________
 |                                                                         |
 |                              Executable Code                            |
 |_________________________________________________________________________|
*/

  /* initialize the best list */
  /* invalidate, ensure they will be dropped */
  for (bstIndex = 0; bstIndex < iNumVectOut; bstIndex++)
  {
    psqlBestOutList[bstIndex].iNum = 1;
    psqlBestOutList[bstIndex].iRCIndex = psqlInList.iRCIndex;
    psqlBestOutList[bstIndex].pswPredErr[0] = 0x7fff;
  }

  /* best list elements replaced in the order:  0,1,2,3... challenger must
   * be < (not <= ) current best */
  for (quantIndex = 0; quantIndex < psqlInList.iNum; quantIndex++)
  {
    bstIndex = 0;
    while (sub(psqlInList.pswPredErr[quantIndex],
               psqlBestOutList[bstIndex].pswPredErr[0]) >= 0 &&
           bstIndex < iNumVectOut)
    {
      bstIndex++;                      /* only increments to next upon
                                        * failure to beat "best" */
    }

    if (bstIndex < iNumVectOut)
    {                                  /* a new value is found */
      /* now add challenger to best list at index bstIndex */
      for (i = iNumVectOut - 1; i > bstIndex; i--)
      {
        psqlBestOutList[i].pswPredErr[0] =
                psqlBestOutList[i - 1].pswPredErr[0];
        psqlBestOutList[i].iRCIndex =
                psqlBestOutList[i - 1].iRCIndex;
      }
      /* get new best value and place in list */
      psqlBestOutList[bstIndex].pswPredErr[0] =
              psqlInList.pswPredErr[quantIndex];
      psqlBestOutList[bstIndex].iRCIndex =
              psqlInList.iRCIndex + quantIndex;
    }
  }
}

/***************************************************************************
 *
 *   FUNCTION NAME: findPeak
 *
 *   PURPOSE:
 *
 *     The purpose of this function is to return the lag
 *     that maximizes CC/G within +- PEAK_VICINITY of the
 *     input lag.  The input lag is an integer lag, and
 *     the search for a peak is done on the surrounding
 *     integer lags.
 *
 *   INPUTS:
 *
 *     swSingleResLag
 *
 *                     Input integer lag, expressed as lag * OS_FCTR
 *
 *     pswCIn[0:127]
 *
 *                     C(k) sequence, k an integer
 *
 *     pswGIn[0:127]
 *
 *                     G(k) sequence, k an integer
 *
 *   OUTPUTS:
 *
 *     none
 *
 *   RETURN VALUE:
 *
 *     Integer lag where peak was found, or zero if no peak was found.
 *     The lag is expressed as lag * OS_FCTR
 *
 *   DESCRIPTION:
 *
 *     This routine is called from pitchLags(), and is used to do the
 *     interpolating CC/G peak search.  This is used in a number of
 *     places in pitchLags().  See description 5.3.1.
 *
 *   REFERENCE:  Sub-clause 4.1.8.2 of GSM Recommendation 06.20
 *
 *   KEYWORDS:
 *
 *************************************************************************/

Shortword findPeak(Shortword swSingleResLag,
                          Shortword pswCIn[],
                          Shortword pswGIn[])
{

/*_________________________________________________________________________
 |                                                                         |
 |                            Automatic Variables                          |
 |_________________________________________________________________________|
*/

  Shortword swCmaxSqr,
         swGmax,
         swFullResPeak;
  short int siUpperBound,
         siLowerBound,
         siRange,
         siPeak;

/*_________________________________________________________________________
 |                                                                         |
 |                              Executable Code                            |
 |_________________________________________________________________________|
*/

  /* get upper and lower bounds for integer lags for peak search */
  /* ----------------------------------------------------------- */

  siUpperBound = add(swSingleResLag, PEAK_VICINITY + 1);
  if (sub(siUpperBound, LMAX + 1) > 0)
  {
    siUpperBound = LMAX + 1;
  }

  siLowerBound = sub(swSingleResLag, PEAK_VICINITY + 1);
  if (sub(siLowerBound, LMIN - 1) < 0)
  {
    siLowerBound = LMIN - 1;
  }

  siRange = sub(siUpperBound, siLowerBound);
  siRange = add(siRange, 1);

  /* do peak search */
  /* -------------- */

  swCmaxSqr = 0;
  swGmax = 0x3f;

  siPeak = fnBest_CG(&pswCIn[siLowerBound - LSMIN],
                     &pswGIn[siLowerBound - LSMIN], &swCmaxSqr, &swGmax,
                     siRange);

  /* if no max found, flag no peak */
  /* ----------------------------- */

  if (add(siPeak, 1) == 0)
  {
    swFullResPeak = 0;
  }

  /* determine peak location      */
  /* if at boundary, flag no peak */
  /* else return lag at peak      */
  /* ---------------------------- */

  else
  {
    siPeak = add(siPeak, siLowerBound);

    if ((sub(siPeak, siLowerBound) == 0) ||
        (sub(siPeak, siUpperBound) == 0))
    {
      swFullResPeak = 0;
    }
    else
    {
      swFullResPeak = shr(extract_l(L_mult(siPeak, OS_FCTR)), 1);
    }
  }
  return (swFullResPeak);
}

/***************************************************************************
 *
 *    FUNCTION NAME: flat
 *
 *    PURPOSE:  Computes the unquantized reflection coefficients from the
 *              input speech using the FLAT algorithm. Also computes the
 *              frame energy, and the index of the element in the R0
 *              quantization table which best represents the frame energy.
 *              Calls function cov32 which computes the F, B, and C
 *              correlation arrays, required by the FLAT algorithm to
 *              compute the reflection coefficients.
 *
 *    INPUT:
 *
 *       pswSpeechIn[0:169]
 *                     A sampled speech vector used to compute
 *                     correlations need for generating the optimal
 *                     reflection coefficients via the FLAT algorithm.
 *
 *    OUTPUTS:
 *
 *       pswRc[NP]     An array of unquantized reflection coefficients.
 *
 *       *piR0Inx      An index of the quantized frame energy value.
 *
 *       Longword pL_VadAcf[4]
 *                     An array with the autocorrelation coefficients to be
 *                     used by the VAD.  Generated by cov16(), a daughter
 *                     function of flat().
 *
 *       Shortword *pswVadScalAuto
 *                     Input scaling factor used by the VAD.
 *                     Generated by cov16(), a daughter function of flat().
 *                     function.
 *
 *    RETURN:          L_R0 normalized frame energy value, required in DTX
 *                     mode.
 *
 *   DESCRIPTION:
 *
 *    An efficient Fixed point LAtice Technique (FLAT) is used to compute
 *    the reflection coefficients, given B, F, and C arrays returned by
 *    function cov32. B, F, and C are backward, forward, and cross
 *    correlations computed from the input speech. The correlations
 *    are spectrally smoothed in cov32.
 *
 *
 *   REFERENCE:  Sub-clause 4.1.3 of GSM Recommendation 06.20
 *
 *   keywords: LPC, FLAT, reflection coefficients, covariance, correlation,
 *   keywords: spectrum, energy, R0, spectral smoothing, SST
 *
 *************************************************************************/

Longword   flat(Shortword pswSpeechIn[],
                   Shortword pswRc[],
                   int *piR0Inx,
                   Longword pL_VadAcf[],
                   Shortword *pswVadScalAuto)
{

/*_________________________________________________________________________
 |                                                                         |
 |                            Automatic Variables                          |
 |_________________________________________________________________________|
*/
  Shortword
         swNum,
         swDen,
         swRcSq,
         swSqrtOut,
         swRShifts,
         swShift,
         swShift1;
  Longword
         pppL_F[NP][NP][2],
         pppL_B[NP][NP][2],
         pppL_C[NP][NP][2],
         L_Num,
         L_TmpA,
         L_TmpB,
         L_temp,
         L_sum,
         L_R0,
         L_Fik,
         L_Bik,
         L_Cik,
         L_Cki;
  short int i,
         j,
         k,
         l,
         j_0,
         j_1;

/*_________________________________________________________________________
 |                                                                         |
 |                              Executable Code                            |
 |_________________________________________________________________________|
*/

  /* Compute from the input speech the elements of the B, F, and C     */
  /* arrays, which form the initial conditions for the FLAT algorithm. */
  /*-------------------------------------------------------------------*/

  swRShifts = cov32(pswSpeechIn, pppL_B, pppL_F, pppL_C, &L_R0,
                    pL_VadAcf, pswVadScalAuto);

  /* Compute the intermediate quantities required by the R0 quantizer */
  /*------------------------------------------------------------------*/

  if (L_R0 != 0)
  {
    swSqrtOut = sqroot(L_R0);          /* If L_R0 > 0, compute sqrt */
  }
  else
  {
    swSqrtOut = 0;                     /* L_R0 = 0, initialize sqrt(0) */
  }

  swRShifts = sub(S_SH + 2, swRShifts);

  /* If odd number of shifts compensate by sqrt(0.5) */
  /*-------------------------------------------------*/

  if (swRShifts & 1)
  {
    L_temp = L_mult(swSqrtOut, 0x5a82);
  }
  else
  {
    L_temp = L_deposit_h(swSqrtOut);
  }
  swRShifts = shr(swRShifts, 1);

  if (swRShifts > 0)
  {
    L_temp = L_shr(L_temp, swRShifts);
  }
  else if (swRShifts < 0)
  {
    L_temp = 0;
  }

  /* Given average energy L_temp, find the index in the R0 quantization */
  /* table which best represents it.                                    */
  /*--------------------------------------------------------------------*/

  *piR0Inx = r0Quant(L_temp);

  L_R0 = L_temp; /* save the unquantized R0 */             /* DTX mode */

  /* Zero out the number of left-shifts to be applied to the  */
  /* F, B, and C matrices.                                    */
  /*----------------------------------------------------------*/

  swShift = 0;

  /* Now compute the NP reflection coefficients  */
  /*---------------------------------------------*/

  for (j = 0; j < NP; j++)
  {

    /* Initialize the modulo indices of the third dimension of arrays  */
    /* F, B, and C, where indices j_0 and j_1 point to:                */
    /* */
    /* j_0 = points to F, B, and C matrix values at stage j-0, which   */
    /* is the current  lattice stage.                                  */
    /* j_1 = points to F, B, and C matrix values at stage j-1, which   */
    /* is the previous lattice stage.                                  */
    /* */
    /* Use of modulo address arithmetic permits to swap values of j_0 and */
    /* and j_1 at each lattice stage, thus eliminating the need to copy   */
    /* the current elements of F, B, and C arrays, into the F, B, and C   */
    /* arrays corresponding to the previous lattice stage, prior to       */
    /* incrementing j, the index of the lattice filter stage.             */
    /*--------------------------------------------------------------------*/

    j_0 = (j + 1) % 2;
    j_1 = j % 2;

    /* Get the numerator for computing the j-th reflection coefficient */
    /*-----------------------------------------------------------------*/

    L_Num = L_add(L_shl(pppL_C[0][0][j_1], swShift),
                  L_shl(pppL_C[NP - j - 1][NP - j - 1][j_1], swShift));

    /* Get the denominator for computing the j-th reflection coefficient */
    /*-------------------------------------------------------------------*/

    L_temp = L_add(L_shl(pppL_F[0][0][j_1], swShift),
                   L_shl(pppL_B[0][0][j_1], swShift));
    L_TmpA = L_add(L_shl(pppL_F[NP - j - 1][NP - j - 1][j_1], swShift),
                   L_shl(pppL_B[NP - j - 1][NP - j - 1][j_1], swShift));
    L_sum = L_add(L_TmpA, L_temp);
    L_sum = L_shr(L_sum, 1);

    /* Normalize the numerator and the denominator terms */
    /*---------------------------------------------------*/

    swShift1 = norm_s(extract_h(L_sum));

    L_sum = L_shl(L_sum, swShift1);
    L_Num = L_shl(L_Num, swShift1);

    swNum = round(L_Num);
    swDen = round(L_sum);

    if (swDen <= 0)
    {

      /* Zero prediction error at the j-th lattice stage, zero */
      /* out remaining reflection coefficients and return.     */
      /*-------------------------------------------------------*/

      for (i = j; i < NP; i++)
      {
        pswRc[i] = 0;
      }

      return (L_R0);
    }
    else
    {

      /* Non-zero prediction error, check if the j-th reflection
       * coefficient */
      /* about to be computed is stable.                           */
      /*-----------------------------------------------------------*/

      if (sub(abs_s(swNum), swDen) >= 0)
      {

        /* Reflection coefficient at j-th lattice stage unstable, so zero  */
        /* out reflection coefficients for lattice stages i=j,...,NP-1, and */
        /* return.                                                         */
        /*-----------------------------------------------------------------*/

        for (i = j; i < NP; i++)
        {
          pswRc[i] = 0;
        }

        return (L_R0);
      }
      else
      {

        /* j-th reflection coefficient is stable, compute it. */
        /*----------------------------------------------------*/

        if (swNum < 0)
        {

          swNum = negate(swNum);
          pswRc[j] = divide_s(swNum, swDen);

        }
        else
        {

          pswRc[j] = divide_s(swNum, swDen);
          pswRc[j] = negate(pswRc[j]);

        }                              /* j-th reflection coefficient
                                        * sucessfully computed. */
        /*----------------------------------------------------*/


      }                                /* End of reflection coefficient
                                        * stability test (and computation) */
      /*------------------------------------------------------------------*/

    }                                  /* End of non-zero prediction error
                                        * case */
    /*----------------------------------------*/



    /* If not at the last lattice stage, update F, B, and C arrays */
    /*-------------------------------------------------------------*/

    if (j != NP - 1)
    {

      /* Compute squared Rc[j] */
      /*-----------------------*/

      swRcSq = mult_r(pswRc[j], pswRc[j]);

      i = 0;
      k = 0;

      /* Compute the common terms used by the FLAT recursion to reduce */
      /* computation.                                                  */
      /*---------------------------------------------------------------*/

      L_Cik = L_shl(pppL_C[i][k][j_1], swShift);

      L_TmpA = L_add(L_Cik, L_Cik);
      L_TmpA = L_mpy_ls(L_TmpA, pswRc[j]);

      /* Update the F array */
      /*--------------------*/

      L_Fik = L_shl(pppL_F[i][k][j_1], swShift);
      L_Bik = L_shl(pppL_B[i][k][j_1], swShift);

      L_temp = L_mpy_ls(L_Bik, swRcSq);
      L_temp = L_add(L_temp, L_Fik);
      pppL_F[i][k][j_0] = L_add(L_temp, L_TmpA);

      for (k = i + 1; k <= NP - j - 2; k++)
      {

        /* Compute the common terms used by the FLAT recursion to reduce */
        /* computation.                                                  */
        /*---------------------------------------------------------------*/

        L_Cik = L_shl(pppL_C[i][k][j_1], swShift);
        L_Cki = L_shl(pppL_C[k][i][j_1], swShift);

        L_TmpA = L_add(L_Cik, L_Cki);
        L_TmpA = L_mpy_ls(L_TmpA, pswRc[j]);

        L_Bik = L_shl(pppL_B[i][k][j_1], swShift);
        L_Fik = L_shl(pppL_F[i][k][j_1], swShift);

        L_TmpB = L_add(L_Bik, L_Fik);
        L_TmpB = L_mpy_ls(L_TmpB, pswRc[j]);

        /* Update the F and C arrays */
        /*---------------------------------*/

        L_temp = L_mpy_ls(L_Bik, swRcSq);
        L_temp = L_add(L_temp, L_Fik);
        pppL_F[i][k][j_0] = L_add(L_temp, L_TmpA);

        L_temp = L_mpy_ls(L_Cki, swRcSq);
        L_temp = L_add(L_temp, L_Cik);
        pppL_C[i][k - 1][j_0] = L_add(L_temp, L_TmpB);

      }

      k = NP - j - 1;

      /* Compute the common terms used by the FLAT recursion to reduce */
      /* computation.                                                  */
      /*---------------------------------------------------------------*/

      L_Bik = L_shl(pppL_B[i][k][j_1], swShift);
      L_Fik = L_shl(pppL_F[i][k][j_1], swShift);

      L_TmpB = L_add(L_Bik, L_Fik);
      L_TmpB = L_mpy_ls(L_TmpB, pswRc[j]);

      /* Update the C array */
      /*-----------------------*/

      L_Cik = L_shl(pppL_C[i][k][j_1], swShift);
      L_Cki = L_shl(pppL_C[k][i][j_1], swShift);

      L_temp = L_mpy_ls(L_Cki, swRcSq);
      L_temp = L_add(L_temp, L_Cik);
      pppL_C[i][k - 1][j_0] = L_add(L_temp, L_TmpB);


      for (i = 1; i <= NP - j - 2; i++)
      {

        k = i;

        /* Compute the common terms used by the FLAT recursion to reduce */
        /* computation.                                                  */
        /*---------------------------------------------------------------*/

        L_Cik = L_shl(pppL_C[i][k][j_1], swShift);

        L_TmpA = L_add(L_Cik, L_Cik);
        L_TmpA = L_mpy_ls(L_TmpA, pswRc[j]);

        L_Bik = L_shl(pppL_B[i][k][j_1], swShift);
        L_Fik = L_shl(pppL_F[i][k][j_1], swShift);

        L_TmpB = L_add(L_Bik, L_Fik);
        L_TmpB = L_mpy_ls(L_TmpB, pswRc[j]);

        /* Update F, B and C arrays */
        /*-----------------------------------*/

        L_temp = L_mpy_ls(L_Bik, swRcSq);
        L_temp = L_add(L_temp, L_Fik);
        pppL_F[i][k][j_0] = L_add(L_temp, L_TmpA);

        L_temp = L_mpy_ls(L_Fik, swRcSq);
        L_temp = L_add(L_temp, L_Bik);
        pppL_B[i - 1][k - 1][j_0] = L_add(L_temp, L_TmpA);

        L_temp = L_mpy_ls(L_Cik, swRcSq);
        L_temp = L_add(L_temp, L_Cik);
        pppL_C[i][k - 1][j_0] = L_add(L_temp, L_TmpB);

        for (k = i + 1; k <= NP - j - 2; k++)
        {

          /* Compute the common terms used by the FLAT recursion to reduce */
          /* computation.                                                  */
          /*---------------------------------------------------------------*/

          L_Cik = L_shl(pppL_C[i][k][j_1], swShift);
          L_Cki = L_shl(pppL_C[k][i][j_1], swShift);

          L_TmpA = L_add(L_Cik, L_Cki);
          L_TmpA = L_mpy_ls(L_TmpA, pswRc[j]);

          L_Bik = L_shl(pppL_B[i][k][j_1], swShift);
          L_Fik = L_shl(pppL_F[i][k][j_1], swShift);

          L_TmpB = L_add(L_Bik, L_Fik);
          L_TmpB = L_mpy_ls(L_TmpB, pswRc[j]);

          /* Update F, B and C arrays */
          /*-----------------------------------*/

          L_temp = L_mpy_ls(L_Bik, swRcSq);
          L_temp = L_add(L_temp, L_Fik);
          pppL_F[i][k][j_0] = L_add(L_temp, L_TmpA);

          L_temp = L_mpy_ls(L_Fik, swRcSq);
          L_temp = L_add(L_temp, L_Bik);
          pppL_B[i - 1][k - 1][j_0] = L_add(L_temp, L_TmpA);

          L_temp = L_mpy_ls(L_Cki, swRcSq);
          L_temp = L_add(L_temp, L_Cik);
          pppL_C[i][k - 1][j_0] = L_add(L_temp, L_TmpB);

          L_temp = L_mpy_ls(L_Cik, swRcSq);
          L_temp = L_add(L_temp, L_Cki);
          pppL_C[k][i - 1][j_0] = L_add(L_temp, L_TmpB);

        }                              /* end of loop indexed by k */
        /*---------------------------*/

        k = NP - j - 1;

        /* Compute the common terms used by the FLAT recursion to reduce */
        /* computation.                                                  */
        /*---------------------------------------------------------------*/

        L_Cik = L_shl(pppL_C[i][k][j_1], swShift);
        L_Cki = L_shl(pppL_C[k][i][j_1], swShift);

        L_TmpA = L_add(L_Cik, L_Cki);
        L_TmpA = L_mpy_ls(L_TmpA, pswRc[j]);

        L_Bik = L_shl(pppL_B[i][k][j_1], swShift);
        L_Fik = L_shl(pppL_F[i][k][j_1], swShift);

        L_TmpB = L_add(L_Bik, L_Fik);
        L_TmpB = L_mpy_ls(L_TmpB, pswRc[j]);

        /* Update B and C arrays */
        /*-----------------------------------*/

        L_temp = L_mpy_ls(L_Fik, swRcSq);
        L_temp = L_add(L_temp, L_Bik);
        pppL_B[i - 1][k - 1][j_0] = L_add(L_temp, L_TmpA);

        L_temp = L_mpy_ls(L_Cki, swRcSq);
        L_temp = L_add(L_temp, L_Cik);
        pppL_C[i][k - 1][j_0] = L_add(L_temp, L_TmpB);

      }                                /* end of loop indexed by i */
      /*---------------------------*/

      i = NP - j - 1;
      for (k = i; k <= NP - j - 1; k++)
      {

        /* Compute the common terms used by the FLAT recursion to reduce */
        /* computation.                                                  */
        /*---------------------------------------------------------------*/

        L_Cik = L_shl(pppL_C[i][k][j_1], swShift);

        L_TmpA = L_add(L_Cik, L_Cik);
        L_TmpA = L_mpy_ls(L_TmpA, pswRc[j]);

        /* Update B array */
        /*-----------------------------------*/

        L_Bik = L_shl(pppL_B[i][k][j_1], swShift);
        L_Fik = L_shl(pppL_F[i][k][j_1], swShift);

        L_temp = L_mpy_ls(L_Fik, swRcSq);
        L_temp = L_add(L_temp, L_Bik);
        pppL_B[i - 1][k - 1][j_0] = L_add(L_temp, L_TmpA);

      }                                /* end of loop indexed by k */
      /*-----------------------------------------------------------*/

      /* OR the F and B matrix diagonals to find maximum for normalization */
      /*********************************************************************/

      L_TmpA = 0;
      for (l = 0; l <= NP - j - 2; l++)
      {
        L_TmpA |= pppL_F[l][l][j_0];
        L_TmpA |= pppL_B[l][l][j_0];
      }
      /* Compute the shift count to be applied to F, B, and C arrays */
      /* at the next lattice stage.                                  */
      /*-------------------------------------------------------------*/

      if (L_TmpA > 0)
      {
        swShift = norm_l(L_TmpA);
        swShift = sub(swShift, CVSHIFT);
      }
      else
      {
        swShift = 0;
      }

    }                                  /* End of update of F, B, and C
                                        * arrays for the next lattice stage */
    /*----------------------------------------------------------------*/

  }                                    /* Finished computation of
                                        * reflection coefficients */
  /*--------------------------------------------------------------*/

  return (L_R0);

}

/**************************************************************************
 *
 *   FUNCTION NAME: fnBest_CG
 *
 *   PURPOSE:
 *     The purpose of this function is to determine the C:G pair from the
 *     input arrays which maximize C*C/G
 *
 *   INPUTS:
 *
 *     pswCframe[0:siNumPairs]
 *
 *                     pointer to start of the C frame vector
 *
 *     pswGframe[0:siNumPairs]
 *
 *                     pointer to start of the G frame vector
 *
 *     pswCmaxSqr
 *
 *                     threshold Cmax**2 or 0 if no threshold
 *
 *     pswGmax
 *
 *                     threshold Gmax, must be > 0
 *
 *     siNumPairs
 *
 *                     number of C:G pairs to test
 *
 *   OUTPUTS:
 *
 *     pswCmaxSqr
 *
 *                     final Cmax**2 value
 *
 *     pswGmax
 *
 *                     final Gmax value
 *
 *   RETURN VALUE:
 *
 *     siMaxLoc
 *
 *                     index of Cmax in the input C matrix or -1 if none
 *
 *   DESCRIPTION:
 *
 *     test the result of (C * C * Gmax) - (Cmax**2 * G)
 *     if the result is > 0 then a new max has been found
 *     the Cmax**2, Gmax and MaxLoc parameters are all updated accordingly.
 *     if no new max is found for all NumPairs then MaxLoc will retain its
 *     original value
 *
 *   REFERENCE:  Sub-clause 4.1.8.1, 4.1.8.2, and 4.1.8.3 of GSM
 *     Recommendation 06.20
 *
 *   KEYWORDS: C_Frame, G_Frame, Cmax, Gmax, DELTA_LAGS, PITCH_LAGS
 *

****************************************************************************/

short int fnBest_CG(Shortword pswCframe[], Shortword pswGframe[],
                           Shortword *pswCmaxSqr, Shortword *pswGmax,
                           short int siNumPairs)
{

/*_________________________________________________________________________
|                                                                           |
|                            Automatic Variables                            |
|___________________________________________________________________________|
*/

  Longword L_Temp2;
  Shortword swCmaxSqr,
         swGmax,
         swTemp;
  short int siLoopCnt,
         siMaxLoc;

/*_________________________________________________________________________
|                                                                           |
|                              Executable Code                              |
|___________________________________________________________________________|
*/

  /* initialize */
  /* ---------- */

  swCmaxSqr = *pswCmaxSqr;
  swGmax = *pswGmax;
  siMaxLoc = -1;

  for (siLoopCnt = 0; siLoopCnt < siNumPairs; siLoopCnt++)
  {

    /* make sure both C and energy > 0 */
    /* ------------------------------- */

    if ((pswGframe[siLoopCnt] > 0) && (pswCframe[siLoopCnt] > 0))
    {

      /* calculate (C * C) */
      /* ----------------- */

      swTemp = mult_r(pswCframe[siLoopCnt], pswCframe[siLoopCnt]);

      /* calculate (C * C * Gmax) */
      /* ------------------------ */

      L_Temp2 = L_mult(swTemp, swGmax);

      /* calculate (C * C * Gmax) - (Cmax**2 * G) */
      /* ----------------------------------------- */

      L_Temp2 = L_msu(L_Temp2, swCmaxSqr, pswGframe[siLoopCnt]);

      /* if new max found, update it and its location */
      /* -------------------------------------------- */

      if (L_Temp2 > 0)
      {
        swCmaxSqr = swTemp;            /* Cmax**2 = current C * C */
        swGmax = pswGframe[siLoopCnt]; /* Gmax */
        siMaxLoc = siLoopCnt;          /* max location = current (C)
                                        * location */
      }
    }
  }

  /* set output */
  /* ---------- */

  *pswCmaxSqr = swCmaxSqr;
  *pswGmax = swGmax;
  return (siMaxLoc);
}

/***************************************************************************
 *
 *   FUNCTION NAME: fnExp2
 *
 *   PURPOSE:
 *     The purpose of this function is to implement a base two exponential
 *     2**(32*x) by polynomial approximation
 *
 *
 *   INPUTS:
 *
 *     L_Input
 *
 *                     unnormalized input exponent (input range constrained
 *                     to be < 0; for input < -0.46 output is 0)
 *
 *   OUTPUTS:
 *
 *     none
 *
 *   RETURN VALUE:
 *
 *     swTemp4
 *
 *                     exponential output
 *
 *   DESCRIPTION:
 *
 *     polynomial approximation is used for the generation of the exponential
 *
 *     2**(32*X) = 0.1713425*X*X + 0.6674432*X + 0.9979554
 *                     c2              c1            c0
 *
 *   REFERENCE:  Sub-clause 4.1.8.2 of GSM Recommendation 06.20, eqn 3.9
 *
 *   KEYWORDS: EXP2, DELTA_LAGS
 *
 *************************************************************************/

Shortword fnExp2(Longword L_Input)
{

/*_________________________________________________________________________
 |                                                                         |
 |                           Local Static Variables                        |
 |_________________________________________________________________________|
*/
  static Shortword pswPCoefE[3] =
  {                                    /* c2,   c1,    c0 */
    0x15ef, 0x556f, 0x7fbd
  };

/*_________________________________________________________________________
 |                                                                         |
 |                            Automatic Variables                          |
 |_________________________________________________________________________|
*/

  Shortword swTemp1,
         swTemp2,
         swTemp3,
         swTemp4;

/*_________________________________________________________________________
 |                                                                         |
 |                              Executable Code                            |
 |_________________________________________________________________________|
*/

  /* initialize */
  /* ---------- */

  swTemp3 = 0x0020;

  /* determine normlization shift count */
  /* ---------------------------------- */

  swTemp1 = extract_h(L_Input);
  L_Input = L_mult(swTemp1, swTemp3);
  swTemp2 = extract_h(L_Input);

  /* determine un-normalized shift count */
  /* ----------------------------------- */

  swTemp3 = -0x0001;
  swTemp4 = sub(swTemp3, swTemp2);

  /* normalize input */
  /* --------------- */

  L_Input = L_Input & LSP_MASK;
  L_Input = L_add(L_Input, L_deposit_h(swTemp3));

  L_Input = L_shr(L_Input, 1);
  swTemp1 = extract_l(L_Input);

  /* calculate x*x*c2 */
  /* ---------------- */

  swTemp2 = mult_r(swTemp1, swTemp1);
  L_Input = L_mult(swTemp2, pswPCoefE[0]);

  /* calculate x*x*c2 + x*c1 */
  /* ----------------------- */

  L_Input = L_mac(L_Input, swTemp1, pswPCoefE[1]);

  /* calculate x*x*c2 + x*c1 + c0 */
  /* --------------------------- */

  L_Input = L_add(L_Input, L_deposit_h(pswPCoefE[2]));

  /* un-normalize exponent if its requires it */
  /* ---------------------------------------- */

  if (swTemp4 > 0)
  {
    L_Input = L_shr(L_Input, swTemp4);
  }

  /* return result */
  /* ------------- */

  swTemp4 = extract_h(L_Input);
  return (swTemp4);
}

/***************************************************************************
 *
 *   FUNCTION NAME: fnLog2
 *
 *   PURPOSE:
 *     The purpose of this function is to take the log base 2 of input and
 *     divide by 32 and return; i.e. output = log2(input)/32
 *
 *   INPUTS:
 *
 *     L_Input
 *
 *                     input
 *
 *   OUTPUTS:
 *
 *     none
 *
 *   RETURN VALUE:
 *
 *     Shortword
 *
 *                     output
 *
 *   DESCRIPTION:
 *
 *     log2(x) = 4.0 * (-.3372223*x*x + .9981958*x -.6626105)
 *                           c0            c1          c2   (includes sign)
 *
 *   REFERENCE:  Sub-clause 4.1.8.2 of GSM Recommendation 06.20, eqn 3.9
 *
 *   KEYWORDS: log, logarithm, logbase2, fnLog2
 *
 *************************************************************************/

Shortword fnLog2(Longword L_Input)
{

/*_________________________________________________________________________
 |                                                                         |
 |                           Static Variables                              |
 |_________________________________________________________________________|
*/

  static Shortword
         swC0 = -0x2b2a,
         swC1 = 0x7fc5,
         swC2 = -0x54d0;

/*_________________________________________________________________________
 |                                                                         |
 |                            Automatic Variables                          |
 |_________________________________________________________________________|
*/

  short int siShiftCnt;
  Shortword swInSqrd,
         swIn;

/*_________________________________________________________________________
 |                                                                         |
 |                              Executable Code                            |
 |_________________________________________________________________________|
*/

  /* normalize input and store shifts required */
  /* ----------------------------------------- */

  siShiftCnt = norm_l(L_Input);
  L_Input = L_shl(L_Input, siShiftCnt);
  siShiftCnt = add(siShiftCnt, 1);
  siShiftCnt = negate(siShiftCnt);

  /* calculate x*x*c0 */
  /* ---------------- */

  swIn = extract_h(L_Input);
  swInSqrd = mult_r(swIn, swIn);
  L_Input = L_mult(swInSqrd, swC0);

  /* add x*c1 */
  /* --------- */

  L_Input = L_mac(L_Input, swIn, swC1);

  /* add c2 */
  /* ------ */

  L_Input = L_add(L_Input, L_deposit_h(swC2));

  /* apply *(4/32) */
  /* ------------- */

  L_Input = L_shr(L_Input, 3);
  L_Input = L_Input & 0x03ffffff;
  siShiftCnt = shl(siShiftCnt, 10);
  L_Input = L_add(L_Input, L_deposit_h(siShiftCnt));

  /* return log */
  /* ---------- */

  return (round(L_Input));
}

/***************************************************************************
 *
 *   FUNCTION NAME: getCCThreshold
 *
 *   PURPOSE:
 *     The purpose of this function is to calculate a threshold for other
 *     correlations (subject to limits), given subframe energy (Rp0),
 *     correlation squared (CC), and energy of delayed sequence (G)
 *
 *   INPUTS:
 *
 *     swRp0
 *
 *                     energy of the subframe
 *
 *     swCC
 *
 *                     correlation (squared) of subframe and delayed sequence
 *
 *     swG
 *
 *                     energy of delayed sequence
 *
 *   OUTPUTS:
 *
 *     none
 *
 *   RETURN VALUE:
 *
 *     swCCThreshold
 *
 *                     correlation (squared) threshold
 *
 *   DESCRIPTION:
 *
 *     CCt/0.5 = R - R(antilog(SCALE*log(max(CLAMP,(RG-CC)/RG))))
 *
 *     The threshold CCt is then applied with an understood value of Gt = 0.5
 *
 *   REFERENCE:  Sub-clause 4.1.8.2 of GSM Recommendation 06.20, eqn 3.9
 *
 *   KEYWORDS: getCCThreshold, getccthreshold, GET_CSQ_THRES
 *
 *************************************************************************/

Shortword getCCThreshold(Shortword swRp0,
                                Shortword swCC,
                                Shortword swG)
{

/*_________________________________________________________________________
 |                                                                         |
 |                            Automatic Variables                          |
 |_________________________________________________________________________|
*/

  Shortword swPGainClamp,
         swPGainScale,
         sw1;
  Longword L_1;

/*_________________________________________________________________________
 |                                                                         |
 |                              Executable Code                            |
 |_________________________________________________________________________|
*/

  /* load CLAMP and SCALE */
  /* -------------------- */

  swPGainClamp = PGAIN_CLAMP;
  swPGainScale = PGAIN_SCALE;

  /* calculate RG-CC */
  /* --------------- */

  L_1 = L_mult(swRp0, swG);
  sw1 = extract_h(L_1);
  L_1 = L_sub(L_1, L_deposit_h(swCC));

  /* if RG - CC > 0 do max(CLAMP, (RG-CC)/RG) */
  /* ---------------------------------------- */

  if (L_1 > 0)
  {

    sw1 = divide_s(extract_h(L_1), sw1);

    L_1 = L_deposit_h(sw1);

    if (sub(sw1, swPGainClamp) <= 0)
    {
      L_1 = L_deposit_h(swPGainClamp);
    }
  }
  /* else max(CLAMP, (RG-CC)/RG) is CLAMP */
  /* ------------------------------------ */

  else
  {
    L_1 = L_deposit_h(swPGainClamp);
  }

  /* L_1 holds max(CLAMP, (RG-CC)/RG)   */
  /* do antilog( SCALE * log( max() ) ) */
  /* ---------------------------------- */

  sw1 = fnLog2(L_1);

  L_1 = L_mult(sw1, swPGainScale);

  sw1 = fnExp2(L_1);


  /* do R - (R * antilog()) */
  /* ---------------------- */

  L_1 = L_deposit_h(swRp0);
  L_1 = L_msu(L_1, swRp0, sw1);

  /* apply Gt value */
  /* -------------- */

  L_1 = L_shr(L_1, 1);

  return (extract_h(L_1));
}


/***************************************************************************
 *
 *   FUNCTION NAME: getNWCoefs
 *
 *   PURPOSE:
 *
 *     Obtains best all-pole fit to various noise weighting
 *     filter combinations
 *
 *   INPUTS:
 *
 *     pswACoefs[0:9] - A(z) coefficient array
 *     psrNWCoefs[0:9] - filter smoothing coefficients
 *
 *   OUTPUTS:
 *
 *     pswHCoefs[0:9] - H(z) coefficient array
 *
 *   RETURN VALUE:
 *
 *     None
 *
 *   DESCRIPTION:
 *
 *     The function getNWCoefs() derives the spectral noise weighting
 *     coefficients W(z)and H(z).  W(z) and H(z) actually consist of
 *     three filters in cascade.  To avoid having such a complicated
 *     filter required for weighting, the filters are reduced to a
 *     single filter.
 *
 *     This is accomplished by passing an impulse through the cascased
 *     filters.  The impulse response of the filters is used to generate
 *     autocorrelation coefficients, which are then are transformed into
 *     a single direct form estimate of W(z) and H(z).  This estimate is
 *     called HHat(z) in the documentation.
 *
 *
 *   REFERENCE:  Sub-clause 4.1.7 of GSM Recommendation 06.20
 *
 *   KEYWORDS: spectral noise weighting, direct form coefficients
 *   KEYWORDS: getNWCoefs
 *
 *************************************************************************/

void   getNWCoefs(Shortword pswACoefs[],
                         Shortword pswHCoefs[])
{

/*_________________________________________________________________________
 |                                                                         |
 |                            Automatic Variables                          |
 |_________________________________________________________________________|
*/

  Shortword pswCoefTmp2[NP],
         pswCoefTmp3[NP],
         pswVecTmp[S_LEN],
         pswVecTmp2[S_LEN],
         pswTempRc[NP];
  Shortword swNormShift,
         iLoopCnt,
         iLoopCnt2;
  Longword pL_AutoCorTmp[NP + 1],
         L_Temp;
  short int siNum,
         k;

/*_________________________________________________________________________
 |                                                                         |
 |                              Executable Code                            |
 |_________________________________________________________________________|
*/

  /* Calculate smoothing parameters for all-zero filter */
  /* -------------------------------------------------- */

  for (iLoopCnt = 0; iLoopCnt < NP; iLoopCnt++)
  {
    pswCoefTmp2[iLoopCnt]
            = mult_r(psrNWCoefs[iLoopCnt], pswACoefs[iLoopCnt]);
  }

  /* Calculate smoothing parameters for all-pole filter */
  /* -------------------------------------------------- */

  for (iLoopCnt = 0; iLoopCnt < NP; iLoopCnt++)
  {
    pswCoefTmp3[iLoopCnt] = msu_r(0, psrNWCoefs[iLoopCnt + NP],
                                  pswACoefs[iLoopCnt]);
  }

  /* Get impulse response of 1st filter                             */
  /* Done by direct form IIR filter of order NP zero input response */
  /* -------------------------------------------------------------- */

  lpcIrZsIir(pswACoefs, pswVecTmp2);

  /* Send impulse response of 1st filter through 2nd filter */
  /* All-zero filter (FIR)                                  */
  /* ------------------------------------------------------ */

  lpcZsFir(pswVecTmp2, pswCoefTmp2, pswVecTmp);

  /* Send impulse response of 2nd filter through 3rd filter */
  /* All-pole filter (IIR)                                  */
  /* ------------------------------------------------------ */

  lpcZsIirP(pswVecTmp, pswCoefTmp3);

  /* Calculate energy in impulse response */
  /* ------------------------------------ */

  swNormShift = g_corr1(pswVecTmp, &L_Temp);

  pL_AutoCorTmp[0] = L_Temp;

  /* Calculate normalized autocorrelation function */
  /* --------------------------------------------- */

  for (k = 1; k <= NP; k++)
  {

    /* Calculate R(k), equation 2.31 */
    /* ----------------------------- */

    L_Temp = L_mult(pswVecTmp[0], pswVecTmp[0 + k]);

    for (siNum = S_LEN - k, iLoopCnt2 = 1; iLoopCnt2 < siNum; iLoopCnt2++)
    {
      L_Temp = L_mac(L_Temp, pswVecTmp[iLoopCnt2],
                     pswVecTmp[iLoopCnt2 + k]);
    }

    /* Normalize R(k) relative to R(0): */
    /* -------------------------------- */

    pL_AutoCorTmp[k] = L_shl(L_Temp, swNormShift);

  }


  /* Convert normalized autocorrelations to direct form coefficients */
  /* --------------------------------------------------------------- */

  aFlatRcDp(pL_AutoCorTmp, pswTempRc);
  rcToADp(ASCALE, pswTempRc, pswHCoefs);

}

/***************************************************************************
 *
 *   FUNCTION NAME: getNextVec
 *
 *   PURPOSE:
 *     The purpose of this function is to get the next vector in the list.
 *
 *   INPUTS: none
 *
 *   OUTPUTS: pswRc
 *
 *   RETURN VALUE: none
 *
 *   DESCRIPTION:
 *
 *     Both the quantizer and pre-quantizer are set concatenated 8 bit
 *     words.  Each of these words represents a reflection coefficient.
 *     The 8 bit words, are actually indices into a reflection
 *     coefficient lookup table.  Memory is organized in 16 bit words, so
 *     there are two reflection coefficients per ROM word.
 *
 *
 *     The full quantizer is subdivided into blocks.  Each of the
 *     pre-quantizers vectors "points" to a full quantizer block.  The
 *     vectors in a block, are comprised of either three or four
 *     elements.  These are concatenated, without leaving any space
 *     between them.
 *
 *     A block of full quantizer elements always begins on an even word.
 *     This may or may not leave a space depending on vector quantizer
 *     size.
 *
 *     getNextVec(), serves to abstract this arcane data format.  Its
 *     function is to simply get the next reflection coefficient vector
 *     in the list, be it a pre or full quantizer list.  This involves
 *     figuring out whether to pick the low or the high part of the 16
 *     bit ROM word.  As well as transforming the 8 bit stored value
 *     into a fractional reflection coefficient.  It also requires a
 *     setup routine to initialize iWordPtr and iWordHalfPtr, two
 *     variables global to this file.
 *
 *
 *
 *   REFERENCE:  Sub-clause 4.1.4.1 of GSM Recommendation 06.20
 *
 *   KEYWORDS: Quant quant vector quantizer
 *
 *************************************************************************/

void   getNextVec(Shortword pswRc[])
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
  i = iLow;

  if (iThree)
  {

    if (iWordHalfPtr == HIGH)
    {
      pswRc[i++] = psrSQuant[high(psrTable[iWordPtr])];
      pswRc[i++] = psrSQuant[low(psrTable[iWordPtr++])];
      pswRc[i] = psrSQuant[high(psrTable[iWordPtr])];
      iWordHalfPtr = LOW;
    }
    else
    {
      pswRc[i++] = psrSQuant[low(psrTable[iWordPtr++])];
      pswRc[i++] = psrSQuant[high(psrTable[iWordPtr])];
      pswRc[i] = psrSQuant[low(psrTable[iWordPtr++])];
      iWordHalfPtr = HIGH;
    }

  }
  else
  {
    pswRc[i++] = psrSQuant[high(psrTable[iWordPtr])];
    pswRc[i++] = psrSQuant[low(psrTable[iWordPtr++])];
    pswRc[i++] = psrSQuant[high(psrTable[iWordPtr])];
    pswRc[i] = psrSQuant[low(psrTable[iWordPtr++])];
  }
}

/***************************************************************************
 *
 *    FUNCTION NAME: getSfrmLpcTx
 *
 *    PURPOSE:
 *       Given frame information from past and present frame, interpolate
 *       (or copy) the frame based lpc coefficients into subframe
 *       lpc coeffs, i.e. the ones which will be used by the subframe
 *       as opposed to those coded and transmitted.
 *
 *    INPUT:
 *       swPrevR0,swNewR0 - Rq0 for the last frame and for this frame.
 *          These are the decoded values, not the codewords.
 *
 *       Previous lpc coefficients from the previous FRAME:
 *          in all filters below array[0] is the t=-1 element array[NP-1]
 *        t=-NP element.
 *       pswPrevFrmKs[NP] - decoded version of the rc's tx'd last frame
 *       pswPrevFrmAs[NP] - the above K's converted to A's.  i.e. direct
 *          form coefficients.
 *       pswPrevFrmSNWCoef[NP] - Coefficients for the Spectral Noise
 *          weighting filter from the previous frame
 *
 *       pswHPFSppech - pointer to High Pass Filtered Input speech
 *
 *       pswSoftInterp - a flag containing the soft interpolation
 *          decision.
 *
 *       Current lpc coefficients from the current frame:
 *       pswNewFrmKs[NP],pswNewFrmAs[NP],
 *       pswNewFrmSNWCoef[NP] - Spectral Noise Weighting Coefficients
 *           for the current frame
 *       ppswSNWCoefAs[1][NP] - pointer into a matrix containing
 *           the interpolated and uninterpolated LP Coefficient
 *           values for the Spectral Noise Weighting Filter.
 *
 *    OUTPUT:
 *       psnsSqrtRs[N_SUB] - a normalized number (struct NormSw)
 *          containing an estimate
 *          of RS for each subframe.  (number and a shift)
 *
 *       ppswSynthAs[N_SUM][NP] - filter coefficients used by the
 *          synthesis filter.
 *
 *    DESCRIPTION:
 *        For interpolated subframes, the direct form coefficients
 *        are converted to reflection coeffiecients to check for
 *        filter stability. If unstable, the uninterpolated coef.
 *        are used for that subframe.
 *
 *
 *    REFERENCE: Sub-clause of 4.1.6 and 4.1.7 of GSM Recommendation
 *        06.20
 *
 *    KEYWORDS: soft interpolation, int_lpc, interpolate, atorc, res_eng
 *
 *************************************************************************/


void   getSfrmLpcTx(Shortword swPrevR0, Shortword swNewR0,
/* last frm*/
                       Shortword pswPrevFrmKs[], Shortword pswPrevFrmAs[],
                           Shortword pswPrevFrmSNWCoef[],
/* this frm*/
                         Shortword pswNewFrmKs[], Shortword pswNewFrmAs[],
                           Shortword pswNewFrmSNWCoef[],
                           Shortword pswHPFSpeech[],
/* output */
                           short *pswSoftInterp,
                           struct NormSw *psnsSqrtRs,
               Shortword ppswSynthAs[][NP], Shortword ppswSNWCoefAs[][NP])
{

/*_________________________________________________________________________
 |                                                                         |
 |                            Automatic Variables                          |
 |_________________________________________________________________________|
*/

  Shortword swSi;
  Longword L_Temp;
  short int siSfrm,
         siStable,
         i;

/*_________________________________________________________________________
 |                                                                         |
 |                              Executable Code                            |
 |_________________________________________________________________________|
*/

  /* perform interpolation - both for synth filter and noise wgt filt */
  /*------------------------------------------------------------------*/
  siSfrm = 0;
  siStable = interpolateCheck(pswPrevFrmKs, pswPrevFrmAs,
                              pswPrevFrmAs, pswNewFrmAs,
                              psrOldCont[siSfrm], psrNewCont[siSfrm],
                              swPrevR0,
                              &psnsSqrtRs[siSfrm],
                              ppswSynthAs[siSfrm]);
  if (siStable)
  {
    for (i = 0; i < NP; i++)
    {
      L_Temp = L_mult(pswNewFrmSNWCoef[i], psrNewCont[siSfrm]);
      ppswSNWCoefAs[siSfrm][i] = mac_r(L_Temp, pswPrevFrmSNWCoef[i],
                                       psrOldCont[siSfrm]);
    }
  }
  else
  {

    /* this subframe is unstable */
    /*---------------------------*/

    for (i = 0; i < NP; i++)
    {
      ppswSNWCoefAs[siSfrm][i] = pswPrevFrmSNWCoef[i];
    }

  }

  /* interpolate subframes one and two */
  /*-----------------------------------*/

  for (siSfrm = 1; siSfrm < N_SUB - 1; siSfrm++)
  {
    siStable = interpolateCheck(pswNewFrmKs, pswNewFrmAs,
                                pswPrevFrmAs, pswNewFrmAs,
                                psrOldCont[siSfrm], psrNewCont[siSfrm],
                                swNewR0,
                                &psnsSqrtRs[siSfrm],
                                ppswSynthAs[siSfrm]);
    if (siStable)
    {
      for (i = 0; i < NP; i++)
      {
        L_Temp = L_mult(pswNewFrmSNWCoef[i], psrNewCont[siSfrm]);
        ppswSNWCoefAs[siSfrm][i] = mac_r(L_Temp, pswPrevFrmSNWCoef[i],
                                         psrOldCont[siSfrm]);
      }
    }
    else
    {
      /* this sfrm has unstable filter coeffs, would like to interp but
       * cant */
      /*--------------------------------------*/

      for (i = 0; i < NP; i++)
      {
        ppswSNWCoefAs[siSfrm][i] = pswNewFrmSNWCoef[i];
      }

    }
  }


  /* the last subframe: never interpolate. */
  /*--------------------------------------*/

  siSfrm = 3;
  for (i = 0; i < NP; i++)
  {
    ppswSNWCoefAs[siSfrm][i] = pswNewFrmSNWCoef[i];
    ppswSynthAs[siSfrm][i] = pswNewFrmAs[i];
  }

  /* calculate the residual energy for the last subframe */
  /*-----------------------------------------------------*/

  res_eng(pswNewFrmKs, swNewR0, &psnsSqrtRs[siSfrm]);



  /* done with interpolation, now compare the two sets of coefs.   */
  /* make the decision whether to interpolate (1) or not (0)       */
  /*---------------------------------------------------------------*/

  swSi = compResidEnergy(pswHPFSpeech,
                         ppswSynthAs, pswPrevFrmAs,
                         pswNewFrmAs, psnsSqrtRs);

  if (swSi == 0)
  {

    /* no interpolation done: copy the frame based data to output
     * coeffiecient arrays */

    siSfrm = 0;
    for (i = 0; i < NP; i++)
    {
      ppswSNWCoefAs[siSfrm][i] = pswPrevFrmSNWCoef[i];
      ppswSynthAs[siSfrm][i] = pswPrevFrmAs[i];
    }

    /* get RS (energy in the residual) for subframe 0 */
    /*------------------------------------------------*/

    res_eng(pswPrevFrmKs, swPrevR0, &psnsSqrtRs[siSfrm]);

    /* for subframe 1 and all subsequent sfrms, use lpc and R0 from new frm */
    /*---------------------------------------------------------------------*/

    res_eng(pswNewFrmKs, swNewR0, &psnsSqrtRs[1]);

    for (siSfrm = 2; siSfrm < N_SUB; siSfrm++)
      psnsSqrtRs[siSfrm] = psnsSqrtRs[1];

    for (siSfrm = 1; siSfrm < N_SUB; siSfrm++)
    {
      for (i = 0; i < NP; i++)
      {
        ppswSNWCoefAs[siSfrm][i] = pswNewFrmSNWCoef[i];
        ppswSynthAs[siSfrm][i] = pswNewFrmAs[i];
      }
    }
  }

  *pswSoftInterp = swSi;
}

/***************************************************************************
 *
 *    FUNCTION NAME: iir_d
 *
 *    PURPOSE:  Performs one second order iir section using double-precision.
 *              feedback,single precision xn and filter coefficients
 *
 *    INPUTS:
 *
 *      pswCoef[0:4]   An array of filter coefficients.
 *
 *      pswIn[0:159]   An array of input samples to be filtered, filtered
 *                     output samples written to the same array.
 *
 *      pswXstate[0:1] An array containing x-state memory.
 *
 *      pswYstate[0:3] An array containing y-state memory.
 *
 *      npts           Number of samples to filter (must be even).
 *
 *      shifts         number of shifts to be made on output y(n) before
 *                     storing to y(n) states.
 *
 *      swPreFirDownSh number of shifts apply to signal before the FIR.
 *
 *      swFinalUpShift number of shifts apply to signal before outputting.
 *
 *    OUTPUTS:
 *
 *       pswIn[0:159]  Output array containing the filtered input samples.
 *
 *    RETURN:
 *
 *       none.
 *
 *    DESCRIPTION:
 *
 *       Transfer function implemented:
 *         (b0 + b1*z-1 + b2*z-2)/(a0 - a1*z-1 - a2*z-2+
 *       data structure:
 *         Coeff array order:  b2,b1,b0,a2,a1
 *         Xstate array order: x(n-2),x(n-1)
 *         Ystate array order: y(n-2)MSB,y(n-2)LSB,y(n-1)MSB,y(n-1)LSB
 *
 *       There is no elaborate discussion of the filter, since it is
 *       trivial.
 *
 *       The filter's cutoff frequency is 120 Hz.
 *
 *    REFERENCE:  Sub-clause 4.1.1 GSM Recommendation 06.20
 *
 *    KEYWORDS: highpass filter, hp, HP, filter
 *
 *************************************************************************/

void   iir_d(Shortword pswCoeff[], Shortword pswIn[], Shortword pswXstate[],
                    Shortword pswYstate[], int npts, int shifts,
                    Shortword swPreFirDownSh, Shortword swFinalUpShift)
{

/*_________________________________________________________________________
 |                                                                         |
 |                            Automatic Variables                          |
 |_________________________________________________________________________|
*/

  int    loop_cnt;
  Longword L_sumA,
         L_sumB;
  Shortword swTemp,
         pswYstate_0,
         pswYstate_1,
         pswYstate_2,
         pswYstate_3,
         pswXstate_0,
         pswXstate_1,
         swx0,
         swx1;

/*_________________________________________________________________________
 |                                                                         |
 |                              Executable Code                            |
 |_________________________________________________________________________|
*/

  /* initialize the temporary state variables */
  /*------------------------------------------*/

  pswYstate_0 = pswYstate[0];
  pswYstate_1 = pswYstate[1];
  pswYstate_2 = pswYstate[2];
  pswYstate_3 = pswYstate[3];

  pswXstate_0 = pswXstate[0];
  pswXstate_1 = pswXstate[1];

  for (loop_cnt = 0; loop_cnt < npts; loop_cnt += 2)
  {

    swx0 = shr(pswIn[loop_cnt], swPreFirDownSh);
    swx1 = shr(pswIn[loop_cnt + 1], swPreFirDownSh);

    L_sumB = L_mult(pswYstate_1, pswCoeff[3]);
    L_sumB = L_mac(L_sumB, pswYstate_3, pswCoeff[4]);
    L_sumB = L_shr(L_sumB, 14);
    L_sumB = L_mac(L_sumB, pswYstate_0, pswCoeff[3]);
    L_sumB = L_mac(L_sumB, pswYstate_2, pswCoeff[4]);


    L_sumA = L_mac(L_sumB, pswCoeff[0], pswXstate_0);
    L_sumA = L_mac(L_sumA, pswCoeff[1], pswXstate_1);
    L_sumA = L_mac(L_sumA, pswCoeff[2], swx0);

    L_sumA = L_shl(L_sumA, shifts);

    pswXstate_0 = swx0;                /* Update X state x(n-1) <- x(n) */

    /* Update double precision Y state temporary variables */
    /*-----------------------------------------------------*/

    pswYstate_0 = extract_h(L_sumA);
    swTemp = extract_l(L_sumA);
    swTemp = shr(swTemp, 2);
    pswYstate_1 = 0x3fff & swTemp;

    /* Round, store output sample and increment to next input sample */
    /*---------------------------------------------------------------*/

    pswIn[loop_cnt] = round(L_shl(L_sumA, swFinalUpShift));

    L_sumB = L_mult(pswYstate_3, pswCoeff[3]);
    L_sumB = L_mac(L_sumB, pswYstate_1, pswCoeff[4]);
    L_sumB = L_shr(L_sumB, 14);
    L_sumB = L_mac(L_sumB, pswYstate_2, pswCoeff[3]);
    L_sumB = L_mac(L_sumB, pswYstate_0, pswCoeff[4]);


    L_sumA = L_mac(L_sumB, pswCoeff[0], pswXstate_1);
    L_sumA = L_mac(L_sumA, pswCoeff[1], pswXstate_0);
    L_sumA = L_mac(L_sumA, pswCoeff[2], swx1);

    L_sumA = L_shl(L_sumA, shifts);

    pswXstate_1 = swx1;                /* Update X state x(n-1) <- x(n) */

    /* Update double precision Y state temporary variables */
    /*-----------------------------------------------------*/

    pswYstate_2 = extract_h(L_sumA);

    swTemp = extract_l(L_sumA);
    swTemp = shr(swTemp, 2);
    pswYstate_3 = 0x3fff & swTemp;

    /* Round, store output sample and increment to next input sample */
    /*---------------------------------------------------------------*/

    pswIn[loop_cnt + 1] = round(L_shl(L_sumA, swFinalUpShift));

  }

  /* update the states for the next frame */
  /*--------------------------------------*/

  pswYstate[0] = pswYstate_0;
  pswYstate[1] = pswYstate_1;
  pswYstate[2] = pswYstate_2;
  pswYstate[3] = pswYstate_3;

  pswXstate[0] = pswXstate_0;
  pswXstate[1] = pswXstate_1;

}

/***************************************************************************
 *
 *    FUNCTION NAME: initPBarVBarFullL
 *
 *    PURPOSE:  Given the Longword normalized correlation sequence, function
 *              initPBarVBarL initializes the Longword initial condition
 *              arrays pL_PBarFull and pL_VBarFull for a full 10-th order LPC
 *              filter. It also shifts down the pL_VBarFull and pL_PBarFull
 *              arrays by a global constant RSHIFT bits. The pL_PBarFull and
 *              pL_VBarFull arrays are used to set the initial Shortword
 *              P and V conditions which are used in the actual search of the
 *              Rc prequantizer and the Rc quantizer.
 *
 *              This is an implementation of equations 4.14 and
 *              4.15.
 *
 *    INPUTS:
 *
 *        pL_CorrelSeq[0:NP]
 *                     A Longword normalized autocorrelation array computed
 *                     from unquantized reflection coefficients.
 *
 *       RSHIFT       The number of right shifts to be applied to the
 *                     input correlations prior to initializing the elements
 *                     of pL_PBarFull and pL_VBarFull output arrays. RSHIFT
 *                     is a global constant.
 *
 *    OUTPUTS:
 *
 *        pL_PBarFull[0:NP-1]
 *                     A Longword output array containing the P initial
 *                     conditions for the full 10-th order LPC filter.
 *                     The address of the 0-th element of  pL_PBarFull
 *                     is passed in when function initPBarVBarFullL is
 *                     called.
 *
 *        pL_VBarFull[-NP+1:NP-1]
 *                     A Longword output array containing the V initial
 *                     conditions for the full 10-th order LPC filter.
 *                     The address of the 0-th element of pL_VBarFull is
 *                     passed in when function initPBarVBarFullL is called.
 *    RETURN:
 *        none.
 *
 *    REFERENCE:  Sub-clause 4.1.4.1 GSM Recommendation 06.20
 *
 *************************************************************************/

void   initPBarFullVBarFullL(Longword pL_CorrelSeq[],
                                    Longword pL_PBarFull[],
                                    Longword pL_VBarFull[])
{

/*_________________________________________________________________________
 |                                                                         |
 |                            Automatic Variables                          |
 |_________________________________________________________________________|
*/

  int    i,
         bound;

/*_________________________________________________________________________
 |                                                                         |
 |                              Executable Code                            |
 |_________________________________________________________________________|
*/

  /* Initialize the AFLAT recursion PBarFull and VBarFull 32 bit arrays */
  /* for a 10-th order LPC filter.                                      */
  /*--------------------------------------------------------------------*/

  bound = NP - 1;

  for (i = 0; i <= bound; i++)
  {
    pL_PBarFull[i] = L_shr(pL_CorrelSeq[i], RSHIFT);
  }

  for (i = -bound; i < 0; i++)
  {
    pL_VBarFull[i] = pL_PBarFull[-i - 1];
  }

  for (i = 0; i < bound; i++)
  {
    pL_VBarFull[i] = pL_PBarFull[i + 1];
  }

  pL_VBarFull[bound] = L_shr(pL_CorrelSeq[bound + 1], RSHIFT);

}

/***************************************************************************
 *
 *    FUNCTION NAME: initPBarVBarL
 *
 *    PURPOSE:  Given the Longword pL_PBarFull array,
 *              function initPBarVBarL initializes the Shortword initial
 *              condition arrays pswPBar and pswVBar for a 3-rd order LPC
 *              filter, since the order of the 1st Rc-VQ segment is 3.
 *              The pswPBar and pswVBar arrays are a Shortword subset
 *              of the initial condition array pL_PBarFull.
 *              pswPBar and pswVBar are the initial conditions for the AFLAT
 *              recursion at a given segment. The AFLAT recursion is used
 *              to evaluate the residual error due to an Rc vector selected
 *              from a prequantizer or a quantizer.
 *
 *              This is an implementation of equation 4.18 and 4.19.
 *
 *    INPUTS:
 *
 *        pL_PBarFull[0:NP-1]
 *                     A Longword input array containing the P initial
 *                     conditions for the full 10-th order LPC filter.
 *                     The address of the 0-th element of  pL_PBarFull
 *                     is passed in when function initPBarVBarL is called.
 *
 *    OUTPUTS:
 *
 *        pswPBar[0:NP_AFLAT-1]
 *                     The output Shortword array containing the P initial
 *                     conditions for the P-V AFLAT recursion, set here
 *                     for the Rc-VQ search at the 1st Rc-VQ segment.
 *                     The address of the 0-th element of pswPBar is
 *                     passed in when function initPBarVBarL is called.
 *
 *        pswVBar[-NP_AFLAT+1:NP_AFLAT-1]
 *                     The output Shortword array containing the V initial
 *                     conditions for the P-V AFLAT recursion, set here
 *                     for the Rc-VQ search at the 1st Rc-VQ segment.
 *                     The address of the 0-th element of pswVBar is
 *                     passed in when function initPBarVBarL is called.
 *
 *    RETURN:
 *
 *        none.
 *
 *    REFERENCE:  Sub-clause 4.1.4.1 GSM Recommendation 06.20
 *
 *
 *************************************************************************/

void   initPBarVBarL(Longword pL_PBarFull[],
                            Shortword pswPBar[],
                            Shortword pswVBar[])
{

/*_________________________________________________________________________
 |                                                                         |
 |                            Automatic Variables                          |
 |_________________________________________________________________________|
*/

  int    bound,
         i;

/*_________________________________________________________________________
 |                                                                         |
 |                              Executable Code                            |
 |_________________________________________________________________________|
*/

  /* Initialize the AFLAT recursion P and V 16 bit arrays for a 3-rd    */
  /* order LPC filter corresponding to the 1-st reflection coefficient  */
  /* VQ segment. The PBar and VBar arrays store the initial conditions  */
  /* for the evaluating the residual error due to Rc vectors being      */
  /* evaluated from the Rc-VQ codebook at the 1-st Rc-VQ segment.       */
  /*--------------------------------------------------------------------*/
  bound = 2;

  for (i = 0; i <= bound; i++)
  {
    pswPBar[i] = round(pL_PBarFull[i]);
  }
  for (i = -bound; i < 0; i++)
  {
    pswVBar[i] = pswPBar[-i - 1];
  }
  for (i = 0; i < bound; i++)
  {
    pswVBar[i] = pswPBar[i + 1];
  }
  pswVBar[bound] = round(pL_PBarFull[bound + 1]);

}

/***************************************************************************
 *
 *   FUNCTION NAME: maxCCOverGWithSign
 *
 *   PURPOSE:
 *
 *     Finds lag which maximizes C^2/G ( C is allowed to be negative ).
 *
 *   INPUTS:
 *
 *     pswCIn[0:swNum-1]
 *
 *                     Array of C values
 *
 *     pswGIn[0:swNum-1]
 *
 *                     Array of G values
 *
 *     pswCCmax
 *
 *                     Initial value of CCmax
 *
 *     pswGmax
 *
 *                     Initial value of Gmax
 *
 *     swNum
 *
 *                     Number of lags to be searched
 *
 *   OUTPUTS:
 *
 *     pswCCmax
 *
 *                     Value of CCmax after search
 *
 *     pswGmax
 *
 *                     Value of Gmax after search
 *
 *   RETURN VALUE:
 *
 *     maxCCGIndex - index for max C^2/G, defaults to zero if all G <= 0
 *
 *   DESCRIPTION:
 *
 *     This routine is called from bestDelta().  The routine is a simple
 *     find the best in a list search.
 *
 *   REFERENCE:  Sub-clause 4.1.8.3 of GSM Recommendation 06.20
 *
 *   KEYWORDS: LTP correlation peak
 *
 *************************************************************************/

Shortword maxCCOverGWithSign(Shortword pswCIn[],
                                    Shortword pswGIn[],
                                    Shortword *pswCCMax,
                                    Shortword *pswGMax,
                                    Shortword swNum)
{

/*_________________________________________________________________________
 |                                                                         |
 |                            Automatic Variables                          |
 |_________________________________________________________________________|
*/

  Shortword swCC,
         i,
         maxCCGIndex,
         swCCMax,
         swGMax;
  Longword L_Temp;

/*_________________________________________________________________________
 |                                                                         |
 |                              Executable Code                            |
 |_________________________________________________________________________|
*/

  /* initialize max c^2/g to be the initial trajectory */
  /*---------------------------------------------------*/

  maxCCGIndex = 0;
  swGMax = pswGIn[0];

  if (pswCIn[0] < 0)
    swCCMax = negate(mult_r(pswCIn[0], pswCIn[0]));
  else
    swCCMax = mult_r(pswCIn[0], pswCIn[0]);

  for (i = 1; i < swNum; i++)
  {

    /* Imperfect interpolation can result in negative energies. */
    /* Check for this                                             */
    /*----------------------------------------------------------*/

    if (pswGIn[i] > 0)
    {

      swCC = mult_r(pswCIn[i], pswCIn[i]);

      if (pswCIn[i] < 0)
        swCC = negate(swCC);

      L_Temp = L_mult(swCC, swGMax);
      L_Temp = L_msu(L_Temp, pswGIn[i], swCCMax);

      /* Check if C^2*Gmax - G*Cmax^2 > 0 */
      /* -------------------------------- */

      if (L_Temp > 0)
      {
        swGMax = pswGIn[i];
        swCCMax = swCC;
        maxCCGIndex = i;
      }
    }
  }
  *pswGMax = swGMax;
  *pswCCMax = swCCMax;

  return (maxCCGIndex);
}                                      /* end of maxCCOverGWithSign */

/***************************************************************************
 *
 *   FUNCTION NAME: openLoopLagSearch
 *
 *   PURPOSE:
 *
 *     Determines voicing level for the frame.  If voiced, obtains list of
 *     lags to be searched in closed-loop lag search; and value of smoothed
 *     pitch and coefficient for harmonic-noise-weighting.
 *
 *   INPUTS:
 *
 *     pswWSpeech[-145:159] ( [-LSMAX:F_LEN-1] )
 *
 *                     W(z) filtered speech frame, with some history.
 *
 *     swPrevR0Index
 *
 *                     Index of R0 from the previous frame.
 *
 *     swCurrR0Index
 *
 *                     Index of R0 for the current frame.
 *
 *     psrLagTbl[0:255]
 *
 *                     Lag quantization table, in global ROM.
 *
 *     ppsrCGIntFilt[0:5][0:5] ( [tap][phase] )
 *
 *                     Interpolation smoothing filter for generating C(k)
 *                     and G(k) terms, where k is fractional.  Global ROM.
 *
 *     swSP        
 *                     speech flag, required for DTX mode
 *
 *   OUTPUTS:
 *
 *     psiUVCode
 *
 *                     (Pointer to) the coded voicing level.
 *
 *     pswLagList[0:?]
 *
 *                     Array of lags to use in the search of the adaptive
 *                     codebook (long-term predictor).  Length determined
 *                     by pswNumLagList[].
 *
 *     pswNumLagList[0:3] ( [0:N_SUB-1] )
 *
 *                     Array of number of lags to use in search of adaptive
 *                     codebook (long-term predictor) for each subframe.
 *
 *     pswPitchBuf[0:3] ( [0:N_SUB-1] )
 *
 *                     Array of estimates of pitch, to be used in harmonic-
 *                     noise-weighting.
 *
 *     pswHNWCoefBuf[0:3] ( [0:N_SUB-1] )
 *
 *                     Array of harmonic-noise-weighting coefficients.
 *
 *     psnsWSfrmEng[-4:3] ( [-N_SUB:N_SUB-1] )
 *
 *                     Array of energies of weighted speech (input speech
 *                     sent through W(z) weighting filter), each stored as
 *                     normalized fraction and shift count.  The zero index
 *                     corresponds to the first subframe of the current
 *                     frame, so there is some history represented.  The
 *                     energies are used for scaling purposes only.
 *
 *     pswVadLags[4]
 *
 *                     An array of Shortwords containing the best open
 *                     loop LTP lags for the four subframes.

 *
 *   DESCRIPTION:
 *
 *     Scaling is done on the input weighted speech, such that the C(k) and
 *     G(k) terms will all be representable.  The amount of scaling is
 *     determined by the maximum energy of any subframe of weighted speech
 *     from the current frame or last frame.  These energies are maintained
 *     in a buffer, and used for scaling when the excitation is determined
 *     later in the analysis.
 *
 *     This function is the main calling program for the open loop lag
 *     search.
 *
 *   REFERENCE:  Sub-clauses 4.1.8.1-4.1.8.4 of GSM Recommendation 06.20
 *
 *   Keywords: openlooplagsearch, openloop, lag, pitch
 *
 **************************************************************************/



void   openLoopLagSearch(Shortword pswWSpeech[],
                                Shortword swPrevR0Index,
                                Shortword swCurrR0Index,
                                Shortword *psiUVCode,
                                Shortword pswLagList[],
                                Shortword pswNumLagList[],
                                Shortword pswPitchBuf[],
                                Shortword pswHNWCoefBuf[],
                                struct NormSw psnsWSfrmEng[],
                                Shortword pswVadLags[],
                                Shortword swSP)
{

/*_________________________________________________________________________
 |                                                                         |
 |                            Automatic Variables                          |
 |_________________________________________________________________________|
*/

  Longword L_WSfrmEng,
         L_G,
         L_C,
         L_Voicing;
  Shortword swBestPG,
         swCCMax,
         swGMax,
         swCCDivG;
  Shortword swTotalCCDivG,
         swCC,
         swG,
         swRG;
  short  i,
         j,
         k,
         siShift,
         siIndex,
         siTrajIndex,
         siAnchorIndex;
  short  siNumPeaks,
         siNumTrajToDo,
         siPeakIndex,
         siFIndex;
  short  siNumDelta,
         siBIndex,
         siBestTrajIndex = 0;
  short  siLowestSoFar,
         siLagsSoFar,
         si1,
         si2,
         si3;
  struct NormSw snsMax;

  Shortword pswGFrame[G_FRAME_LEN];
  Shortword *ppswGSfrm[N_SUB];
  Shortword pswSfrmEng[N_SUB];
  Shortword pswCFrame[C_FRAME_LEN];
  Shortword *ppswCSfrm[N_SUB];
  Shortword pswLIntBuf[N_SUB];
  Shortword pswCCThresh[N_SUB];
  Shortword pswScaledWSpeechBuffer[W_F_BUFF_LEN];
  Shortword *pswScaledWSpeech;
  Shortword ppswTrajLBuf[N_SUB * NUM_TRAJ_MAX][N_SUB];
  Shortword ppswTrajCCBuf[N_SUB * NUM_TRAJ_MAX][N_SUB];
  Shortword ppswTrajGBuf[N_SUB * NUM_TRAJ_MAX][N_SUB];
  Shortword pswLPeaks[2 * LMAX / LMIN];
  Shortword pswCPeaks[2 * LMAX / LMIN];
  Shortword pswGPeaks[2 * LMAX / LMIN];
  Shortword pswLArray[DELTA_LEVELS];

  pswScaledWSpeech = pswScaledWSpeechBuffer + LSMAX;

/*_________________________________________________________________________
 |                                                                         |
 |                            Executable Code                              |
 |_________________________________________________________________________|
*/

  /* Scale the weighted speech so that all correlations and energies */
  /* will be less than 1.0 in magnitude.  The scale factor is        */
  /* determined by the maximum energy of any subframe contained in   */
  /* the weighted speech buffer                                      */
  /*-----------------------------------------------------------------*/

  /* Perform one frame of delay on the subframe energy array */
  /*---------------------------------------------------------*/

  for (i = 0; i < N_SUB; i++)
    psnsWSfrmEng[i - N_SUB] = psnsWSfrmEng[i];

  /* Calculate the subframe energies of the current weighted speech frame. */
  /* Overflow protection is done based on the energy in the LPC analysis  */
  /* window (previous or current) which is closest to the subframe.       */
  /*----------------------------------------------------------------------*/

  psnsWSfrmEng[0].sh = g_corr1s(&pswWSpeech[0],
                                r0BasedEnergyShft(swPrevR0Index),
                                &L_WSfrmEng);
  psnsWSfrmEng[0].man = round(L_WSfrmEng);

  psnsWSfrmEng[1].sh = g_corr1s(&pswWSpeech[S_LEN],
                                r0BasedEnergyShft(swCurrR0Index),
                                &L_WSfrmEng);
  psnsWSfrmEng[1].man = round(L_WSfrmEng);

  psnsWSfrmEng[2].sh = g_corr1s(&pswWSpeech[2 * S_LEN],
                                r0BasedEnergyShft(swCurrR0Index),
                                &L_WSfrmEng);
  psnsWSfrmEng[2].man = round(L_WSfrmEng);

  psnsWSfrmEng[3].sh = g_corr1s(&pswWSpeech[3 * S_LEN],
                                r0BasedEnergyShft(swCurrR0Index),
                                &L_WSfrmEng);
  psnsWSfrmEng[3].man = round(L_WSfrmEng);

  /* Find the maximum weighted speech subframe energy from all values */
  /* in the array (the array includes the previous frame's subframes, */
  /* and the current frame's subframes)                               */
  /*------------------------------------------------------------------*/

  snsMax.man = 0;
  snsMax.sh = 0;
  for (i = -N_SUB; i < N_SUB; i++)
  {

    if (psnsWSfrmEng[i].man > 0)
    {

      if (snsMax.man == 0)
        snsMax = psnsWSfrmEng[i];

      if (sub(psnsWSfrmEng[i].sh, snsMax.sh) < 0)
        snsMax = psnsWSfrmEng[i];

      if (sub(psnsWSfrmEng[i].sh, snsMax.sh) == 0 &&
          sub(psnsWSfrmEng[i].man, snsMax.man) > 0)
        snsMax = psnsWSfrmEng[i];
    }
  }

  /* Now scale speech up or down, such that the maximum subframe */
  /* energy value will be in range [0.125, 0.25).  This gives a  */
  /* little room for other maxima, and interpolation filtering   */
  /*-------------------------------------------------------------*/

  siShift = sub(shr(snsMax.sh, 1), 1);

  for (i = 0; i < W_F_BUFF_LEN; i++)
    pswScaledWSpeech[i - LSMAX] = shl(pswWSpeech[i - LSMAX], siShift);

  /* Calculate the G(k) (k an integer) terms for all subframes.  (A note */
  /* on the organization of the G buffer: G(k) for a given subframe is   */
  /* the energy in the weighted speech sequence of length S_LEN (40)     */
  /* which begins k back from the beginning of the given subframe-- that */
  /* is, it begins at a lag of k.  These sequences overlap from one      */
  /* subframe to the next, so it is only necessary to compute and store  */
  /* the unique energies.  The unique energies are computed and stored   */
  /* in this buffer, and pointers are assigned for each subframe to make */
  /* array indexing for each subframe easier.                            */
  /* */
  /* (Terms in the G buffer are in order of increasing k, so the energy  */
  /* of the first sequence-- that is, the oldest sequence-- in the       */
  /* weighted speech buffer appears at the end of the G buffer.          */
  /* */
  /* (The subframe pointers are assigned so that they point to the first */
  /* k for their respective subframes, k = LSMIN.)                       */
  /*---------------------------------------------------------------------*/

  L_G = 0;
  for (i = -LSMAX; i < -LSMAX + S_LEN; i++)
    L_G = L_mac(L_G, pswScaledWSpeech[i], pswScaledWSpeech[i]);

  pswGFrame[G_FRAME_LEN - 1] = extract_h(L_G);

  for (i = -LSMAX; i < G_FRAME_LEN - LSMAX - 1; i++)
  {

    L_G = L_msu(L_G, pswScaledWSpeech[i], pswScaledWSpeech[i]);
    L_G = L_mac(L_G, pswScaledWSpeech[i + S_LEN],
                pswScaledWSpeech[i + S_LEN]);
    pswGFrame[G_FRAME_LEN - LSMAX - 2 - i] = extract_h(L_G);
  }

  ppswGSfrm[0] = pswGFrame + 3 * S_LEN;
  ppswGSfrm[1] = pswGFrame + 2 * S_LEN;
  ppswGSfrm[2] = pswGFrame + S_LEN;
  ppswGSfrm[3] = pswGFrame;

  /* Copy the G(k) terms which also happen to be the subframe energies; */
  /* calculate the 4th subframe energy, which is not a G(k)             */
  /*--------------------------------------------------------------------*/

  pswSfrmEng[0] = pswGFrame[G_FRAME_LEN - 1 - LSMAX];
  pswSfrmEng[1] = pswGFrame[G_FRAME_LEN - 1 - LSMAX - S_LEN];
  pswSfrmEng[2] = pswGFrame[G_FRAME_LEN - 1 - LSMAX - 2 * S_LEN];

  L_WSfrmEng = 0;
  for (i = F_LEN - S_LEN; i < F_LEN; i++)
    L_WSfrmEng = L_mac(L_WSfrmEng, pswScaledWSpeech[i], pswScaledWSpeech[i]);

  pswSfrmEng[3] = extract_h(L_WSfrmEng);

  /* Calculate the C(k) (k an integer) terms for all subframes. */
  /* (The C(k) terms are all unique, so there is no overlapping */
  /* as in the G buffer.)                                       */
  /*------------------------------------------------------------*/

  for (i = 0; i < N_SUB; i++)
  {

    for (j = LSMIN; j <= LSMAX; j++)
    {

      L_C = 0;
      for (k = 0; k < S_LEN; k++)
      {

        L_C = L_mac(L_C, pswScaledWSpeech[i * S_LEN + k],
                    pswScaledWSpeech[i * S_LEN - j + k]);
      }

      pswCFrame[i * CG_TERMS + j - LSMIN] = extract_h(L_C);
    }
  }

  ppswCSfrm[0] = pswCFrame;
  ppswCSfrm[1] = pswCFrame + CG_TERMS;
  ppswCSfrm[2] = pswCFrame + 2 * CG_TERMS;
  ppswCSfrm[3] = pswCFrame + 3 * CG_TERMS;

  /* For each subframe: find the max C(k)*C(k)/G(k) where C(k) > 0 and */
  /* k is integer; save the corresponding k; calculate the             */
  /* threshold against which other peaks in the interpolated CC/G      */
  /* sequence will be checked.  Meanwhile, accumulate max CC/G over    */
  /* the frame for the voiced/unvoiced determination.                  */
  /*-------------------------------------------------------------------*/

  swBestPG = 0;
  for (i = 0; i < N_SUB; i++)
  {

    /* Find max CC/G (C > 0), store corresponding k */
    /*----------------------------------------------*/

    swCCMax = 0;
    swGMax = 0x003f;
    siIndex = fnBest_CG(&ppswCSfrm[i][LMIN - LSMIN],
                        &ppswGSfrm[i][LMIN - LSMIN],
                        &swCCMax, &swGMax,
                        LMAX - LMIN + 1);

    if (siIndex == -1)
    {
      pswLIntBuf[i] = 0;
      pswVadLags[i] = LMIN;            /* store lag value for VAD algorithm */
    }
    else
    {
      pswLIntBuf[i] = add(LMIN, (Shortword) siIndex);
      pswVadLags[i] = pswLIntBuf[i];   /* store lag value for VAD algorithm */
    }

    if (pswLIntBuf[i] > 0)
    {

      /* C > 0 was found: accumulate CC/G, get threshold */
      /*-------------------------------------------------*/

      if (sub(swCCMax, swGMax) < 0)
        swCCDivG = divide_s(swCCMax, swGMax);
      else
        swCCDivG = SW_MAX;

      swBestPG = add(swCCDivG, swBestPG);

      pswCCThresh[i] = getCCThreshold(pswSfrmEng[i], swCCMax, swGMax);
    }
    else
      pswCCThresh[i] = 0;
  }

  /* Make voiced/unvoiced decision */
  /*-------------------------------*/

  L_Voicing = 0;
  for (i = 0; i < N_SUB; i++)
    L_Voicing = L_mac(L_Voicing, pswSfrmEng[i], UV_SCALE0);

  L_Voicing = L_add(L_Voicing, L_deposit_h(swBestPG));

  if ( (L_Voicing <= 0) || (swSP == 0) )
  {

    /* Unvoiced: set return values to zero */
    /*-------------------------------------*/

    for (i = 0; i < N_SUB; i++)
    {

      pswNumLagList[i] = 0;
      pswLagList[i] = 0;
      pswPitchBuf[i] = 0;
      pswHNWCoefBuf[i] = 0;
    }

    *psiUVCode = 0;
  }
  else
  {

    /* Voiced: get best delta-codeable lag trajectory; find pitch and */
    /* harmonic-noise-weighting coefficients for each subframe        */
    /*----------------------------------------------------------------*/

    siTrajIndex = 0;
    swBestPG = SW_MIN;
    for (siAnchorIndex = 0; siAnchorIndex < N_SUB; siAnchorIndex++)
    {

      /* Call pitchLags: for the current subframe, find peaks in the */
      /* C(k)**2/G(k) (k fractional) function which exceed the       */
      /* threshold set by the maximum C(k)**2/G(k) (k integer)       */
      /* (also get the smoothed pitch and harmonic-noise-weighting   */
      /* coefficient).                                               */
      /* */
      /* If there is no C(k) > 0 (k integer), set the smoothed pitch */
      /* to its minimum value and set the harmonic-noise-weighting   */
      /* coefficient to zero.                                        */
      /*-------------------------------------------------------------*/

      if (pswLIntBuf[siAnchorIndex] != 0)
      {

        pitchLags(pswLIntBuf[siAnchorIndex],
                  ppswCSfrm[siAnchorIndex],
                  ppswGSfrm[siAnchorIndex],
                  pswCCThresh[siAnchorIndex],
                  pswLPeaks,
                  pswCPeaks,
                  pswGPeaks,
                  &siNumPeaks,
                  &pswPitchBuf[siAnchorIndex],
                  &pswHNWCoefBuf[siAnchorIndex]);

        siPeakIndex = 0;
      }
      else
      {

        /* No C(k) > 0 (k integer): set pitch to min, coef to zero, */
        /* go to next subframe.                                     */
        /*----------------------------------------------------------*/

        pswPitchBuf[siAnchorIndex] = LMIN_FR;
        pswHNWCoefBuf[siAnchorIndex] = 0;
        continue;
      }

      /* It is possible that by interpolating, the only positive */
      /* C(k) was made negative.  Check for this here            */
      /*---------------------------------------------------------*/

      if (siNumPeaks == 0)
      {

        pswPitchBuf[siAnchorIndex] = LMIN_FR;
        pswHNWCoefBuf[siAnchorIndex] = 0;
        continue;
      }

      /* Peak(s) found in C**2/G function: find the best delta-codeable */
      /* trajectory through each peak (unless that peak has already     */
      /* paritcipated in a trajectory) up to a total of NUM_TRAJ_MAX    */
      /* peaks.  After each trajectory is found, check to see if that   */
      /* trajectory is the best one so far.                             */
      /*----------------------------------------------------------------*/

      if (siNumPeaks > NUM_TRAJ_MAX)
        siNumTrajToDo = NUM_TRAJ_MAX;
      else
        siNumTrajToDo = siNumPeaks;

      while (siNumTrajToDo)
      {

        /* Check if this peak has already participated in a trajectory. */
        /* If so, skip it, decrement the number of trajectories yet to  */
        /* be evaluated, and go on to the next best peak                */
        /*--------------------------------------------------------------*/

        si1 = 0;
        for (i = 0; i < siTrajIndex; i++)
        {

          if (sub(pswLPeaks[siPeakIndex],
                  ppswTrajLBuf[i][siAnchorIndex]) == 0)
            si1 = 1;
        }

        if (si1)
        {

          siNumTrajToDo--;
          if (siNumTrajToDo)
          {

            siPeakIndex++;
            continue;
          }
          else
            break;
        }

        /* The current peak is not in a previous trajectory: find */
        /* the best trajectory through it.                        */
        /* */
        /* First, store the lag, C**2, and G for the peak in the  */
        /* trajectory storage buffer                              */
        /*--------------------------------------------------------*/

        ppswTrajLBuf[siTrajIndex][siAnchorIndex] = pswLPeaks[siPeakIndex];
        ppswTrajGBuf[siTrajIndex][siAnchorIndex] = pswGPeaks[siPeakIndex];
        ppswTrajCCBuf[siTrajIndex][siAnchorIndex] =
                mult_r(pswCPeaks[siPeakIndex], pswCPeaks[siPeakIndex]);

        /* Complete the part of the trajectory that extends forward */
        /* from the anchor subframe                                 */
        /*----------------------------------------------------------*/

        for (siFIndex = siAnchorIndex + 1; siFIndex < N_SUB; siFIndex++)
        {

          /* Get array of lags which are delta-codeable */
          /* */
          /* First, get code for largest lag in array   */
          /* (limit it)                                 */
          /*--------------------------------------------*/

          quantLag(ppswTrajLBuf[siTrajIndex][siFIndex - 1],
                   &si1);
          si2 = add(si1, (DELTA_LEVELS / 2 - 1) - (NUM_CLOSED - 1));
          if (sub(si2, (1 << L_BITS) - 1) > 0)
            si2 = (1 << L_BITS) - 1;

          /* Get code for smallest lag in array (limit it) */
          /*-----------------------------------------------*/

          si3 = sub(si1, (DELTA_LEVELS / 2) - (NUM_CLOSED - 1));
          if (si3 < 0)
            si3 = 0;

          /* Generate array of lags */
          /*------------------------*/

          for (i = si3, j = 0; i <= si2; i++, j++)
            pswLArray[j] = psrLagTbl[i];

          siNumDelta = add(sub(si2, si3), 1);

          /* Search array of delta-codeable lags for one which maximizes */
          /* C**2/G                                                      */
          /*-------------------------------------------------------------*/

          bestDelta(pswLArray, ppswCSfrm[siFIndex], ppswGSfrm[siFIndex],
                    siNumDelta, siFIndex,
                    ppswTrajLBuf[siTrajIndex], ppswTrajCCBuf[siTrajIndex],
                    ppswTrajGBuf[siTrajIndex]);
        }

        /* Complete the part of the trajectory that extends backward */
        /* from the anchor subframe                                  */
        /*-----------------------------------------------------------*/

        for (siBIndex = siAnchorIndex - 1; siBIndex >= 0; siBIndex--)
        {

          /* Get array of lags which are delta-codeable */
          /* */
          /* First, get code for largest lag in array   */
          /* (limit it)                                 */
          /*--------------------------------------------*/

          quantLag(ppswTrajLBuf[siTrajIndex][siBIndex + 1],
                   &si1);
          si2 = add(si1, (DELTA_LEVELS / 2) - (NUM_CLOSED - 1));
          if (sub(si2, (1 << L_BITS) - 1) > 0)
            si2 = (1 << L_BITS) - 1;

          /* Get code for smallest lag in array (limit it) */
          /*-----------------------------------------------*/

          si3 = sub(si1, (DELTA_LEVELS / 2 - 1) - (NUM_CLOSED - 1));
          if (si3 < 0)
            si3 = 0;

          /* Generate array of lags */
          /*------------------------*/

          for (i = si3, j = 0; i <= si2; i++, j++)
            pswLArray[j] = psrLagTbl[i];

          siNumDelta = add(sub(si2, si3), 1);

          /* Search array of delta-codeable lags for one which maximizes */
          /* C**2/G                                                      */
          /*-------------------------------------------------------------*/

          bestDelta(pswLArray, ppswCSfrm[siBIndex], ppswGSfrm[siBIndex],
                    siNumDelta, siBIndex,
                    ppswTrajLBuf[siTrajIndex], ppswTrajCCBuf[siTrajIndex],
                    ppswTrajGBuf[siTrajIndex]);
        }

        /* This trajectory done: check total C**2/G for this trajectory */
        /* against current best trajectory                              */
        /* */
        /* Get total C**2/G for this trajectory                         */
        /*--------------------------------------------------------------*/

        swTotalCCDivG = 0;
        for (i = 0; i < N_SUB; i++)
        {

          swCC = ppswTrajCCBuf[siTrajIndex][i];
          swG = ppswTrajGBuf[siTrajIndex][i];

          if (swG <= 0)
          {

            /* Negative G (imperfect interpolation): do not include in */
            /* total                                                   */
            /*---------------------------------------------------------*/

            swCCDivG = 0;
          }
          else if (sub(abs_s(swCC), swG) > 0)
          {

            /* C**2/G > 0: limit quotient, add to total */
            /*------------------------------------------*/

            if (swCC > 0)
              swCCDivG = SW_MAX;
            else
              swCCDivG = SW_MIN;

            swTotalCCDivG = add(swTotalCCDivG, swCCDivG);
          }
          else
          {

            /* accumulate C**2/G */
            /*-------------------*/

            if (swCC < 0)
            {

              swCCDivG = divide_s(negate(swCC), swG);
              swTotalCCDivG = sub(swTotalCCDivG, swCCDivG);
            }
            else
            {

              swCCDivG = divide_s(swCC, swG);
              swTotalCCDivG = add(swTotalCCDivG, swCCDivG);
            }
          }
        }

        /* Compare this trajectory with current best, update if better */
        /*-------------------------------------------------------------*/

        if (sub(swTotalCCDivG, swBestPG) > 0)
        {

          swBestPG = swTotalCCDivG;
          siBestTrajIndex = siTrajIndex;
        }

        /* Update trajectory index, peak index, decrement the number */
        /* of trajectories left to do.                               */
        /*-----------------------------------------------------------*/

        siTrajIndex++;
        siPeakIndex++;
        siNumTrajToDo--;
      }
    }

    if (siTrajIndex == 0)
    {

      /* No trajectories searched despite voiced designation: change */
      /* designation to unvoiced.                                    */
      /*-------------------------------------------------------------*/

      for (i = 0; i < N_SUB; i++)
      {

        pswNumLagList[i] = 0;
        pswLagList[i] = 0;
        pswPitchBuf[i] = 0;
        pswHNWCoefBuf[i] = 0;
      }

      *psiUVCode = 0;
    }
    else
    {

      /* Best trajectory determined: get voicing level, generate the */
      /* constrained list of lags to search in the adaptive codebook */
      /* for each subframe                                           */
      /* */
      /* First, get voicing level                                    */
      /*-------------------------------------------------------------*/

      *psiUVCode = 3;
      siLowestSoFar = 2;
      for (i = 0; i < N_SUB; i++)
      {

        /* Check this subframe against highest voicing threshold */
        /*-------------------------------------------------------*/

        swCC = ppswTrajCCBuf[siBestTrajIndex][i];
        swG = ppswTrajGBuf[siBestTrajIndex][i];

        swRG = mult_r(swG, pswSfrmEng[i]);
        L_Voicing = L_deposit_h(swCC);
        L_Voicing = L_mac(L_Voicing, swRG, UV_SCALE2);

        if (L_Voicing < 0)
        {

          /* Voicing for this subframe failed to meet 2/3 threshold:  */
          /* therefore, voicing level for entire frame can only be as */
          /* high as 2                                                */
          /*----------------------------------------------------------*/

          *psiUVCode = siLowestSoFar;

          L_Voicing = L_deposit_h(swCC);
          L_Voicing = L_mac(L_Voicing, swRG, UV_SCALE1);

          if (L_Voicing < 0)
          {

            /* Voicing for this subframe failed to meet 1/2 threshold: */
            /* therefore, voicing level for entire frame can only be   */
            /* as high as 1                                            */
            /*---------------------------------------------------------*/

            *psiUVCode = siLowestSoFar = 1;
          }
        }
      }

      /* Generate list of lags to be searched in closed-loop */
      /*-----------------------------------------------------*/

      siLagsSoFar = 0;
      for (i = 0; i < N_SUB; i++)
      {

        quantLag(ppswTrajLBuf[siBestTrajIndex][i], &si1);

        si2 = add(si1, NUM_CLOSED / 2);
        if (sub(si2, (1 << L_BITS) - 1) > 0)
          si2 = (1 << L_BITS) - 1;

        si3 = sub(si1, NUM_CLOSED / 2);
        if (si3 < 0)
          si3 = 0;

        pswNumLagList[i] = add(sub(si2, si3), 1);

        for (j = siLagsSoFar; j < siLagsSoFar + pswNumLagList[i]; j++)
          pswLagList[j] = psrLagTbl[si3++];

        siLagsSoFar += pswNumLagList[i];
      }
    }
  }
}                                      /* end of openLoopLagSearch */

/***************************************************************************
 *
 *   FUNCTION NAME: pitchLags
 *
 *   PURPOSE:
 *
 *     Locates peaks in the interpolated C(k)*C(k)/G(k) sequence for a
 *     subframe which exceed a given threshold. Also determines the
 *     fundamental pitch, and a harmonic-noise-weighting coefficient.
 *
 *   INPUTS:
 *
 *     swBestIntLag
 *
 *                     The integer lag for which CC/G is maximum.
 *
 *     pswIntCs[0:127]
 *
 *                     The C(k) sequence for the subframe, with k an integer.
 *
 *     pswIntGs[0:127]
 *
 *                     The G(k) sequence for the subframe, with k an integer.
 *
 *     swCCThreshold
 *
 *                     The CC/G threshold which peaks must exceed (G is
 *                     understood to be 0.5).
 *
 *     psrLagTbl[0:255]
 *
 *                     The lag quantization table.
 *
 *
 *   OUTPUTS:
 *
 *     pswLPeaksSorted[0:10(max)]
 *
 *                     List of fractional lags where CC/G peaks, highest
 *                     peak first.
 *
 *     pswCPeaksSorted[0:10(max)]
 *
 *                     List of C's where CC/G peaks.
 *
 *     pswGPeaksSorted[0:10(max)]
 *
 *                     List of G's where CC/G peaks.
 *
 *     psiNumSorted
 *
 *                     Number of peaks found.
 *
 *     pswPitch
 *
 *                     The fundamental pitch for current subframe.
 *
 *     pswHNWCoef
 *
 *                     The harmonic-noise-weighting coefficient for the
 *                     current subframe.
 *
 *   RETURN VALUE:
 *
 *     None
 *
 *   DESCRIPTION:
 *
 *
 *   REFERENCE:  Sub-clauses 4.1.9, 4.1.8.2 of GSM Recommendation 06.20
 *
 *   KEYWORDS: pitchLags, pitchlags, PITCH_LAGS
 *
 *************************************************************************/

void   pitchLags(Shortword swBestIntLag,
                        Shortword pswIntCs[],
                        Shortword pswIntGs[],
                        Shortword swCCThreshold,
                        Shortword pswLPeaksSorted[],
                        Shortword pswCPeaksSorted[],
                        Shortword pswGPeaksSorted[],
                        Shortword *psiNumSorted,
                        Shortword *pswPitch,
                        Shortword *pswHNWCoef)
{

/*_________________________________________________________________________
 |                                                                         |
 |                            Automatic Variables                          |
 |_________________________________________________________________________|
*/

  Shortword pswLBuf[2 * OS_FCTR - 1],
         pswCBuf[2 * OS_FCTR - 1];
  Shortword pswGBuf[2 * OS_FCTR - 1],
         pswLPeaks[2 * LMAX / LMIN];
  Shortword pswCPeaks[2 * LMAX / LMIN],
         pswGPeaks[2 * LMAX / LMIN];
  short  siLPeakIndex,
         siCPeakIndex,
         siGPeakIndex,
         siPeakIndex;
  short  siSortedIndex,
         siLPeaksSortedInit,
         swWorkingLag;
  Shortword swSubMult,
         swFullResPeak,
         swCPitch,
         swGPitch,
         swMult;
  Shortword swMultInt,
         sw1,
         sw2,
         si1,
         si2;
  Longword L_1;
  short  siNum,
         siUpperBound,
         siLowerBound,
         siSMFIndex;
  short  siNumPeaks,
         siRepeat,
         i,
         j;

  static ShortwordRom psrSubMultFactor[] = {0x0aab,     /* 1.0/12.0 */
    0x071c,                            /* 1.0/18.0 */
    0x0555,                            /* 1.0/24.0 */
    0x0444,                            /* 1.0/30.0 */
  0x038e};                             /* 1.0/36.0 */


/*_________________________________________________________________________
 |                                                                         |
 |                              Executable Code                            |
 |_________________________________________________________________________|
*/

  /* Get array of valid lags within one integer lag of the best open-loop */
  /* integer lag; get the corresponding interpolated C and G arrays;      */
  /* find the best among these; store the info corresponding to this peak */
  /* in the interpolated CC/G sequence                                    */
  /*----------------------------------------------------------------------*/

  sw1 = shr(extract_l(L_mult(swBestIntLag, OS_FCTR)), 1);

  siNum = CGInterpValid(sw1, pswIntCs, pswIntGs,
                        pswLBuf, pswCBuf, pswGBuf);

  sw1 = 0;
  sw2 = 0x003f;
  siPeakIndex = fnBest_CG(pswCBuf, pswGBuf, &sw1, &sw2, siNum);
  if (siPeakIndex == -1)
  {

    /* It is possible, although rare, that the interpolated sequence */
    /* will not have a peak where the original sequence did.         */
    /* Indicate this on return                                       */
    /*---------------------------------------------------------------*/

    *psiNumSorted = 0;
    return;
  }

  siLPeakIndex = 0;
  siCPeakIndex = 0;
  siGPeakIndex = 0;
  siSortedIndex = 0;
  pswCPeaks[siCPeakIndex] = pswCBuf[siPeakIndex];
  siCPeakIndex = add(siCPeakIndex, 1);
  pswLPeaks[siLPeakIndex] = pswLBuf[siPeakIndex];
  siLPeakIndex = add(siLPeakIndex, 1);
  pswGPeaks[siGPeakIndex] = pswGBuf[siPeakIndex];
  siGPeakIndex = add(siGPeakIndex, 1);

  /* Find all peaks at submultiples of the first peak */
  /*--------------------------------------------------*/

  siSMFIndex = 0;
  swSubMult = mult_r(pswLPeaks[0], psrSubMultFactor[siSMFIndex++]);

  while (sub(swSubMult, LMIN) >= 0 && siSMFIndex <= 5)
  {

    /* Check if there is peak in the integer CC/G sequence within */
    /* PEAK_VICINITY of the submultiple                           */
    /*------------------------------------------------------------*/

    swFullResPeak = findPeak(swSubMult, pswIntCs, pswIntGs);

    if (swFullResPeak)
    {

      /* Peak found at submultiple: interpolate to get C's and G's */
      /* corresponding to valid lags within one of the new found   */
      /* peak; get best C**2/G from these;  if it meets threshold, */
      /* store the info corresponding to this peak                 */
      /*-----------------------------------------------------------*/

      siNum = CGInterpValid(swFullResPeak, pswIntCs, pswIntGs,
                            pswLBuf, pswCBuf, pswGBuf);

      sw1 = swCCThreshold;
      sw2 = 0x4000;
      siPeakIndex = fnBest_CG(pswCBuf, pswGBuf, &sw1, &sw2, siNum);
      if (siPeakIndex != -1)
      {

        pswCPeaks[siCPeakIndex] = pswCBuf[siPeakIndex];
        siCPeakIndex = add(siCPeakIndex, 1);
        pswLPeaks[siLPeakIndex] = pswLBuf[siPeakIndex];
        siLPeakIndex = add(siLPeakIndex, 1);
        pswGPeaks[siGPeakIndex] = pswGBuf[siPeakIndex];
        siGPeakIndex = add(siGPeakIndex, 1);

      }
    }


    if (siSMFIndex < 5)
    {

      /* Get next submultiple */
      /*----------------------*/
      swSubMult = mult_r(pswLPeaks[0], psrSubMultFactor[siSMFIndex]);

    }

    siSMFIndex++;
  }

  /* Get pitch from fundamental peak: first, build array of fractional */
  /* lags around which to search for peak.  Note that these lags are   */
  /* NOT restricted to those in the lag table, but may take any value  */
  /* in the range LMIN_FR to LMAX_FR                                   */
  /*-------------------------------------------------------------------*/

  siUpperBound = add(pswLPeaks[siLPeakIndex - 1], OS_FCTR);
  siUpperBound = sub(siUpperBound, 1);
  if (sub(siUpperBound, LMAX_FR) > 0)
    siUpperBound = LMAX_FR;

  siLowerBound = sub(pswLPeaks[siLPeakIndex - 1], OS_FCTR);
  siLowerBound = add(siLowerBound, 1);
  if (sub(siLowerBound, LMIN_FR) < 0)
    siLowerBound = LMIN_FR;
  for (i = siLowerBound, j = 0; i <= siUpperBound; i++, j++)
    pswLBuf[j] = i;

  /* Second, find peak in the interpolated CC/G sequence. */
  /* The corresponding lag is the fundamental pitch.  The */
  /* interpolated C(pitch) and G(pitch) values are stored */
  /* for later use in calculating the Harmonic-Noise-     */
  /* Weighting coefficient                                */
  /*------------------------------------------------------*/

  siNum = sub(siUpperBound, siLowerBound);
  siNum = add(siNum, 1);
  CGInterp(pswLBuf, siNum, pswIntCs, pswIntGs, LSMIN,
           pswCBuf, pswGBuf);
  sw1 = 0;
  sw2 = 0x003f;
  siPeakIndex = fnBest_CG(pswCBuf, pswGBuf, &sw1, &sw2, siNum);
  if (siPeakIndex == -1)
  {
    swCPitch = 0;
    *pswPitch = LMIN_FR;
    swGPitch = 0x003f;
  }
  else
  {
    swCPitch = pswCBuf[siPeakIndex];
    *pswPitch = pswLBuf[siPeakIndex];
    swGPitch = pswGBuf[siPeakIndex];
  }

  /* Find all peaks which are multiples of fundamental pitch */
  /*---------------------------------------------------------*/

  swMult = add(*pswPitch, *pswPitch);
  swMultInt = mult_r(swMult, INV_OS_FCTR);

  while (sub(swMultInt, LMAX) <= 0)
  {

    /* Check if there is peak in the integer CC/G sequence within */
    /* PEAK_VICINITY of the multiple                              */
    /*------------------------------------------------------------*/

    swFullResPeak = findPeak(swMultInt, pswIntCs, pswIntGs);

    if (swFullResPeak)
    {

      /* Peak found at multiple: interpolate to get C's and G's    */
      /* corresponding to valid lags within one of the new found   */
      /* peak; get best C**2/G from these;  if it meets threshold, */
      /* store the info corresponding to this peak                 */
      /*-----------------------------------------------------------*/

      siNum = CGInterpValid(swFullResPeak, pswIntCs, pswIntGs,
                            pswLBuf, pswCBuf, pswGBuf);

      sw1 = swCCThreshold;
      sw2 = 0x4000;
      siPeakIndex = fnBest_CG(pswCBuf, pswGBuf, &sw1, &sw2, siNum);
      if (siPeakIndex != -1)
      {

        pswCPeaks[siCPeakIndex] = pswCBuf[siPeakIndex];
        siCPeakIndex = add(siCPeakIndex, 1);
        pswLPeaks[siLPeakIndex] = pswLBuf[siPeakIndex];
        siLPeakIndex = add(siLPeakIndex, 1);
        pswGPeaks[siGPeakIndex] = pswGBuf[siPeakIndex];
        siGPeakIndex = add(siGPeakIndex, 1);

      }
    }

    /* Get next multiple */
    /*-------------------*/

    swMult = add(*pswPitch, swMult);
    swMultInt = mult_r(swMult, INV_OS_FCTR);
  }

  /* Get Harmonic-Noise-Weighting coefficient = 0.4 * C(pitch) / G(pitch) */
  /* Note: a factor of 0.5 is applied the the HNW coeffcient              */
  /*----------------------------------------------------------------------*/

  si2 = norm_s(swCPitch);
  sw1 = shl(swCPitch, si2);
  L_1 = L_mult(sw1, PW_FRAC);

  si1 = norm_s(swGPitch);
  sw1 = shl(swGPitch, si1);

  sw2 = round(L_shr(L_1, 1));
  sw2 = divide_s(sw2, sw1);

  if (sub(si1, si2) > 0)
    sw2 = shl(sw2, sub(si1, si2));

  if (sub(si1, si2) < 0)
    sw2 = shift_r(sw2, sub(si1, si2));

  *pswHNWCoef = sw2;

  /* Sort peaks into output arrays, largest first */
  /*----------------------------------------------*/

  siLPeaksSortedInit = siSortedIndex;
  *psiNumSorted = 0;
  siNumPeaks = siLPeakIndex;
  for (i = 0; i < siNumPeaks; i++)
  {

    sw1 = 0;
    sw2 = 0x003f;
    siPeakIndex = fnBest_CG(pswCPeaks, pswGPeaks, &sw1, &sw2, siNumPeaks);

    siRepeat = 0;
    swWorkingLag = pswLPeaks[siPeakIndex];
    for (j = 0; j < *psiNumSorted; j++)
    {

      if (sub(swWorkingLag, pswLPeaksSorted[j + siLPeaksSortedInit]) == 0)
        siRepeat = 1;
    }

    if (!siRepeat)
    {

      pswLPeaksSorted[siSortedIndex] = swWorkingLag;
      pswCPeaksSorted[siSortedIndex] = pswCPeaks[siPeakIndex];
      pswGPeaksSorted[siSortedIndex] = pswGPeaks[siPeakIndex];
      siSortedIndex = add(siSortedIndex, 1);
      *psiNumSorted = add(*psiNumSorted, 1);
    }

    pswCPeaks[siPeakIndex] = 0x0;
  }
}                                      /* end of pitchLags */



/***************************************************************************
 *
 *   FUNCTION NAME: quantLag
 *
 *   PURPOSE:
 *
 *     Quantizes a given fractional lag according to the provided table
 *     of allowable fractional lags.
 *
 *   INPUTS:
 *
 *     swRawLag
 *
 *                     Raw lag value: a fractional lag value*OS_FCTR.
 *
 *     psrLagTbl[0:255]
 *
 *                     Fractional lag table.
 *
 *   OUTPUTS:
 *
 *     pswCode
 *
 *                     Index in lag table of the quantized lag-- that is,
 *                     the coded value of the lag.
 *
 *   RETURN VALUE:
 *
 *     Quantized lag value.
 *
 *
 *   REFERENCE:  Sub-clause 4.1.8.3 of GSM Recommendation 06.20
 *
 *   KEYWORDS: quantization, LTP, adaptive codebook
 *
 *************************************************************************/

Shortword quantLag(Shortword swRawLag,
                          Shortword *pswCode)
{

/*_________________________________________________________________________
 |                                                                         |
 |                            Automatic Variables                          |
 |_________________________________________________________________________|
*/
  Shortword siOffset,
         swIndex1,
         swIndex2;

/*_________________________________________________________________________
 |                                                                         |
 |                              Executable Code                            |
 |_________________________________________________________________________|
*/

  swIndex1 = 0;
  siOffset = shr(LAG_TABLE_LEN, 1);
  swIndex2 = siOffset;
  siOffset = shr(siOffset, 1);

  while (sub(swIndex2, swIndex1) != 0)
  {
    if (sub(swRawLag, psrLagTbl[swIndex2]) >= 0)
      swIndex1 = swIndex2;
    swIndex2 = add(swIndex1, siOffset);
    siOffset = shr(siOffset, 1);
  }
  *pswCode = swIndex2;

  return (psrLagTbl[swIndex2]);

}                                      /* end of quantLag */

/***************************************************************************
 *
 *   FUNCTION NAME: r0Quant
 *
 *   PURPOSE:
 *
 *      Quantize the unquantized R(0).  Returns r0 codeword (index).
 *
 *   INPUTS:
 *
 *     L_UnqntzdR0
 *                     The average frame energy R(0)
 *
 *   OUTPUTS: none
 *
 *   RETURN VALUE:
 *
 *                     A 16 bit number representing the codeword whose
 *                     associated R(0) is closest to the input frame energy.
 *
 *   DESCRIPTION:
 *
 *     Returns r0 codeword (index) not the actual Rq(0).
 *
 *     Level compare input with boundary value (the boundary
 *     above ,louder) of candidate r0Index i.e.
 *     lowerBnd[i] <= inputR(0) < upperBnd[i+1]
 *
 *     The compare in the routine is:
 *     inputR(0) < upperBnd[i+1] if false return i as codeword
 *
 *   REFERENCE:  Sub-clause 4.1.3 of GSM Recommendation 06.20
 *
 *
 *************************************************************************/

Shortword r0Quant(Longword L_UnqntzdR0)
{

/*_________________________________________________________________________
 |                                                                         |
 |                            Automatic Variables                          |
 |_________________________________________________________________________|
*/

  Shortword swR0Index,
         swUnqntzdR0;

/*_________________________________________________________________________
 |                                                                         |
 |                            Executable Code                              |
 |_________________________________________________________________________|
*/

  swUnqntzdR0 = round(L_UnqntzdR0);

  for (swR0Index = 0; swR0Index < (1 << R0BITS) - 1; swR0Index++)
  {
    if (sub(swUnqntzdR0, psrR0DecTbl[2 * swR0Index + 1]) < 0)
    {
      return (swR0Index);
    }
  }
  return ((1 << R0BITS) - 1);          /* return the maximum */
}

/***************************************************************************
 *
 *   FUNCTION NAME: setupPreQ
 *
 *   PURPOSE:
 *     The purpose of this function is to setup the internal pointers so that
 *     getNextVec knows where to start.
 *
 *   INPUTS: iSeg, iVector
 *
 *   OUTPUTS: None
 *
 *   RETURN VALUE: none
 *
 *   DESCRIPTION:
 *
 *   REFERENCE:  Sub-clause  4.1.4.1 of GSM Recommendation 06.20
 *
 *   KEYWORDS:  vector quantizer preQ
 *
 *************************************************************************/

void   setupPreQ(int iSeg, int iVector)
{

/*_________________________________________________________________________
 |                                                                         |
 |                              Executable Code                            |
 |_________________________________________________________________________|
*/

  iLow = psvqIndex[iSeg - 1].l;
  iThree = ((psvqIndex[iSeg - 1].h - iLow) == 2);

  switch (iSeg)
  {
    case 1:
      {
        psrTable = psrPreQ1;
        iWordPtr = (iVector * 3) >> 1;
        if (odd(iVector))
          iWordHalfPtr = LOW;
        else
          iWordHalfPtr = HIGH;
        break;
      }

    case 2:
      {
        psrTable = psrPreQ2;
        iWordPtr = (iVector * 3) >> 1;
        if (odd(iVector))
          iWordHalfPtr = LOW;
        else
          iWordHalfPtr = HIGH;
        break;
      }

    case 3:
      {
        psrTable = psrPreQ3;
        iWordPtr = iVector * 2;
        iWordHalfPtr = HIGH;
        break;
      }
  }
}

/***************************************************************************
 *
 *   FUNCTION NAME: setupQuant
 *
 *   PURPOSE:
 *     The purpose of this function is to setup the internal pointers so that
 *     getNextVec knows where to start.
 *
 *
 *   INPUTS: iSeg, iVector
 *
 *   OUTPUTS: None
 *
 *   RETURN VALUE: none
 *
 *   DESCRIPTION:
 *
 *   REFERENCE:  Sub-clause 4.1.4.1 of GSM Recommendation 06.20
 *
 *   KEYWORDS:  vector quantizer Quant
 *
 *************************************************************************/

void   setupQuant(int iSeg, int iVector)
{

/*_________________________________________________________________________
 |                                                                         |
 |                              Executable Code                            |
 |_________________________________________________________________________|
*/

  iLow = psvqIndex[iSeg - 1].l;
  iThree = ((psvqIndex[iSeg - 1].h - iLow) == 2);

  switch (iSeg)
  {
    case 1:
      {
        psrTable = psrQuant1;
        iWordPtr = (iVector * 3) >> 1;
        if (odd(iVector))
          iWordHalfPtr = LOW;
        else
          iWordHalfPtr = HIGH;
        break;
      }

    case 2:
      {
        psrTable = psrQuant2;
        iWordPtr = (iVector * 3) >> 1;
        if (odd(iVector))
          iWordHalfPtr = LOW;
        else
          iWordHalfPtr = HIGH;
        break;
      }

    case 3:
      {
        psrTable = psrQuant3;
        iWordPtr = iVector * 2;
        iWordHalfPtr = HIGH;
        break;
      }
  }
}

/***************************************************************************
 *
 *   FUNCTION NAME: weightSpeechFrame
 *
 *   PURPOSE:
 *
 *     The purpose of this function is to perform the spectral
 *     weighting filter  (W(z)) of the input speech frame.
 *
 *   INPUTS:
 *
 *     pswSpeechFrm[0:F_LEN]
 *
 *                     high pass filtered input speech frame
 *
 *     pswWNumSpace[0:NP*N_SUB]
 *
 *                     W(z) numerator coefficients
 *
 *     pswWDenomSpace[0:NP*N_SUB]
 *
 *                     W(z) denominator coefficients
 *
 *     pswWSpeechBuffBase[0:F_LEN+LMAX+CG_INT_MACS/2]
 *
 *                     previous W(z) filtered speech
 *
 *   OUTPUTS:
 *
 *     pswWSpeechBuffBase[0:F_LEN+LMAX+CG_INT_MACS/2]
 *
 *                     W(z) filtered speech frame
 *
 *   RETURN VALUE:
 *
 *     none
 *
 *   DESCRIPTION:
 *
 *   REFERENCE:  Sub-clause 4.1.7 of GSM Recommendation 06.20
 *
 *   KEYWORDS:
 *
 *************************************************************************/

void   weightSpeechFrame(Shortword pswSpeechFrm[],
                                Shortword pswWNumSpace[],
                                Shortword pswWDenomSpace[],
                                Shortword pswWSpeechBuffBase[])
{

/*_________________________________________________________________________
 |                                                                         |
 |                            Automatic Variables                          |
 |_________________________________________________________________________|
*/

  Shortword pswScratch0[W_F_BUFF_LEN],
        *pswWSpeechFrm;
  short int siI;

/*_________________________________________________________________________
 |                                                                         |
 |                              Executable Code                            |
 |_________________________________________________________________________|
*/

  /* Delay the weighted speech buffer by one frame */
  /* --------------------------------------------- */

  for (siI = 0; siI < LSMAX; siI++)
  {
    pswWSpeechBuffBase[siI] = pswWSpeechBuffBase[siI + F_LEN];
  }

  /* pass speech frame through W(z) */
  /* ------------------------------ */

  pswWSpeechFrm = pswWSpeechBuffBase + LSMAX;

  for (siI = 0; siI < N_SUB; siI++)
  {
    lpcFir(&pswSpeechFrm[siI * S_LEN], &pswWNumSpace[siI * NP],
           pswWStateNum, &pswScratch0[siI * S_LEN]);
  }

  for (siI = 0; siI < N_SUB; siI++)
  {
    lpcIir(&pswScratch0[siI * S_LEN], &pswWDenomSpace[siI * NP],
           pswWStateDenom, &pswWSpeechFrm[siI * S_LEN]);
  }
}
