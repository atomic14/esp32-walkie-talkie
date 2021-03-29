/***************************************************************************
 *
 *   File Name:  sp_dec.c
 *
 *   Purpose:
 *      Contains all functions for decoding speech.  It does not
 *      include those routines needed to decode channel information.
 *
 *      Since the GSM half-rate speech coder is an analysis-by-synthesis
 *      coder, many of the routines in this file are also called by the
 *      encoder.  Functions are included for coded-parameter lookup,
 *      LPC filter coefficient interpolation, excitation vector lookup
 *      and construction, vector quantized gain lookup, and LPC synthesis
 *      filtering.  In addition, some post-processing functions are
 *      included.
 *
 *     Below is a listing of all the functions appearing in the file.
 *     The functions are arranged according to their purpose.  Under
 *     each heading, the ordering is hierarchical.
 *
 *     The entire speech decoder, under which all these routines fall,
 *     except were noted:
 *     speechDecoder()
 *
 *     Spectral Smoothing of LPC:
 *       a_sst()
 *         aFlatRcDp()
 *         rcToCorrDpL()
 *         aToRc()
 *         rcToADp()
 *     VSELP codevector construction:
 *       b_con()
 *       v_con()
 *     LTP vector contruction:
 *       fp_ex()
 *         get_ipjj()
 *       lagDecode()
 *     LPC contruction
 *       getSfrmLpc()
 *         interpolateCheck()
 *         res_eng()
 *       lookupVq()
 *     Excitation scaling:
 *       rs_rr()
 *         g_corr1() (no scaling)
 *       rs_rrNs()
 *         g_corr1s() (g_corr1 with scaling)
 *       scaleExcite()
 *     Post filtering:
 *       pitchPreFilt()
 *         agcGain()
 *         lpcIir()
 *       r0BasedEnergyShft()
 *       spectralPostFilter()
 *         lpcFir()
 *
 *
 *     Routines not referenced by speechDecoder()
 *     Filtering routines:
 *       lpcIrZsIir()
 *       lpcZiIir()
 *       lpcZsFir()
 *       lpcZsIir()
 *       lpcZsIirP()
 *     Square root:
 *       sqroot()
 *
 **************************************************************************/

/*_________________________________________________________________________
 |                                                                         |
 |                            Include Files                                |
 |_________________________________________________________________________|
*/

#include "typedefs.h"
#include "mathhalf.h"
#include "sp_rom.h"
#include "sp_dec.h"
#include "err_conc.h"
#include "dtx.h"


/*_________________________________________________________________________
 |                                                                         |
 |            Local Functions (scope is limited to this file)              |
 |_________________________________________________________________________|
*/

static void a_sst(Shortword swAshift, Shortword swAscale,
                         Shortword pswDirectFormCoefIn[],
                         Shortword pswDirectFormCoefOut[]);

  static short aToRc(Shortword swAshift, Shortword pswAin[],
                            Shortword pswRc[]);

  static Shortword agcGain(Shortword pswStateCurr[],
                                  struct NormSw snsInSigEnergy,
                                  Shortword swEngyRShft);

  static Shortword lagDecode(Shortword swDeltaLag);

  static void lookupVq(Shortword pswVqCodeWds[], Shortword pswRCOut[]);

  static void pitchPreFilt(Shortword pswExcite[],
                                  Shortword swRxGsp0,
                                  Shortword swRxLag,
                                  Shortword swUvCode,
                                  Shortword swSemiBeta,
                                  struct NormSw snsSqrtRs,
                                  Shortword pswExciteOut[],
                                  Shortword pswPPreState[]);

  static void spectralPostFilter(Shortword pswSPFIn[],
                           Shortword pswNumCoef[], Shortword pswDenomCoef[],
                                        Shortword pswSPFOut[]);

/*_________________________________________________________________________
 |                                                                         |
 |                              Local Defines                              |
 |_________________________________________________________________________|
*/

#define  P_INT_MACS   10
#define  ASCALE       0x0800
#define  ASHIFT       4
#define  DELTA_LEVELS 16
#define  GSP0_SCALE   1
#define  C_BITS_V     9                /* number of bits in any voiced VSELP
                                        * codeword */
#define  C_BITS_UV    7                /* number of bits in a unvoiced VSELP
                                        * codeword */
#define  MAXBITS      C_BITS_V         /* max number of bits in any VSELP
                                        * codeword */
#define  LTP_LEN      147              /* 147==0x93 length of LTP history */
#define  SQRT_ONEHALF 0x5a82           /* the 0.5 ** 0.5 */
#define  LPC_ROUND    0x00000800L      /* 0x8000 >> ASHIFT */
#define  AFSHIFT      2                /* number of right shifts to be
                                        * applied to the autocorrelation
                                        * sequence in aFlatRcDp     */

/*_________________________________________________________________________
 |                                                                         |
 |                         State variables (globals)                       |
 |_________________________________________________________________________|
*/

  Shortword gswPostFiltAgcGain,
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


  Shortword swMuteFlagOld;             /* error concealment */


 /* DTX state variables */
 /* ------------------- */

  Shortword swRxDTXState = CNINTPER - 1;        /* DTX State at the rx.
                                                 * Modulo */

 /* counter [0,11].             */

  Shortword swDecoMode = SPEECH;
  Shortword swDtxMuting = 0;
  Shortword swDtxBfiCnt = 0;

  Shortword swOldR0IndexDec = 0;

  Shortword swRxGsHistPtr = 0;
  Longword pL_RxGsHist[(OVERHANG - 1) * N_SUB];


/*_________________________________________________________________________
 |                                                                         |
 |                               Global Data                               |
 |                     (scope is global to this file)                      |
 |_________________________________________________________________________|
*/

  Shortword swR0Dec;

  Shortword swVoicingMode,             /* MODE */
         pswVq[3],                     /* LPC1, LPC2, LPC3 */
         swSi,                         /* INT_LPC */
         swEngyRShift;                 /* for use by spectral postfilter */


  Shortword swR0NewCN;                 /* DTX mode */

  extern LongwordRom ppLr_gsTable[4][32];       /* DTX mode */


/***************************************************************************
 *
 *   FUNCTION NAME: aFlatRcDp
 *
 *   PURPOSE:
 *
 *     Given a Longword autocorrelation sequence, representing LPC
 *     information, aFlatRcDp converts the vector to one of NP
 *     Shortword reflection coefficients.
 *
 *   INPUT:
 *
 *
 *     pL_R[0:NP]    - An input Longword autocorrelation sequence, (pL_R[0] =
 *                     not necessarily 0x7fffffffL).  pL_R is altered in the
 *                     call, by being right shifted by global constant
 *                     AFSHIFT bits.
 *
 *                     The input array pL_R[] should be shifted left as much
 *                     as possible to improve precision.
 *
 *     AFSHIFT       - The number of right shifts to be applied to the
 *                     normalized autocorrelation sequence pL_R.
 *
 *   OUTPUT:
 *
 *     pswRc[0:NP-1] - A Shortword output vector of NP reflection
 *                     coefficients.
 *
 *   RETURN VALUE:
 *
 *     None
 *
 *   DESCRIPTION:
 *
 *     This routine transforms LPC information from one set of
 *     parameters to another.  It is better suited for fixed point
 *     implementations than the Levinson-Dubin recursion.
 *
 *     The function is called by a_sst(), and getNWCoefs().  In a_sst()
 *     direct form coefficients are converted to autocorrelations,
 *     and smoothed in that domain.  Conversion back to direct form
 *     coefficients is done by calling aFlatRc(), followed by rcToADp().
 *
 *     In getNwCoefs() again the conversion back to direct form
 *     coefficients is done by calling aFlatRc(), followed by rcToADp().
 *     In getNwCoefs() an autocorrelation sequence is generated from the
 *     impulse response of the weighting filters.
 *
 *     The fundamental recursion is derived from AFLAT, which is
 *     described in section 4.1.4.1.
 *
 *     Unlike in AFLAT where the reflection coefficients are known, here
 *     they are the unknowns.  PBar and VBar for j==0 are initially
 *     known, as is rSub1.  From this information the next set of P's
 *     and V's are generated.  At the end of the recursion the next,
 *     reflection coefficient rSubj (pswRc[j]) can be calcluated by
 *     dividing Vsubj by Psubj.
 *
 *     Precision is crucial in this routine.  At each stage, a
 *     normalization is performed prior to the reflection coefficient
 *     calculation.  In addition, to prevent overflow, the
 *     autocorrelation sequence is scaled down by ASHIFT (4) right
 *     shifts.
 *
 *
 *   REFERENCES: Sub_Clause 4.1.9 and 4.2.1  of GSM Recomendation 06.20
 *
 *   KEYWORDS: reflection coefficients, AFLAT, aflat, recursion, LPC
 *   KEYWORDS: autocorrelation
 *
 *************************************************************************/

  void   aFlatRcDp(Longword *pL_R, Shortword *pswRc)
{

/*_________________________________________________________________________
 |                                                                         |
 |                            Automatic Variables                          |
 |_________________________________________________________________________|
*/

  Longword pL_pjNewSpace[NP];
  Longword pL_pjOldSpace[NP];
  Longword pL_vjNewSpace[2 * NP - 1];
  Longword pL_vjOldSpace[2 * NP - 1];

  Longword *pL_pjOld;
  Longword *pL_pjNew;
  Longword *pL_vjOld;
  Longword *pL_vjNew;
  Longword *pL_swap;

  Longword L_temp;
  Longword L_sum;
  Shortword swRc,
         swRcSq,
         swTemp,
         swTemp1,
         swAbsTemp1,
         swTemp2;
  int    i,
         j;


/*_________________________________________________________________________
 |                                                                         |
 |                              Executable Code                            |
 |_________________________________________________________________________|
*/

  pL_pjOld = pL_pjOldSpace;
  pL_pjNew = pL_pjNewSpace;
  pL_vjOld = pL_vjOldSpace + NP - 1;
  pL_vjNew = pL_vjNewSpace + NP - 1;


  /* Extract the 0-th reflection coefficient */
  /*-----------------------------------------*/

  swTemp1 = round(pL_R[1]);
  swTemp2 = round(pL_R[0]);
  swAbsTemp1 = abs_s(swTemp1);
  if (swTemp2 <= 0 || sub(swAbsTemp1, swTemp2) >= 0)
  {
    j = 0;
    for (i = j; i < NP; i++)
    {
      pswRc[i] = 0;
    }
    return;
  }

  swRc = divide_s(swAbsTemp1, swTemp2);/* return division result */

  if (sub(swTemp1, swAbsTemp1) == 0)
    swRc = negate(swRc);               /* negate reflection Rc[j] */

  pswRc[0] = swRc;                     /* copy into the output Rc array */

  for (i = 0; i <= NP; i++)
  {
    pL_R[i] = L_shr(pL_R[i], AFSHIFT);
  }

  /* Initialize the pjOld and vjOld recursion arrays */
  /*-------------------------------------------------*/

  for (i = 0; i < NP; i++)
  {
    pL_pjOld[i] = pL_R[i];
    pL_vjOld[i] = pL_R[i + 1];
  }
  for (i = -1; i > -NP; i--)
    pL_vjOld[i] = pL_R[-(i + 1)];


  /* Compute the square of the j=0 reflection coefficient */
  /*------------------------------------------------------*/

  swRcSq = mult_r(swRc, swRc);

  /* Update pjNew and vjNew arrays for lattice stage j=1 */
  /*-----------------------------------------------------*/

  /* Updating pjNew: */
  /*-------------------*/

  for (i = 0; i <= NP - 2; i++)
  {
    L_temp = L_mpy_ls(pL_vjOld[i], swRc);
    L_sum = L_add(L_temp, pL_pjOld[i]);
    L_temp = L_mpy_ls(pL_pjOld[i], swRcSq);
    L_sum = L_add(L_temp, L_sum);
    L_temp = L_mpy_ls(pL_vjOld[-i], swRc);
    pL_pjNew[i] = L_add(L_sum, L_temp);
  }

  /* Updating vjNew: */
  /*-------------------*/

  for (i = -NP + 2; i <= NP - 2; i++)
  {
    L_temp = L_mpy_ls(pL_vjOld[-i - 1], swRcSq);
    L_sum = L_add(L_temp, pL_vjOld[i + 1]);
    L_temp = L_mpy_ls(pL_pjOld[(((i + 1) >= 0) ? i + 1 : -(i + 1))], swRc);
    L_temp = L_shl(L_temp, 1);
    pL_vjNew[i] = L_add(L_temp, L_sum);
  }



  j = 0;

  /* Compute reflection coefficients Rc[1],...,Rc[9] */
  /*-------------------------------------------------*/

  for (j = 1; j < NP; j++)
  {

    /* Swap pjNew and pjOld buffers */
    /*------------------------------*/

    pL_swap = pL_pjNew;
    pL_pjNew = pL_pjOld;
    pL_pjOld = pL_swap;

    /* Swap vjNew and vjOld buffers */
    /*------------------------------*/

    pL_swap = pL_vjNew;
    pL_vjNew = pL_vjOld;
    pL_vjOld = pL_swap;

    /* Compute the j-th reflection coefficient */
    /*-----------------------------------------*/

    swTemp = norm_l(pL_pjOld[0]);      /* get shift count */
    swTemp1 = round(L_shl(pL_vjOld[0], swTemp));        /* normalize num.  */
    swTemp2 = round(L_shl(pL_pjOld[0], swTemp));        /* normalize den.  */

    /* Test for invalid divide conditions: a) devisor < 0 b) abs(divident) >
     * abs(devisor) If either of these conditions is true, zero out
     * reflection coefficients for i=j,...,NP-1 and return. */

    swAbsTemp1 = abs_s(swTemp1);
    if (swTemp2 <= 0 || sub(swAbsTemp1, swTemp2) >= 0)
    {
      i = j;
      for (i = j; i < NP; i++)
      {
        pswRc[i] = 0;
      }
      return;
    }

    swRc = divide_s(swAbsTemp1, swTemp2);       /* return division result */
    if (sub(swTemp1, swAbsTemp1) == 0)
      swRc = negate(swRc);             /* negate reflection Rc[j] */
    swRcSq = mult_r(swRc, swRc);       /* compute Rc^2 */
    pswRc[j] = swRc;                   /* copy Rc[j] to output array */

    /* Update pjNew and vjNew arrays for the next lattice stage if j < NP-1 */
    /*---------------------------------------------------------------------*/

    /* Updating pjNew: */
    /*-----------------*/

    for (i = 0; i <= NP - j - 2; i++)
    {
      L_temp = L_mpy_ls(pL_vjOld[i], swRc);
      L_sum = L_add(L_temp, pL_pjOld[i]);
      L_temp = L_mpy_ls(pL_pjOld[i], swRcSq);
      L_sum = L_add(L_temp, L_sum);
      L_temp = L_mpy_ls(pL_vjOld[-i], swRc);
      pL_pjNew[i] = L_add(L_sum, L_temp);
    }

    /* Updating vjNew: */
    /*-----------------*/

    for (i = -NP + j + 2; i <= NP - j - 2; i++)
    {
      L_temp = L_mpy_ls(pL_vjOld[-i - 1], swRcSq);
      L_sum = L_add(L_temp, pL_vjOld[i + 1]);
      L_temp = L_mpy_ls(pL_pjOld[(((i + 1) >= 0) ? i + 1 : -(i + 1))], swRc);
      L_temp = L_shl(L_temp, 1);
      pL_vjNew[i] = L_add(L_temp, L_sum);
    }
  }
  return;
}

/***************************************************************************
 *
 *   FUNCTION NAME: aToRc
 *
 *   PURPOSE:
 *
 *     This subroutine computes a vector of reflection coefficients, given
 *     an input vector of direct form LPC filter coefficients.
 *
 *   INPUTS:
 *
 *     NP
 *                     order of the LPC filter (global constant)
 *
 *     swAshift
 *                     The number of right shifts applied externally
 *                     to the direct form filter coefficients.
 *
 *     pswAin[0:NP-1]
 *                     The input vector of direct form LPC filter
 *                     coefficients.
 *
 *   OUTPUTS:
 *
 *     pswRc[0:NP-1]
 *                     Array containing the reflection coefficients.
 *
 *   RETURN VALUE:
 *
 *     siUnstableFlt
 *                     If stable reflection coefficients 0, 1 if unstable.
 *
 *
 *   DESCRIPTION:
 *
 *     This function performs the conversion from direct form
 *     coefficients to reflection coefficients. It is used in a_sst()
 *     and interpolateCheck().  In a_sst() reflection coefficients used
 *     as a transitional data format.  aToRc() is used for this
 *     conversion.
 *
 *     When performing interpolation, a stability check must be
 *     performed. interpolateCheck() does this by calling aToRc().
 *
 *     First coefficients are shifted down by iAshift. NP, the filter
 *     order is 10. The a's and rc's each have NP elements in them. An
 *     elaborate algorithm description can be found on page 443, of
 *     "Digital Processing of Speech Signals" by L.R. Rabiner and R.W.
 *     Schafer; Prentice-Hall; Englewood Cliffs, NJ (USA).  1978.
 *
 *   REFERENCES: Sub_Clause 4.1.6, and 4.2.3 of GSM Recomendation 06.20
 *
 *   KEYWORDS: reflectioncoefficients, parcors, conversion, atorc, ks, as
 *   KEYWORDS: parcorcoefficients, lpc, flat, vectorquantization
 *
 *************************************************************************/

static short aToRc(Shortword swAshift, Shortword pswAin[],
                          Shortword pswRc[])
{

/*_________________________________________________________________________
 |                                                                         |
 |                            Constants                                    |
 |_________________________________________________________________________|
*/

/*_________________________________________________________________________
 |                                                                         |
 |                            Automatic Variables                          |
 |_________________________________________________________________________|
*/

  Shortword pswTmpSpace[NP],
         pswASpace[NP],
         swNormShift,
         swActShift,
         swNormProd,
         swRcOverE,
         swDiv,
        *pswSwap,
        *pswTmp,
        *pswA;

  Longword L_temp;

  short int siUnstableFlt,
         i,
         j;                            /* Loop control variables */

/*_________________________________________________________________________
 |                                                                         |
 |                              Executable Code                            |
 |_________________________________________________________________________|
*/

  /* Initialize starting addresses for temporary buffers */
  /*-----------------------------------------------------*/

  pswA = pswASpace;
  pswTmp = pswTmpSpace;

  /* Copy the direct form filter coefficients to a temporary array */
  /*---------------------------------------------------------------*/

  for (i = 0; i < NP; i++)
  {
    pswA[i] = pswAin[i];
  }

  /* Initialize the flag for filter stability check */
  /*------------------------------------------------*/

  siUnstableFlt = 0;

  /* Start computation of the reflection coefficients, Rc[9],...,Rc[1] */
  /*-------------------------------------------------------------------*/

  for (i = NP - 1; i >= 1; i--)
  {

    pswRc[i] = shl(pswA[i], swAshift); /* write Rc[i] to output array */

    /* Check the stability of i-th reflection coefficient */
    /*----------------------------------------------------*/

    siUnstableFlt = siUnstableFlt | isSwLimit(pswRc[i]);

    /* Precompute intermediate variables for needed for the computation */
    /* of direct form filter of order i-1                               */
    /*------------------------------------------------------------------*/

    if (sub(pswRc[i], SW_MIN) == 0)
    {
      siUnstableFlt = 1;
      swRcOverE = 0;
      swDiv = 0;
      swActShift = 2;
    }
    else
    {
      L_temp = LW_MAX;                 /* Load ~1.0 into accum */
      L_temp = L_msu(L_temp, pswRc[i], pswRc[i]);       /* 1.-Rc[i]*Rc[i]  */
      swNormShift = norm_l(L_temp);
      L_temp = L_shl(L_temp, swNormShift);
      swNormProd = extract_h(L_temp);
      swActShift = add(2, swNormShift);
      swDiv = divide_s(0x2000, swNormProd);
      swRcOverE = mult_r(pswRc[i], swDiv);
    }
    /* Check stability   */
    /*---------------------*/
    siUnstableFlt = siUnstableFlt | isSwLimit(swRcOverE);

    /* Compute direct form filter coefficients corresponding to */
    /* a direct form filter of order i-1                        */
    /*----------------------------------------------------------*/

    for (j = 0; j <= i - 1; j++)
    {
      L_temp = L_mult(pswA[j], swDiv);
      L_temp = L_msu(L_temp, pswA[i - j - 1], swRcOverE);
      L_temp = L_shl(L_temp, swActShift);
      pswTmp[j] = round(L_temp);
      siUnstableFlt = siUnstableFlt | isSwLimit(pswTmp[j]);
    }

    /* Swap swA and swTmp buffers */
    /*----------------------------*/

    pswSwap = pswA;
    pswA = pswTmp;
    pswTmp = pswSwap;
  }

  /* Compute reflection coefficient Rc[0] */
  /*--------------------------------------*/

  pswRc[0] = shl(pswA[0], swAshift);   /* write Rc[0] to output array */

  /* Check the stability of 0-th reflection coefficient */
  /*----------------------------------------------------*/

  siUnstableFlt = siUnstableFlt | isSwLimit(pswRc[0]);

  return (siUnstableFlt);
}

/***************************************************************************
 *
 *   FUNCTION NAME: a_sst
 *
 *   PURPOSE:
 *
 *     The purpose of this function is to perform spectral smoothing of the
 *     direct form filter coefficients
 *
 *   INPUTS:
 *
 *     swAshift
 *                     number of shift for coefficients
 *
 *     swAscale
 *                     scaling factor for coefficients
 *
 *     pswDirectFormCoefIn[0:NP-1]
 *
 *                     array of input direct form coefficients
 *
 *   OUTPUTS:
 *
 *     pswDirectFormCoefOut[0:NP-1]
 *
 *                     array of output direct form coefficients
 *
 *   RETURN VALUE:
 *
 *     none
 *
 *   DESCRIPTION:
 *
 *     In a_sst() direct form coefficients are converted to
 *     autocorrelations, and smoothed in that domain.  The function is
 *     used in the spectral postfilter.  A description can be found in
 *     section 3.2.4 as well as in the reference by Y. Tohkura et al.
 *     "Spectral Smoothing Technique in PARCOR Speech
 *     Analysis-Synthesis", IEEE Trans. ASSP, vol. ASSP-26, pp. 591-596,
 *     Dec. 1978.
 *
 *     After smoothing is performed conversion back to direct form
 *     coefficients is done by calling aFlatRc(), followed by rcToADp().
 *
 *     The spectral smoothing filter coefficients with bandwidth set to 300
 *     and a sampling rate of 8000 be :
 *     static ShortwordRom psrSST[NP+1] = { 0x7FFF,
 *         0x7F5C, 0x7D76, 0x7A5B, 0x7622, 0x70EC,
 *         0x6ADD, 0x641F, 0x5CDD, 0x5546, 0x4D86
 *     }
 *
 *   REFERENCES: Sub_Clause 4.2.4 of GSM Recomendation 06.20
 *
 *   KEYWORDS: spectral smoothing, direct form coef, sst, atorc, atocor
 *   KEYWORDS: levinson
 *
 *************************************************************************/

static void a_sst(Shortword swAshift, Shortword swAscale,
                         Shortword pswDirectFormCoefIn[],
                         Shortword pswDirectFormCoefOut[])
{

/*_________________________________________________________________________
 |                                                                         |
 |                           Local Static Variables                        |
 |_________________________________________________________________________|
*/

  static ShortwordRom psrSST[NP + 1] = {0x7FFF,
    0x7F5C, 0x7D76, 0x7A5B, 0x7622, 0x70EC,
    0x6ADD, 0x641F, 0x5CDD, 0x5546, 0x4D86,
  };

/*_________________________________________________________________________
 |                                                                         |
 |                            Automatic Variables                          |
 |_________________________________________________________________________|
*/

  Longword pL_CorrTemp[NP + 1];

  Shortword pswRCNum[NP],
         pswRCDenom[NP];

  short int siLoopCnt;

/*_________________________________________________________________________
 |                                                                         |
 |                              Executable Code                            |
 |_________________________________________________________________________|
*/

  /* convert direct form coefs to reflection coefs */
  /* --------------------------------------------- */

  aToRc(swAshift, pswDirectFormCoefIn, pswRCDenom);

  /* convert to autocorrelation coefficients */
  /* --------------------------------------- */

  rcToCorrDpL(swAshift, swAscale, pswRCDenom, pL_CorrTemp);

  /* do spectral smoothing technique */
  /* ------------------------------- */

  for (siLoopCnt = 1; siLoopCnt <= NP; siLoopCnt++)
  {
    pL_CorrTemp[siLoopCnt] = L_mpy_ls(pL_CorrTemp[siLoopCnt],
                                      psrSST[siLoopCnt]);
  }

  /* Compute the reflection coefficients via AFLAT */
  /*-----------------------------------------------*/

  aFlatRcDp(pL_CorrTemp, pswRCNum);


  /* Convert reflection coefficients to direct form filter coefficients */
  /*-------------------------------------------------------------------*/

  rcToADp(swAscale, pswRCNum, pswDirectFormCoefOut);
}

/**************************************************************************
 *
 *   FUNCTION NAME: agcGain
 *
 *   PURPOSE:
 *
 *     Figure out what the agc gain should be to make the energy in the
 *     output signal match that of the input signal.  Used in the post
 *     filters.
 *
 *   INPUT:
 *
 *      pswStateCurr[0:39]
 *                     Input signal into agc block whose energy is
 *                     to be modified using the gain returned. Signal is not
 *                     modified in this routine.
 *
 *      snsInSigEnergy
 *                     Normalized number with shift count - the energy in
 *                     the input signal.
 *
 *      swEngyRShft
 *                     Number of right shifts to apply to the vectors energy
 *                     to ensure that it remains less than 1.0
 *                     (swEngyRShft is always positive or zero)
 *
 *   OUTPUT:
 *
 *      none
 *
 *   RETURN:
 *
 *      the agc's gain/2 note DIVIDED by 2
 *
 *
 *   REFERENCES: Sub_Clause 4.2.2 and 4.2.4 of GSM Recomendation 06.20
 *
 *   KEYWORDS: postfilter, agc, automaticgaincontrol, leveladjust
 *
 *************************************************************************/

static Shortword agcGain(Shortword pswStateCurr[],
                        struct NormSw snsInSigEnergy, Shortword swEngyRShft)
{

/*_________________________________________________________________________
 |                                                                         |
 |                            Automatic Variables                          |
 |_________________________________________________________________________|
*/

  Longword L_OutEnergy,
         L_AgcGain;

  struct NormSw snsOutEnergy,
         snsAgc;

  Shortword swAgcOut,
         swAgcShftCnt;

/*_________________________________________________________________________
 |                                                                         |
 |                              Executable Code                            |
 |_________________________________________________________________________|
*/

  /* Calculate the energy in the output vector divided by 2 */
  /*--------------------------------------------------------*/

  snsOutEnergy.sh = g_corr1s(pswStateCurr, swEngyRShft, &L_OutEnergy);

  /* reduce energy by a factor of 2 */
  snsOutEnergy.sh = add(snsOutEnergy.sh, 1);

  /* if waveform has nonzero energy, find AGC gain */
  /*-----------------------------------------------*/

  if (L_OutEnergy == 0)
  {
    swAgcOut = 0;
  }
  else
  {

    snsOutEnergy.man = round(L_OutEnergy);

    /* divide input energy by 2 */
    snsInSigEnergy.man = shr(snsInSigEnergy.man, 1);


    /* Calculate AGC gain squared */
    /*----------------------------*/

    snsAgc.man = divide_s(snsInSigEnergy.man, snsOutEnergy.man);
    swAgcShftCnt = norm_s(snsAgc.man);
    snsAgc.man = shl(snsAgc.man, swAgcShftCnt);

    /* find shift count for G^2 */
    /*--------------------------*/

    snsAgc.sh = add(sub(snsInSigEnergy.sh, snsOutEnergy.sh),
                    swAgcShftCnt);
    L_AgcGain = L_deposit_h(snsAgc.man);


    /* Calculate AGC gain */
    /*--------------------*/

    snsAgc.man = sqroot(L_AgcGain);


    /* check if 1/2 sqrt(G^2) >= 1.0                      */
    /* This is equivalent to checking if shiftCnt/2+1 < 0 */
    /*----------------------------------------------------*/

    if (add(snsAgc.sh, 2) < 0)
    {
      swAgcOut = SW_MAX;
    }
    else
    {

      if (0x1 & snsAgc.sh)
      {
        snsAgc.man = mult(snsAgc.man, SQRT_ONEHALF);
      }

      snsAgc.sh = shr(snsAgc.sh, 1);   /* shiftCnt/2 */
      snsAgc.sh = add(snsAgc.sh, 1);   /* shiftCnt/2 + 1 */

      if (snsAgc.sh > 0)
      {
        snsAgc.man = shr(snsAgc.man, snsAgc.sh);
      }
      swAgcOut = snsAgc.man;
    }
  }

  return (swAgcOut);
}

/***************************************************************************
 *
 *   FUNCTION NAME: b_con
 *
 *   PURPOSE:
 *     Expands codeword into an one dimensional array. The 0/1 input is
 *     changed to an element with magnitude +/- 0.5.
 *
 *     input  output
 *
 *       0    -0.5
 *       1    +0.5
 *
 *   INPUT:
 *
 *      swCodeWord
 *                     Input codeword, information in the LSB's
 *
 *      siNumBits
 *                     number of bits in the input codeword and number
 *                     of elements in output vector
 *
 *      pswVectOut[0:siNumBits]
 *
 *                     pointer to bit array
 *
 *   OUTPUT:
 *
 *      pswVectOut[0:siNumBits]
 *
 *                     signed bit array
 *
 *   RETURN:
 *
 *      none
 *
 *   REFERENCES: Sub_Clause 4.1.10 and 4.2.1 of GSM Recomendation 06.20
 *
 *   KEYWORDS: b_con, codeword, expansion
 *
 *************************************************************************/

void   b_con(Shortword swCodeWord, short siNumBits,
                    Shortword pswVectOut[])
{

/*_________________________________________________________________________
 |                                                                         |
 |                            Automatic Variables                          |
 |_________________________________________________________________________|
*/

  short int siLoopCnt;

/*_________________________________________________________________________
 |                                                                         |
 |                              Executable Code                            |
 |_________________________________________________________________________|
*/

  for (siLoopCnt = 0; siLoopCnt < siNumBits; siLoopCnt++)
  {

    if (swCodeWord & 1)                /* temp accumulator get 0.5 */
      pswVectOut[siLoopCnt] = (Shortword) 0x4000;
    else                               /* temp accumulator gets -0.5 */
      pswVectOut[siLoopCnt] = (Shortword) 0xc000;

    swCodeWord = shr(swCodeWord, 1);
  }
}

/***************************************************************************
 *
 *   FUNCTION NAME: fp_ex
 *
 *   PURPOSE:
 *
 *     Looks up a vector in the adaptive excitation codebook (long-term
 *     predictor).
 *
 *   INPUTS:
 *
 *     swOrigLagIn
 *
 *                     Extended resolution lag (lag * oversampling factor)
 *
 *     pswLTPState[-147:39]
 *
 *                     Adaptive codebook (with space at end for looked up
 *                     vector).  Indicies [-147:-1] are the history, [0:39]
 *                     are for the looked up vector.
 *
 *     psrPitchIntrpFirBase[0:59]
 *     ppsrPVecIntFilt[0:9][0:5] ([tap][phase])
 *
 *                     Interpolating FIR filter coefficients.
 *
 *   OUTPUTS:
 *
 *     pswLTPState[0:39]
 *
 *                     Array containing the contructed output vector
 *
 *   RETURN VALUE:
 *     none
 *
 *   DESCRIPTION:
 *
 *     The adaptive codebook consists of the history of the excitation.
 *     The vector is looked up by going back into this history
 *     by the amount of the input lag.  If the input lag is fractional,
 *     then the samples to be looked up are interpolated from the existing
 *     samples in the history.
 *
 *     If the lag is less than the length of the vector to be generated
 *     (i.e. less than the subframe length), then the lag is doubled
 *     after the first n samples have been looked up (n = input lag).
 *     In this way, the samples being generated are not part of the
 *     codebook.  This is described in section 4.1.8.
 *
 *   REFERENCES: Sub_Clause 4.1.8.5 and 4.2.1  of GSM Recomendation 06.20
 *
 *   Keywords: pitch, excitation vector, long term filter, history,
 *   Keywords: fractional lag, get_ipjj
 *
 *************************************************************************/



void   fp_ex(Shortword swOrigLagIn,
                    Shortword pswLTPState[])
{

/*_________________________________________________________________________
 |                                                                         |
 |                            Automatic Variables                          |
 |_________________________________________________________________________|
*/

  Longword L_Temp;
  Shortword swIntLag,
         swRemain,
         swRunningLag;
  short int siSampsSoFar,
         siSampsThisPass,
         i,
         j;

/*_________________________________________________________________________
 |                                                                         |
 |                              Executable Code                            |
 |_________________________________________________________________________|
*/

  /* Loop: execute until all samples in the vector have been looked up */
  /*-------------------------------------------------------------------*/

  swRunningLag = swOrigLagIn;
  siSampsSoFar = 0;
  while (siSampsSoFar < S_LEN)
  {

    /* Get integer lag and remainder.  These are used in addressing */
    /* the LTP state and the interpolating filter, respectively     */
    /*--------------------------------------------------------------*/

    get_ipjj(swRunningLag, &swIntLag, &swRemain);


    /* Get the number of samples to look up in this pass */
    /*---------------------------------------------------*/

    if (sub(swIntLag, S_LEN) < 0)
      siSampsThisPass = swIntLag - siSampsSoFar;
    else
      siSampsThisPass = S_LEN - siSampsSoFar;

    /* Look up samples by interpolating (fractional lag), or copying */
    /* (integer lag).                                                */
    /*---------------------------------------------------------------*/

    if (swRemain == 0)
    {

      /* Integer lag: copy samples from history */
      /*----------------------------------------*/

      for (i = siSampsSoFar; i < siSampsSoFar + siSampsThisPass; i++)
        pswLTPState[i] = pswLTPState[i - swIntLag];
    }
    else
    {

      /* Fractional lag: interpolate to get samples */
      /*--------------------------------------------*/

      for (i = siSampsSoFar; i < siSampsSoFar + siSampsThisPass; i++)
      {

        /* first tap with rounding offset */
        /*--------------------------------*/
        L_Temp = L_mac((long) 32768,
                       pswLTPState[i - swIntLag - P_INT_MACS / 2],
                       ppsrPVecIntFilt[0][swRemain]);

        for (j = 1; j < P_INT_MACS - 1; j++)
        {

          L_Temp = L_mac(L_Temp,
                         pswLTPState[i - swIntLag - P_INT_MACS / 2 + j],
                         ppsrPVecIntFilt[j][swRemain]);

        }

        pswLTPState[i] = extract_h(L_mac(L_Temp,
                             pswLTPState[i - swIntLag + P_INT_MACS / 2 - 1],
                                ppsrPVecIntFilt[P_INT_MACS - 1][swRemain]));
      }
    }

    /* Done with this pass: update loop controls */
    /*-------------------------------------------*/

    siSampsSoFar += siSampsThisPass;
    swRunningLag = add(swRunningLag, swOrigLagIn);
  }
}

/***************************************************************************
 *
 *    FUNCTION NAME: g_corr1 (no scaling)
 *
 *    PURPOSE:
 *
 *     Calculates energy in subframe vector.  Differs from g_corr1s,
 *     in that there the estimate of the maximum possible
 *     energy is < 1.0
 *
 *
 *    INPUT:
 *
 *       pswIn[0:39]
 *                     A subframe vector.
 *
 *
 *    OUTPUT:
 *
 *       *pL_out
 *                     A Longword containing the normalized energy
 *                     in the input vector.
 *
 *    RETURN:
 *
 *       swOut
 *                     Number of right shifts which the accumulator was
 *                     shifted to normalize it.  Negative number implies
 *                     a left shift, and therefore an energy larger than
 *                     1.0.
 *
 *    REFERENCES: Sub_Clause 4.1.10.2 and 4.2.1 of GSM Recomendation 06.20
 *
 *    KEYWORDS: energy, autocorrelation, correlation, g_corr1
 *
 *
 *************************************************************************/

Shortword g_corr1(Shortword *pswIn, Longword *pL_out)
{


/*_________________________________________________________________________
 |                                                                         |
 |                            Automatic Variables                          |
 |_________________________________________________________________________|
*/

  Longword L_sum;
  Shortword swEngyLShft;
  int    i;


/*_________________________________________________________________________
 |                                                                         |
 |                              Executable Code                            |
 |_________________________________________________________________________|
*/


  /* Calculate energy in subframe vector (40 samples) */
  /*--------------------------------------------------*/

  L_sum = L_mult(pswIn[0], pswIn[0]);
  for (i = 1; i < S_LEN; i++)
  {
    L_sum = L_mac(L_sum, pswIn[i], pswIn[i]);
  }



  if (L_sum != 0)
  {

    /* Normalize the energy in the output Longword */
    /*---------------------------------------------*/

    swEngyLShft = norm_l(L_sum);
    *pL_out = L_shl(L_sum, swEngyLShft);        /* normalize output
                                                 * Longword */
  }
  else
  {

    /* Special case: energy is zero */
    /*------------------------------*/

    *pL_out = L_sum;
    swEngyLShft = 0;
  }

  return (swEngyLShft);
}

/***************************************************************************
 *
 *    FUNCTION NAME: g_corr1s (g_corr1 with scaling)
 *
 *    PURPOSE:
 *
 *     Calculates energy in subframe vector.  Differs from g_corr1,
 *     in that there is an estimate of the maximum possible energy in the
 *     vector.
 *
 *    INPUT:
 *
 *       pswIn[0:39]
 *                     A subframe vector.
 *
 *       swEngyRShft
 *
 *                     Number of right shifts to apply to the vectors energy
 *                     to ensure that it remains less than 1.0
 *                     (swEngyRShft is always positive or zero)
 *
 *    OUTPUT:
 *
 *       *pL_out
 *                     A Longword containing the normalized energy
 *                     in the input vector.
 *
 *    RETURN:
 *
 *       swOut
 *                     Number of right shifts which the accumulator was
 *                     shifted to normalize it.  Negative number implies
 *                     a left shift, and therefore an energy larger than
 *                     1.0.
 *
 *    REFERENCES: Sub-Clause 4.1.8, 4.2.1, 4.2.2, and 4.2.4
 *                of GSM Recomendation 06.20
 *
 *    keywords: energy, autocorrelation, correlation, g_corr1
 *
 *
 *************************************************************************/

Shortword g_corr1s(Shortword pswIn[], Shortword swEngyRShft,
                          Longword *pL_out)
{


/*_________________________________________________________________________
 |                                                                         |
 |                            Automatic Variables                          |
 |_________________________________________________________________________|
*/

  Longword L_sum;
  Shortword swTemp,
         swEngyLShft;
  Shortword swInputRShft;

  int    i;


/*_________________________________________________________________________
 |                                                                         |
 |                              Executable Code                            |
 |_________________________________________________________________________|
*/


  /* Calculate energy in subframe vector (40 samples) */
  /*--------------------------------------------------*/

  if (sub(swEngyRShft, 1) <= 0)
  {

    /* use the energy shift factor, although it is an odd shift count */
    /*----------------------------------------------------------------*/

    swTemp = shr(pswIn[0], swEngyRShft);
    L_sum = L_mult(pswIn[0], swTemp);
    for (i = 1; i < S_LEN; i++)
    {
      swTemp = shr(pswIn[i], swEngyRShft);
      L_sum = L_mac(L_sum, pswIn[i], swTemp);
    }

  }
  else
  {

    /* convert energy shift factor to an input shift factor */
    /*------------------------------------------------------*/

    swInputRShft = shift_r(swEngyRShft, -1);
    swEngyRShft = shl(swInputRShft, 1);

    swTemp = shr(pswIn[0], swInputRShft);
    L_sum = L_mult(swTemp, swTemp);
    for (i = 1; i < S_LEN; i++)
    {
      swTemp = shr(pswIn[i], swInputRShft);
      L_sum = L_mac(L_sum, swTemp, swTemp);
    }
  }

  if (L_sum != 0)
  {

    /* Normalize the energy in the output Longword */
    /*---------------------------------------------*/

    swTemp = norm_l(L_sum);
    *pL_out = L_shl(L_sum, swTemp);    /* normalize output Longword */
    swEngyLShft = sub(swTemp, swEngyRShft);
  }
  else
  {

    /* Special case: energy is zero */
    /*------------------------------*/

    *pL_out = L_sum;
    swEngyLShft = 0;
  }

  return (swEngyLShft);
}

/***************************************************************************
 *
 *   FUNCTION NAME: getSfrmLpc
 *
 *   PURPOSE:
 *
 *     Given frame information from past and present frame, interpolate
 *     (or copy) the frame based LPC coefficients into subframe
 *     lpc coeffs, i.e. the ones which will be used by the subframe
 *     as opposed to those coded and transmitted.
 *
 *   INPUTS:
 *
 *     siSoftInterpolation
 *
 *                     interpolate 1/0, a coded parameter.
 *
 *     swPrevR0,swNewR0
 *
 *                     Rq0 for the last frame and for this frame.
 *                     These are the decoded values, not the codewords.
 *
 *     Previous lpc coefficients from the previous frame:
 *       in all filters below array[0] is the t=-1 element array[9]
 *       t=-10 element.
 *
 *     pswPrevFrmKs[0:9]
 *
 *                     decoded version of the rc's tx'd last frame
 *
 *     pswPrevFrmAs[0:9]
 *
 *                     the above K's converted to A's.  i.e. direct
 *                     form coefficients.
 *
 *     pswPrevFrmPFNum[0:9], pswPrevFrmPFDenom[0:9]
 *
 *                     numerator and denominator coefficients used in the
 *                     postfilter
 *
 *     Current lpc coefficients from the current frame:
 *
 *     pswNewFrmKs[0:9], pswNewFrmAs[0:9],
 *     pswNewFrmPFNum[0:9], pswNewFrmPFDenom[0:9] same as above.
 *
 *   OUTPUTS:
 *
 *     psnsSqrtRs[0:3]
 *
 *                      a normalized number (struct NormSw)
 *                      containing an estimate of RS for each subframe.
 *                      (number and a shift)
 *
 *     ppswSynthAs[0:3][0:9]
 *
 *                      filter coefficients used by the synthesis filter.
 *
 *     ppswPFNumAs[0:3][0:9]
 *
 *                      filter coefficients used by the postfilters
 *                      numerator.
 *
 *     ppswPFDenomAs[0:3][0:9]
 *
 *                      filter coefficients used by postfilters denominator.
 *
 *   RETURN VALUE:
 *
 *     None
 *
 *   DESCRIPTION:
 *
 *     For interpolated subframes, the direct form coefficients
 *     are converted to reflection coeffiecients to check for
 *     filter stability. If unstable, the uninterpolated coef.
 *     are used for that subframe.
 *
 *     Interpolation is described in section 4.1.6, "Soft Interpolation
 *     of the Spectral Parameters"
 *
 *    REFERENCES: Sub_clause 4.2.1 of GSM Recomendation 06.20
 *
 *   KEYWORDS: soft interpolation, int_lpc, interpolate, atorc,res_eng,i_mov
 *
 *************************************************************************/

void   getSfrmLpc(short int siSoftInterpolation,
                         Shortword swPrevR0, Shortword swNewR0,
          /* last frm */ Shortword pswPrevFrmKs[], Shortword pswPrevFrmAs[],
                         Shortword pswPrevFrmPFNum[],
                         Shortword pswPrevFrmPFDenom[],

            /* this frm */ Shortword pswNewFrmKs[], Shortword pswNewFrmAs[],
                         Shortword pswNewFrmPFNum[],
                         Shortword pswNewFrmPFDenom[],

                   /* output */ struct NormSw *psnsSqrtRs,
                         Shortword *ppswSynthAs[], Shortword *ppswPFNumAs[],
                         Shortword *ppswPFDenomAs[])
{

/*_________________________________________________________________________
 |                                                                         |
 |                           Local Static Variables                        |
 |_________________________________________________________________________|
*/


/*_________________________________________________________________________
 |                                                                         |
 |                            Automatic Variables                          |
 |_________________________________________________________________________|
*/

  short int siSfrm,
         siStable,
         i;

  Longword L_Temp1,
         L_Temp2;

/*_________________________________________________________________________
 |                                                                         |
 |                              Executable Code                            |
 |_________________________________________________________________________|
*/

  if (siSoftInterpolation)
  {
    /* yes, interpolating */
    /* ------------------ */

    siSfrm = 0;

    siStable = interpolateCheck(pswPrevFrmKs, pswPrevFrmAs,
                                pswPrevFrmAs, pswNewFrmAs,
                                psrOldCont[siSfrm], psrNewCont[siSfrm],
                                swPrevR0,
                                &psnsSqrtRs[siSfrm],
                                ppswSynthAs[siSfrm]);
    if (siStable)
    {

      /* interpolate between direct form coefficient sets */
      /* for both numerator and denominator coefficients  */
      /* assume output will be stable                     */
      /* ------------------------------------------------ */

      for (i = 0; i < NP; i++)
      {
        L_Temp1 = L_mult(pswNewFrmPFNum[i], psrNewCont[siSfrm]);
        ppswPFNumAs[siSfrm][i] = mac_r(L_Temp1, pswPrevFrmPFNum[i],
                                       psrOldCont[siSfrm]);
        L_Temp2 = L_mult(pswNewFrmPFDenom[i], psrNewCont[siSfrm]);
        ppswPFDenomAs[siSfrm][i] = mac_r(L_Temp2, pswPrevFrmPFDenom[i],
                                         psrOldCont[siSfrm]);
      }
    }
    else
    {
      /* this subframe is unstable */
      /* ------------------------- */
      for (i = 0; i < NP; i++)
      {
        ppswPFNumAs[siSfrm][i] = pswPrevFrmPFNum[i];
        ppswPFDenomAs[siSfrm][i] = pswPrevFrmPFDenom[i];
      }
    }
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

        /* interpolate between direct form coefficient sets */
        /* for both numerator and denominator coefficients  */
        /* assume output will be stable                     */
        /* ------------------------------------------------ */

        for (i = 0; i < NP; i++)
        {
          L_Temp1 = L_mult(pswNewFrmPFNum[i], psrNewCont[siSfrm]);
          ppswPFNumAs[siSfrm][i] = mac_r(L_Temp1, pswPrevFrmPFNum[i],
                                         psrOldCont[siSfrm]);
          L_Temp2 = L_mult(pswNewFrmPFDenom[i], psrNewCont[siSfrm]);
          ppswPFDenomAs[siSfrm][i] = mac_r(L_Temp2, pswPrevFrmPFDenom[i],
                                           psrOldCont[siSfrm]);
        }
      }
      else
      {
        /* this subframe has unstable filter coeffs, would like to
         * interpolate but can not  */
        /* -------------------------------------- */
        for (i = 0; i < NP; i++)
        {
          ppswPFNumAs[siSfrm][i] = pswNewFrmPFNum[i];
          ppswPFDenomAs[siSfrm][i] = pswNewFrmPFDenom[i];
        }
      }
    }
    /* the last subframe never interpolate */
    /* ----------------------------------- */
    siSfrm = 3;
    for (i = 0; i < NP; i++)
    {
      ppswPFNumAs[siSfrm][i] = pswNewFrmPFNum[i];
      ppswPFDenomAs[siSfrm][i] = pswNewFrmPFDenom[i];
      ppswSynthAs[siSfrm][i] = pswNewFrmAs[i];
    }

    res_eng(pswNewFrmKs, swNewR0, &psnsSqrtRs[siSfrm]);

  }
  /* SoftInterpolation == 0  - no interpolation */
  /* ------------------------------------------ */
  else
  {
    siSfrm = 0;
    for (i = 0; i < NP; i++)
    {
      ppswPFNumAs[siSfrm][i] = pswPrevFrmPFNum[i];
      ppswPFDenomAs[siSfrm][i] = pswPrevFrmPFDenom[i];
      ppswSynthAs[siSfrm][i] = pswPrevFrmAs[i];
    }

    res_eng(pswPrevFrmKs, swPrevR0, &psnsSqrtRs[siSfrm]);

    /* for subframe 1 and all subsequent sfrms, use result from new frm */
    /* ---------------------------------------------------------------- */


    res_eng(pswNewFrmKs, swNewR0, &psnsSqrtRs[1]);

    for (siSfrm = 1; siSfrm < N_SUB; siSfrm++)
    {


      psnsSqrtRs[siSfrm].man = psnsSqrtRs[1].man;
      psnsSqrtRs[siSfrm].sh = psnsSqrtRs[1].sh;

      for (i = 0; i < NP; i++)
      {
        ppswPFNumAs[siSfrm][i] = pswNewFrmPFNum[i];
        ppswPFDenomAs[siSfrm][i] = pswNewFrmPFDenom[i];
        ppswSynthAs[siSfrm][i] = pswNewFrmAs[i];
      }
    }
  }
}

/***************************************************************************
 *
 *   FUNCTION NAME: get_ipjj
 *
 *   PURPOSE:
 *
 *     This subroutine calculates IP, the single-resolution lag rounded
 *     down to the nearest integer, and JJ, the remainder when the
 *     extended resolution lag is divided by the oversampling factor
 *
 *   INPUTS:
 *
 *     swLagIn
 *                     extended resolution lag as an integer, i.e.
 *                     fractional lag x oversampling factor
 *
 *   OUTPUTS:
 *
 *     *pswIp
 *                     fractional lag rounded down to nearest integer, IP
 *
 *     *pswJj
 *                     the remainder JJ
 *
 *   RETURN VALUE:
 *
 *     none
 *
 *   DESCRIPTION:
 *
 *     ip = integer[lag/OS_FCTR]
 *     jj = integer_round[((lag/OS_FCTR)-ip)*(OS_FCTR)]
 *          if the rounding caused an 'overflow'
 *            set remainder jj to 0 and add 'carry' to ip
 *
 *     This routine is involved in the mechanics of fractional and
 *     integer LTP searchs.  The LTP is described in section 5.
 *
 *   REFERENCES: Sub-clause 4.1.8 and 4.2.2 of GSM Recomendation 06.20
 *
 *   KEYWORDS: lag, fractional, remainder, ip, jj, get_ipjj
 *
 *************************************************************************/

void   get_ipjj(Shortword swLagIn,
                       Shortword *pswIp, Shortword *pswJj)
{

/*_________________________________________________________________________
 |                                                                         |
 |                              Local Constants                            |
 |_________________________________________________________________________|
*/

#define  OS_FCTR_INV  (Shortword)0x1555/* SW_MAX/OS_FCTR */

/*_________________________________________________________________________
 |                                                                         |
 |                            Automatic Variables                          |
 |_________________________________________________________________________|
*/

  Longword L_Temp;

  Shortword swTemp,
         swTempIp,
         swTempJj;

/*_________________________________________________________________________
 |                                                                         |
 |                              Executable Code                            |
 |_________________________________________________________________________|
*/

  /* calculate ip */
  /* ------------ */

  L_Temp = L_mult(OS_FCTR_INV, swLagIn);        /* lag/OS_FCTR */
  swTempIp = extract_h(L_Temp);

  /* calculate jj */
  /* ------------ */

  swTemp = extract_l(L_Temp);          /* loose ip */
  swTemp = shr(swTemp, 1);             /* isolate jj fraction */
  swTemp = swTemp & SW_MAX;
  L_Temp = L_mult(swTemp, OS_FCTR);    /* ((lag/OS_FCTR)-ip))*(OS_FCTR) */
  swTemp = round(L_Temp);              /* round and pick-off jj */
  if (sub(swTemp, OS_FCTR) == 0)
  {                                    /* if 'overflow ' */
    swTempJj = 0;                      /* set remainder,jj to 0 */
    swTempIp = add(swTempIp, 1);       /* 'carry' overflow into ip */
  }
  else
  {
    swTempJj = swTemp;                 /* read-off remainder,jj */
  }

  /* return ip and jj */
  /* ---------------- */

  *pswIp = swTempIp;
  *pswJj = swTempJj;
}

/***************************************************************************
 *
 *   FUNCTION NAME: interpolateCheck
 *
 *   PURPOSE:
 *
 *     Interpolates between direct form coefficient sets.
 *     Before releasing the interpolated coefficients, they are checked.
 *     If unstable, the "old" parameters are used.
 *
 *   INPUTS:
 *
 *     pswRefKs[0:9]
 *                     decoded version of the rc's tx'd last frame
 *
 *     pswRefCoefsA[0:9]
 *                     above K's converted to direct form coefficients
 *
 *     pswOldCoefsA[0:9]
 *                     array of old Coefseters
 *
 *     pswNewCoefsA[0:9]
 *                     array of new Coefseters
 *
 *     swOldPer
 *                     amount old coefs supply to the output
 *
 *     swNewPer
 *                     amount new coefs supply to the output
 *
 *     ASHIFT
 *                     shift for reflection coef. conversion
 *
 *     swRq
 *                     quantized energy to use for subframe
 * *
 *   OUTPUTS:
 *
 *     psnsSqrtRsOut
 *                     output pointer to sqrt(RS) normalized
 *
 *     pswCoefOutA[0:9]
 *                     output coefficients
 *
 *   RETURN VALUE:
 *
 *     siInterp_flg
 *                     temporary subframe interpolation flag
 *                     0 - coef. interpolated, 1 -coef. not interpolated
 *
 *   DESCRIPTION:
 *
 *     For interpolated subframes, the direct form coefficients
 *     are converted to reflection coefficients to check for
 *     filter stability. If unstable, the uninterpolated coef.
 *     are used for that subframe.  Section 4.1.6 describes
 *     interpolation.
 *
 *   REFERENCES: Sub-clause 4.1.6 and 4.2.3 of GSM Recomendation 06.20
 *
 *   KEYWORDS: soft interpolation, int_lpc, interpolate, atorc,res_eng,i_mov
 *
 *************************************************************************/

short int interpolateCheck(Shortword pswRefKs[],
                                  Shortword pswRefCoefsA[],
                         Shortword pswOldCoefsA[], Shortword pswNewCoefsA[],
                                  Shortword swOldPer, Shortword swNewPer,
                                  Shortword swRq,
                                  struct NormSw *psnsSqrtRsOut,
                                  Shortword pswCoefOutA[])
{

/*_________________________________________________________________________
 |                                                                         |
 |                            Automatic Variables                          |
 |_________________________________________________________________________|
*/

  Shortword pswRcTemp[NP];

  Longword L_Temp;

  short int siInterp_flg,
         i;

/*_________________________________________________________________________
 |                                                                         |
 |                              Executable Code                            |
 |_________________________________________________________________________|
*/

  /* Interpolation loop, NP is order of LPC filter */
  /* --------------------------------------------- */

  for (i = 0; i < NP; i++)
  {
    L_Temp = L_mult(pswNewCoefsA[i], swNewPer);
    pswCoefOutA[i] = mac_r(L_Temp, pswOldCoefsA[i], swOldPer);
  }

  /* Convert to reflection coefficients and check stability */
  /* ------------------------------------------------------ */

  if (aToRc(ASHIFT, pswCoefOutA, pswRcTemp) != 0)
  {

    /* Unstable, use uninterpolated parameters and compute RS update the
     * state with the frame data closest to this subfrm */
    /* --------------------------------------------------------- */

    res_eng(pswRefKs, swRq, psnsSqrtRsOut);

    for (i = 0; i < NP; i++)
    {
      pswCoefOutA[i] = pswRefCoefsA[i];
    }
    siInterp_flg = 0;
  }
  else
  {

    /* Stable, compute RS */
    /* ------------------ */
    res_eng(pswRcTemp, swRq, psnsSqrtRsOut);

    /* Set temporary subframe interpolation flag */
    /* ----------------------------------------- */
    siInterp_flg = 1;
  }

  /* Return subframe interpolation flag */
  /* ---------------------------------- */
  return (siInterp_flg);
}

/***************************************************************************
 *
 *   FUNCTION NAME: lagDecode
 *
 *   PURPOSE:
 *
 *     The purpose of this function is to decode the lag received from the
 *     speech encoder into a full resolution lag for the speech decoder
 *
 *   INPUTS:
 *
 *     swDeltaLag
 *
 *                     lag received from channel decoder
 *
 *     giSfrmCnt
 *
 *                     current sub-frame count
 *
 *     swLastLag
 *
 *                     previous lag to un-delta this sub-frame's lag
 *
 *     psrLagTbl[0:255]
 *
 *                     table used to look up full resolution lag
 *
 *   OUTPUTS:
 *
 *     swLastLag
 *
 *                     new previous lag for next sub-frame
 *
 *   RETURN VALUE:
 *
 *     swLag
 *
 *                     decoded full resolution lag
 *
 *   DESCRIPTION:
 *
 *     If first subframe, use lag as index to look up table directly.
 *
 *     If it is one of the other subframes, the codeword represents a
 *     delta offset.  The previously decoded lag is used as a starting
 *     point for decoding the current lag.
 *
 *   REFERENCES: Sub-clause 4.2.1 of GSM Recomendation 06.20
 *
 *   KEYWORDS: deltalags, lookup lag
 *
 *************************************************************************/

static Shortword lagDecode(Shortword swDeltaLag)
{

/*_________________________________________________________________________
 |                                                                         |
 |                              Local Constants                            |
 |_________________________________________________________________________|
*/

#define  DELTA_LEVELS_D2  DELTA_LEVELS/2
#define  MAX_LAG          0x00ff
#define  MIN_LAG          0x0000

/*_________________________________________________________________________
 |                                                                         |
 |                           Local Static Variables                        |
 |_________________________________________________________________________|
*/

  static Shortword swLastLag;

/*_________________________________________________________________________
 |                                                                         |
 |                            Automatic Variables                          |
 |_________________________________________________________________________|
*/

  Shortword swLag;

/*_________________________________________________________________________
 |                                                                         |
 |                              Executable Code                            |
 |_________________________________________________________________________|
*/

  /* first sub-frame */
  /* --------------- */

  if (giSfrmCnt == 0)
  {
    swLastLag = swDeltaLag;
  }

  /* remaining sub-frames */
  /* -------------------- */

  else
  {

    /* get lag biased around 0 */
    /* ----------------------- */

    swLag = sub(swDeltaLag, DELTA_LEVELS_D2);

    /* get real lag relative to last */
    /* ----------------------------- */

    swLag = add(swLag, swLastLag);

    /* clip to max or min */
    /* ------------------ */

    if (sub(swLag, MAX_LAG) > 0)
    {
      swLastLag = MAX_LAG;
    }
    else if (sub(swLag, MIN_LAG) < 0)
    {
      swLastLag = MIN_LAG;
    }
    else
    {
      swLastLag = swLag;
    }
  }

  /* return lag after look up */
  /* ------------------------ */

  swLag = psrLagTbl[swLastLag];
  return (swLag);
}

/***************************************************************************
 *
 *   FUNCTION NAME: lookupVq
 *
 *   PURPOSE:
 *
 *     The purpose of this function is to recover the reflection coeffs from
 *     the received LPC codewords.
 *
 *   INPUTS:
 *
 *     pswVqCodeWds[0:2]
 *
 *                         the codewords for each of the segments
 *
 *   OUTPUTS:
 *
 *     pswRCOut[0:NP-1]
 *
 *                        the decoded reflection coefficients
 *
 *   RETURN VALUE:
 *
 *     none.
 *
 *   DESCRIPTION:
 *
 *     For each segment do the following:
 *       setup the retrieval pointers to the correct vector
 *       get that vector
 *
 *   REFERENCES: Sub-clause 4.2.3 of GSM Recomendation 06.20
 *
 *   KEYWORDS: vq, vectorquantizer, lpc
 *
 *************************************************************************/

static void lookupVq(Shortword pswVqCodeWds[], Shortword pswRCOut[])
{
/*_________________________________________________________________________
 |                                                                         |
 |                              Local Constants                            |
 |_________________________________________________________________________|
*/

#define  LSP_MASK  0x00ff

/*_________________________________________________________________________
 |                                                                         |
 |                            Automatic Variables                          |
 |_________________________________________________________________________|
*/

  short int siSeg,
         siIndex,
         siVector,
         siVector1,
         siVector2,
         siWordPtr;

  ShortwordRom *psrQTable;

/*_________________________________________________________________________
 |                                                                         |
 |                              Executable Code                            |
 |_________________________________________________________________________|
*/

  /* for each segment */
  /* ---------------- */

  for (siSeg = 0; siSeg < QUANT_NUM_OF_TABLES; siSeg++)
  {

    siVector = pswVqCodeWds[siSeg];
    siIndex = psvqIndex[siSeg].l;

    if (sub(siSeg, 2) == 0)
    {                                  /* segment 3 */

      /* set table */
      /* --------- */

      psrQTable = psrQuant3;

      /* set offset into table */
      /* ---------------------- */

      siWordPtr = add(siVector, siVector);

      /* look up coeffs */
      /* -------------- */

      siVector1 = psrQTable[siWordPtr];
      siVector2 = psrQTable[siWordPtr + 1];

      pswRCOut[siIndex - 1] = psrSQuant[shr(siVector1, 8) & LSP_MASK];
      pswRCOut[siIndex] = psrSQuant[siVector1 & LSP_MASK];
      pswRCOut[siIndex + 1] = psrSQuant[shr(siVector2, 8) & LSP_MASK];
      pswRCOut[siIndex + 2] = psrSQuant[siVector2 & LSP_MASK];
    }
    else
    {                                  /* segments 1 and 2 */

      /* set tables */
      /* ---------- */

      if (siSeg == 0)
      {
        psrQTable = psrQuant1;
      }
      else
      {
        psrQTable = psrQuant2;

      }

      /* set offset into table */
      /* --------------------- */

      siWordPtr = add(siVector, siVector);
      siWordPtr = add(siWordPtr, siVector);
      siWordPtr = shr(siWordPtr, 1);

      /* look up coeffs */
      /* -------------- */

      siVector1 = psrQTable[siWordPtr];
      siVector2 = psrQTable[siWordPtr + 1];

      if ((siVector & 0x0001) == 0)
      {
        pswRCOut[siIndex - 1] = psrSQuant[shr(siVector1, 8) & LSP_MASK];
        pswRCOut[siIndex] = psrSQuant[siVector1 & LSP_MASK];
        pswRCOut[siIndex + 1] = psrSQuant[shr(siVector2, 8) & LSP_MASK];
      }
      else
      {
        pswRCOut[siIndex - 1] = psrSQuant[siVector1 & LSP_MASK];
        pswRCOut[siIndex] = psrSQuant[shr(siVector2, 8) & LSP_MASK];
        pswRCOut[siIndex + 1] = psrSQuant[siVector2 & LSP_MASK];
      }
    }
  }
}

/***************************************************************************
 *
 *   FUNCTION NAME: lpcFir
 *
 *   PURPOSE:
 *
 *     The purpose of this function is to perform direct form fir filtering
 *     assuming a NP order filter and given state, coefficients, and input.
 *
 *   INPUTS:
 *
 *     NP
 *                     order of the lpc filter
 *
 *     S_LEN
 *                     number of samples to filter
 *
 *     pswInput[0:S_LEN-1]
 *
 *                     input array of points to be filtered.
 *                     pswInput[0] is the oldest point (first to be filtered)
 *                     pswInput[siLen-1] is the last point filtered (newest)
 *
 *     pswCoef[0:NP-1]
 *
 *                     array of direct form coefficients
 *                     pswCoef[0] = coeff for delay n = -1
 *                     pswCoef[NP-1] = coeff for delay n = -NP
 *
 *     ASHIFT
 *                     number of shifts input A's have been shifted down by
 *
 *     LPC_ROUND
 *                     rounding constant
 *
 *     pswState[0:NP-1]
 *
 *                     array of the filter state following form of pswCoef
 *                     pswState[0] = state of filter for delay n = -1
 *                     pswState[NP-1] = state of filter for delay n = -NP
 *
 *   OUTPUTS:
 *
 *     pswState[0:NP-1]
 *
 *                     updated filter state, ready to filter
 *                     pswInput[siLen], i.e. the next point
 *
 *     pswFiltOut[0:S_LEN-1]
 *
 *                     the filtered output
 *                     same format as pswInput, pswFiltOut[0] is oldest point
 *
 *   RETURN VALUE:
 *
 *     none
 *
 *   DESCRIPTION:
 *
 *     because of the default sign of the coefficients the
 *     formula for the filter is :
 *     i=0, i < S_LEN
 *        out[i] = rounded(state[i]*coef[0])
 *        j=1, j < NP
 *           out[i] += state[j]*coef[j] (state is taken from either input
 *                                       state[] or input in[] arrays)
 *        rescale(out[i])
 *        out[i] += in[i]
 *     update final state array using in[]
 *
 *   REFERENCES: Sub-clause 4.1.7 and 4.2.4 of GSM Recomendation 06.20
 *
 *   KEYWORDS: lpc, directform, fir, lpcFir, inversefilter, lpcFilt
 *   KEYWORDS: dirForm, dir_mod, dir_clr, dir_neg, dir_set, i_dir_mod
 *
 *************************************************************************/

void   lpcFir(Shortword pswInput[], Shortword pswCoef[],
                     Shortword pswState[], Shortword pswFiltOut[])
{

/*_________________________________________________________________________
 |                                                                         |
 |                            Automatic Variables                          |
 |_________________________________________________________________________|
*/

  Longword L_Sum;
  short int siStage,
         siSmp;

/*_________________________________________________________________________
 |                                                                         |
 |                              Executable Code                            |
 |_________________________________________________________________________|
*/

  /* filter 1st sample */
  /* ----------------- */

  /* sum past state outputs */
  /* ---------------------- */
  /* 0th coef, with rounding */
  L_Sum = L_mac(LPC_ROUND, pswState[0], pswCoef[0]);

  for (siStage = 1; siStage < NP; siStage++)
  {                                    /* remaining coefs */
    L_Sum = L_mac(L_Sum, pswState[siStage], pswCoef[siStage]);
  }

  /* add input to partial output */
  /* --------------------------- */

  L_Sum = L_shl(L_Sum, ASHIFT);
  L_Sum = L_msu(L_Sum, pswInput[0], 0x8000);

  /* save 1st output sample */
  /* ---------------------- */

  pswFiltOut[0] = extract_h(L_Sum);

  /* filter remaining samples */
  /* ------------------------ */

  for (siSmp = 1; siSmp < S_LEN; siSmp++)
  {

    /* sum past outputs */
    /* ---------------- */
    /* 0th coef, with rounding */
    L_Sum = L_mac(LPC_ROUND, pswInput[siSmp - 1], pswCoef[0]);
    /* remaining coefs */
    for (siStage = 1; ((0 < (siSmp - siStage)) && siStage < NP); siStage++)
    {
      L_Sum = L_mac(L_Sum, pswInput[siSmp - siStage - 1], pswCoef[siStage]);
    }

    /* sum past states, if any */
    /* ----------------------- */

    for (siStage = siSmp; siStage < NP; siStage++)
    {
      L_Sum = L_mac(L_Sum, pswState[siStage - siSmp], pswCoef[siStage]);
    }

    /* add input to partial output */
    /* --------------------------- */

    L_Sum = L_shl(L_Sum, ASHIFT);
    L_Sum = L_msu(L_Sum, pswInput[siSmp], 0x8000);

    /* save current output sample */
    /* -------------------------- */

    pswFiltOut[siSmp] = extract_h(L_Sum);
  }

  /* save final state */
  /* ---------------- */

  for (siStage = 0; siStage < NP; siStage++)
  {
    pswState[siStage] = pswInput[S_LEN - siStage - 1];
  }

}

/***************************************************************************
 *
 *   FUNCTION NAME: lpcIir
 *
 *   PURPOSE:
 *
 *     The purpose of this function is to perform direct form IIR filtering
 *     assuming a NP order filter and given state, coefficients, and input
 *
 *   INPUTS:
 *
 *     NP
 *                     order of the lpc filter
 *
 *     S_LEN
 *                     number of samples to filter
 *
 *     pswInput[0:S_LEN-1]
 *
 *                     input array of points to be filtered
 *                     pswInput[0] is the oldest point (first to be filtered)
 *                     pswInput[siLen-1] is the last point filtered (newest)
 *
 *     pswCoef[0:NP-1]
 *                     array of direct form coefficients
 *                     pswCoef[0] = coeff for delay n = -1
 *                     pswCoef[NP-1] = coeff for delay n = -NP
 *
 *     ASHIFT
 *                     number of shifts input A's have been shifted down by
 *
 *     LPC_ROUND
 *                     rounding constant
 *
 *     pswState[0:NP-1]
 *
 *                     array of the filter state following form of pswCoef
 *                     pswState[0] = state of filter for delay n = -1
 *                     pswState[NP-1] = state of filter for delay n = -NP
 *
 *   OUTPUTS:
 *
 *     pswState[0:NP-1]
 *
 *                     updated filter state, ready to filter
 *                     pswInput[siLen], i.e. the next point
 *
 *     pswFiltOut[0:S_LEN-1]
 *
 *                     the filtered output
 *                     same format as pswInput, pswFiltOut[0] is oldest point
 *
 *   RETURN VALUE:
 *
 *     none
 *
 *   DESCRIPTION:
 *
 *     because of the default sign of the coefficients the
 *     formula for the filter is :
 *     i=0, i < S_LEN
 *        out[i] = rounded(state[i]*coef[0])
 *        j=1, j < NP
 *           out[i] -= state[j]*coef[j] (state is taken from either input
 *                                       state[] or prior out[] arrays)
 *        rescale(out[i])
 *        out[i] += in[i]
 *     update final state array using out[]
 *
 *   REFERENCES: Sub-clause 4.1.7 and 4.2.4 of GSM Recomendation 06.20
 *
 *   KEYWORDS: lpc, directform, iir, synthesisfilter, lpcFilt
 *   KEYWORDS: dirForm, dir_mod, dir_clr, dir_neg, dir_set
 *
 *************************************************************************/

void   lpcIir(Shortword pswInput[], Shortword pswCoef[],
                     Shortword pswState[], Shortword pswFiltOut[])
{

/*_________________________________________________________________________
 |                                                                         |
 |                            Automatic Variables                          |
 |_________________________________________________________________________|
*/

  Longword L_Sum;
  short int siStage,
         siSmp;

/*_________________________________________________________________________
 |                                                                         |
 |                              Executable Code                            |
 |_________________________________________________________________________|
*/

  /* filter 1st sample */
  /* ----------------- */

  /* sum past state outputs */
  /* ---------------------- */
  /* 0th coef, with rounding */
  L_Sum = L_msu(LPC_ROUND, pswState[0], pswCoef[0]);

  for (siStage = 1; siStage < NP; siStage++)
  {                                    /* remaining coefs */
    L_Sum = L_msu(L_Sum, pswState[siStage], pswCoef[siStage]);
  }

  /* add input to partial output */
  /* --------------------------- */

  L_Sum = L_shl(L_Sum, ASHIFT);
  L_Sum = L_msu(L_Sum, pswInput[0], 0x8000);

  /* save 1st output sample */
  /* ---------------------- */

  pswFiltOut[0] = extract_h(L_Sum);

  /* filter remaining samples */
  /* ------------------------ */

  for (siSmp = 1; siSmp < S_LEN; siSmp++)
  {

    /* sum past outputs */
    /* ---------------- */
    /* 0th coef, with rounding */
    L_Sum = L_msu(LPC_ROUND, pswFiltOut[siSmp - 1], pswCoef[0]);
    /* remaining coefs */
    for (siStage = 1; ((0 < (siSmp - siStage)) && siStage < NP); siStage++)
    {
      L_Sum = L_msu(L_Sum, pswFiltOut[siSmp - siStage - 1],
                    pswCoef[siStage]);
    }

    /* sum past states, if any */
    /* ----------------------- */

    for (siStage = siSmp; siStage < NP; siStage++)
    {
      L_Sum = L_msu(L_Sum, pswState[siStage - siSmp], pswCoef[siStage]);
    }

    /* add input to partial output */
    /* --------------------------- */

    L_Sum = L_shl(L_Sum, ASHIFT);
    L_Sum = L_msu(L_Sum, pswInput[siSmp], 0x8000);

    /* save current output sample */
    /* -------------------------- */

    pswFiltOut[siSmp] = extract_h(L_Sum);
  }

  /* save final state */
  /* ---------------- */

  for (siStage = 0; siStage < NP; siStage++)
  {
    pswState[siStage] = pswFiltOut[S_LEN - siStage - 1];
  }
}

/***************************************************************************
 *
 *   FUNCTION NAME: lpcIrZsIir
 *
 *   PURPOSE:
 *
 *     The purpose of this function is to calculate the impulse response
 *     via direct form IIR filtering with zero state assuming a NP order
 *     filter and given coefficients
 *
 *   INPUTS:
 *
 *     NP
 *                     order of the lpc filter
 *
 *     S_LEN
 *                     number of samples to filter
 *
 *     pswCoef[0:NP-1]
 *                     array of direct form coefficients
 *                     pswCoef[0] = coeff for delay n = -1
 *                     pswCoef[NP-1] = coeff for delay n = -NP
 *
 *     ASHIFT
 *                     number of shifts input A's have been shifted down by
 *
 *     LPC_ROUND
 *                     rounding constant
 *
 *   OUTPUTS:
 *
 *     pswFiltOut[0:S_LEN-1]
 *
 *                     the filtered output
 *                     same format as pswInput, pswFiltOut[0] is oldest point
 *
 *   RETURN VALUE:
 *
 *     none
 *
 *   DESCRIPTION:
 *
 *     This routine is called by getNWCoefs().
 *
 *     Because of the default sign of the coefficients the
 *     formula for the filter is :
 *     i=0, i < S_LEN
 *        out[i] = rounded(state[i]*coef[0])
 *        j=1, j < NP
 *           out[i] -= state[j]*coef[j] (state taken from prior output[])
 *        rescale(out[i])
 *
 *   REFERENCES: Sub-clause 4.1.8 of GSM Recomendation 06.20
 *
 *   KEYWORDS: lpc, directform, iir, synthesisfilter, lpcFilt
 *   KEYWORDS: dirForm, dir_mod, dir_clr, dir_neg, dir_set
 *
 *************************************************************************/

void   lpcIrZsIir(Shortword pswCoef[], Shortword pswFiltOut[])
{

/*_________________________________________________________________________
 |                                                                         |
 |                            Automatic Variables                          |
 |_________________________________________________________________________|
*/

  Longword L_Sum;
  short int siStage,
         siSmp;

/*_________________________________________________________________________
 |                                                                         |
 |                              Executable Code                            |
 |_________________________________________________________________________|
*/

  /* output 1st sample */
  /* ----------------- */

  pswFiltOut[0] = 0x0400;

  /* filter remaining samples */
  /* ------------------------ */

  for (siSmp = 1; siSmp < S_LEN; siSmp++)
  {

    /* sum past outputs */
    /* ---------------- */
    /* 0th coef, with rounding */
    L_Sum = L_msu(LPC_ROUND, pswFiltOut[siSmp - 1], pswCoef[0]);
    /* remaining coefs */
    for (siStage = 1; ((0 < (siSmp - siStage)) && siStage < NP); siStage++)
    {
      L_Sum = L_msu(L_Sum, pswFiltOut[siSmp - siStage - 1],
                    pswCoef[siStage]);
    }

    /* scale output */
    /* ------------ */

    L_Sum = L_shl(L_Sum, ASHIFT);

    /* save current output sample */
    /* -------------------------- */

    pswFiltOut[siSmp] = extract_h(L_Sum);
  }
}

/***************************************************************************
 *
 *   FUNCTION NAME: lpcZiIir
 *
 *   PURPOSE:
 *     The purpose of this function is to perform direct form iir filtering
 *     with zero input assuming a NP order filter, and given state and
 *     coefficients
 *
 *   INPUTS:
 *
 *     NP
 *                     order of the lpc filter
 *
 *     S_LEN
 *                     number of samples to filter MUST be <= MAX_ZIS
 *
 *     pswCoef[0:NP-1]
 *
 *                     array of direct form coefficients.
 *                     pswCoef[0] = coeff for delay n = -1
 *                     pswCoef[NP-1] = coeff for delay n = -NP
 *
 *     ASHIFT
 *                     number of shifts input A's have been shifted down by
 *
 *     LPC_ROUND
 *                     rounding constant
 *
 *     pswState[0:NP-1]
 *
 *                     array of the filter state following form of pswCoef
 *                     pswState[0] = state of filter for delay n = -1
 *                     pswState[NP-1] = state of filter for delay n = -NP
 *
 *   OUTPUTS:
 *
 *     pswFiltOut[0:S_LEN-1]
 *
 *                     the filtered output
 *                     same format as pswIn, pswFiltOut[0] is oldest point
 *
 *   RETURN VALUE:
 *
 *     none
 *
 *   DESCRIPTION:
 *
 *     The routine is called from sfrmAnalysis, and is used to let the
 *     LPC filters ring out.
 *
 *     because of the default sign of the coefficients the
 *     formula for the filter is :
 *     i=0, i < S_LEN
 *        out[i] = rounded(state[i]*coef[0])
 *        j=1, j < NP
 *           out[i] -= state[j]*coef[j] (state is taken from either input
 *                                       state[] or prior output[] arrays)
 *        rescale(out[i])
 *
 *   REFERENCES: Sub-clause 4.1.7 of GSM Recomendation 06.20
 *
 *   KEYWORDS: lpc, directform, iir, synthesisfilter, lpcFilt
 *   KEYWORDS: dirForm, dir_mod, dir_clr, dir_neg, dir_set
 *
 *************************************************************************/

void   lpcZiIir(Shortword pswCoef[], Shortword pswState[],
                       Shortword pswFiltOut[])
{

/*_________________________________________________________________________
 |                                                                         |
 |                            Automatic Variables                          |
 |_________________________________________________________________________|
*/

  Longword L_Sum;
  short int siStage,
         siSmp;

/*_________________________________________________________________________
 |                                                                         |
 |                              Executable Code                            |
 |_________________________________________________________________________|
*/

  /* filter 1st sample */
  /* ----------------- */

  /* sum past state outputs */
  /* ---------------------- */
  /* 0th coef, with rounding */
  L_Sum = L_msu(LPC_ROUND, pswState[0], pswCoef[0]);

  for (siStage = 1; siStage < NP; siStage++)
  {                                    /* remaining coefs */
    L_Sum = L_msu(L_Sum, pswState[siStage], pswCoef[siStage]);
  }

  /* scale output */
  /* ------------ */

  L_Sum = L_shl(L_Sum, ASHIFT);

  /* save 1st output sample */
  /* ---------------------- */

  pswFiltOut[0] = extract_h(L_Sum);

  /* filter remaining samples */
  /* ------------------------ */

  for (siSmp = 1; siSmp < S_LEN; siSmp++)
  {

    /* sum past outputs */
    /* ---------------- */
    /* 0th coef, with rounding */
    L_Sum = L_msu(LPC_ROUND, pswFiltOut[siSmp - 1], pswCoef[0]);
    /* remaining coefs */
    for (siStage = 1; ((0 < (siSmp - siStage)) && siStage < NP); siStage++)
    {
      L_Sum = L_msu(L_Sum, pswFiltOut[siSmp - siStage - 1],
                    pswCoef[siStage]);
    }

    /* sum past states, if any */
    /* ----------------------- */

    for (siStage = siSmp; siStage < NP; siStage++)
    {
      L_Sum = L_msu(L_Sum, pswState[siStage - siSmp], pswCoef[siStage]);
    }

    /* scale output */
    /* ------------ */

    L_Sum = L_shl(L_Sum, ASHIFT);

    /* save current output sample */
    /* -------------------------- */

    pswFiltOut[siSmp] = extract_h(L_Sum);
  }
}

/***************************************************************************
 *
 *   FUNCTION NAME: lpcZsFir
 *
 *   PURPOSE:
 *     The purpose of this function is to perform direct form fir filtering
 *     with zero state, assuming a NP order filter and given coefficients
 *     and non-zero input.
 *
 *   INPUTS:
 *
 *     NP
 *                     order of the lpc filter
 *
 *     S_LEN
 *                     number of samples to filter
 *
 *     pswInput[0:S_LEN-1]
 *
 *                     input array of points to be filtered.
 *                     pswInput[0] is the oldest point (first to be filtered)
 *                     pswInput[siLen-1] is the last point filtered (newest)
 *
 *     pswCoef[0:NP-1]
 *
 *                     array of direct form coefficients
 *                     pswCoef[0] = coeff for delay n = -1
 *                     pswCoef[NP-1] = coeff for delay n = -NP
 *
 *     ASHIFT
 *                     number of shifts input A's have been shifted down by
 *
 *     LPC_ROUND
 *                     rounding constant
 *
 *   OUTPUTS:
 *
 *     pswFiltOut[0:S_LEN-1]
 *
 *                     the filtered output
 *                     same format as pswInput, pswFiltOut[0] is oldest point
 *
 *   RETURN VALUE:
 *
 *     none
 *
 *   DESCRIPTION:
 *
 *     This routine is used in getNWCoefs().  See section 4.1.7.
 *
 *     because of the default sign of the coefficients the
 *     formula for the filter is :
 *     i=0, i < S_LEN
 *        out[i] = rounded(state[i]*coef[0])
 *        j=1, j < NP
 *           out[i] += state[j]*coef[j] (state taken from in[])
 *        rescale(out[i])
 *        out[i] += in[i]
 *
 *   REFERENCES: Sub-clause 4.1.7 of GSM Recomendation 06.20
 *
 *   KEYWORDS: lpc, directform, fir, lpcFir, inversefilter, lpcFilt
 *   KEYWORDS: dirForm, dir_mod, dir_clr, dir_neg, dir_set, i_dir_mod
 *
 *************************************************************************/

void   lpcZsFir(Shortword pswInput[], Shortword pswCoef[],
                       Shortword pswFiltOut[])
{

/*_________________________________________________________________________
 |                                                                         |
 |                            Automatic Variables                          |
 |_________________________________________________________________________|
*/

  Longword L_Sum;
  short int siStage,
         siSmp;

/*_________________________________________________________________________
 |                                                                         |
 |                              Executable Code                            |
 |_________________________________________________________________________|
*/

  /* output 1st sample */
  /* ----------------- */

  pswFiltOut[0] = pswInput[0];

  /* filter remaining samples */
  /* ------------------------ */

  for (siSmp = 1; siSmp < S_LEN; siSmp++)
  {

    /* sum past outputs */
    /* ---------------- */
    /* 0th coef, with rounding */
    L_Sum = L_mac(LPC_ROUND, pswInput[siSmp - 1], pswCoef[0]);
    /* remaining coefs */
    for (siStage = 1; ((0 < (siSmp - siStage)) && siStage < NP); siStage++)
    {
      L_Sum = L_mac(L_Sum, pswInput[siSmp - siStage - 1],
                    pswCoef[siStage]);
    }

    /* add input to partial output */
    /* --------------------------- */

    L_Sum = L_shl(L_Sum, ASHIFT);
    L_Sum = L_msu(L_Sum, pswInput[siSmp], 0x8000);

    /* save current output sample */
    /* -------------------------- */

    pswFiltOut[siSmp] = extract_h(L_Sum);
  }
}

/***************************************************************************
 *
 *   FUNCTION NAME: lpcZsIir
 *
 *   PURPOSE:
 *
 *     The purpose of this function is to perform direct form IIR filtering
 *     with zero state, assuming a NP order filter and given coefficients
 *     and non-zero input.
 *
 *   INPUTS:
 *
 *     NP
 *                     order of the lpc filter
 *
 *     S_LEN
 *                     number of samples to filter
 *
 *     pswInput[0:S_LEN-1]
 *
 *                     input array of points to be filtered
 *                     pswInput[0] is the oldest point (first to be filtered)
 *                     pswInput[siLen-1] is the last point filtered (newest)
 *
 *     pswCoef[0:NP-1]
 *                     array of direct form coefficients
 *                     pswCoef[0] = coeff for delay n = -1
 *                     pswCoef[NP-1] = coeff for delay n = -NP
 *
 *     ASHIFT
 *                     number of shifts input A's have been shifted down by
 *
 *     LPC_ROUND
 *                     rounding constant
 *
 *   OUTPUTS:
 *
 *     pswFiltOut[0:S_LEN-1]
 *
 *                     the filtered output
 *                     same format as pswInput, pswFiltOut[0] is oldest point
 *
 *   RETURN VALUE:
 *
 *     none
 *
 *   DESCRIPTION:
 *
 *     This routine is used in the subframe analysis process.  It is
 *     called by sfrmAnalysis() and fnClosedLoop().  It is this function
 *     which performs the weighting of the excitation vectors.
 *
 *     because of the default sign of the coefficients the
 *     formula for the filter is :
 *     i=0, i < S_LEN
 *        out[i] = rounded(state[i]*coef[0])
 *        j=1, j < NP
 *           out[i] -= state[j]*coef[j] (state taken from prior out[])
 *        rescale(out[i])
 *        out[i] += in[i]
 *
 *   REFERENCES: Sub-clause 4.1.8.5 of GSM Recomendation 06.20
 *
 *   KEYWORDS: lpc, directform, iir, synthesisfilter, lpcFilt
 *   KEYWORDS: dirForm, dir_mod, dir_clr, dir_neg, dir_set
 *
 *************************************************************************/

void   lpcZsIir(Shortword pswInput[], Shortword pswCoef[],
                       Shortword pswFiltOut[])
{

/*_________________________________________________________________________
 |                                                                         |
 |                            Automatic Variables                          |
 |_________________________________________________________________________|
*/

  Longword L_Sum;
  short int siStage,
         siSmp;

/*_________________________________________________________________________
 |                                                                         |
 |                              Executable Code                            |
 |_________________________________________________________________________|
*/

  /* output 1st sample */
  /* ----------------- */

  pswFiltOut[0] = pswInput[0];

  /* filter remaining samples */
  /* ------------------------ */

  for (siSmp = 1; siSmp < S_LEN; siSmp++)
  {

    /* sum past outputs */
    /* ---------------- */
    /* 0th coef, with rounding */
    L_Sum = L_msu(LPC_ROUND, pswFiltOut[siSmp - 1], pswCoef[0]);
    /* remaining coefs */
    for (siStage = 1; ((0 < (siSmp - siStage)) && siStage < NP); siStage++)
    {
      L_Sum = L_msu(L_Sum, pswFiltOut[siSmp - siStage - 1],
                    pswCoef[siStage]);
    }

    /* add input to partial output */
    /* --------------------------- */

    L_Sum = L_shl(L_Sum, ASHIFT);
    L_Sum = L_msu(L_Sum, pswInput[siSmp], 0x8000);

    /* save current output sample */
    /* -------------------------- */

    pswFiltOut[siSmp] = extract_h(L_Sum);
  }
}

/***************************************************************************
 *
 *   FUNCTION NAME: lpcZsIirP
 *
 *   PURPOSE:
 *
 *     The purpose of this function is to perform direct form iir filtering
 *     with zero state, assuming a NP order filter and given coefficients
 *     and input
 *
 *   INPUTS:
 *
 *     NP
 *                     order of the lpc filter
 *
 *     S_LEN
 *                     number of samples to filter
 *
 *     pswCommonIO[0:S_LEN-1]
 *
 *                     input array of points to be filtered
 *                     pswCommonIO[0] is oldest point (first to be filtered)
 *                     pswCommonIO[siLen-1] is last point filtered (newest)
 *
 *     pswCoef[0:NP-1]
 *                     array of direct form coefficients
 *                     pswCoef[0] = coeff for delay n = -1
 *                     pswCoef[NP-1] = coeff for delay n = -NP
 *
 *     ASHIFT
 *                     number of shifts input A's have been shifted down by
 *
 *     LPC_ROUND
 *                     rounding constant
 *
 *   OUTPUTS:
 *
 *     pswCommonIO[0:S_LEN-1]
 *
 *                     the filtered output
 *                     pswCommonIO[0] is oldest point
 *
 *   RETURN VALUE:
 *
 *     none
 *
 *   DESCRIPTION:
 *
 *     This function is called by geNWCoefs().  See section 4.1.7.
 *
 *     because of the default sign of the coefficients the
 *     formula for the filter is :
 *     i=0, i < S_LEN
 *        out[i] = rounded(state[i]*coef[0])
 *        j=1, j < NP
 *           out[i] += state[j]*coef[j] (state taken from prior out[])
 *        rescale(out[i])
 *        out[i] += in[i]
 *
 *   REFERENCES: Sub-clause 4.1.7 of GSM Recomendation 06.20
 *
 *   KEYWORDS: lpc, directform, iir, synthesisfilter, lpcFilt
 *   KEYWORDS: dirForm, dir_mod, dir_clr, dir_neg, dir_set
 *
 *************************************************************************/

void   lpcZsIirP(Shortword pswCommonIO[], Shortword pswCoef[])
{

/*_________________________________________________________________________
 |                                                                         |
 |                            Automatic Variables                          |
 |_________________________________________________________________________|
*/

  Longword L_Sum;
  short int siStage,
         siSmp;

/*_________________________________________________________________________
 |                                                                         |
 |                              Executable Code                            |
 |_________________________________________________________________________|
*/

  /* filter remaining samples */
  /* ------------------------ */

  for (siSmp = 1; siSmp < S_LEN; siSmp++)
  {

    /* sum past outputs */
    /* ---------------- */
    /* 0th coef, with rounding */
    L_Sum = L_mac(LPC_ROUND, pswCommonIO[siSmp - 1], pswCoef[0]);
    /* remaining coefs */
    for (siStage = 1; ((0 < (siSmp - siStage)) && siStage < NP); siStage++)
    {
      L_Sum = L_mac(L_Sum, pswCommonIO[siSmp - siStage - 1],
                    pswCoef[siStage]);
    }

    /* add input to partial output */
    /* --------------------------- */

    L_Sum = L_shl(L_Sum, ASHIFT);
    L_Sum = L_msu(L_Sum, pswCommonIO[siSmp], 0x8000);

    /* save current output sample */
    /* -------------------------- */

    pswCommonIO[siSmp] = extract_h(L_Sum);
  }
}

/**************************************************************************
 *
 *   FUNCTION NAME: pitchPreFilt
 *
 *   PURPOSE:
 *
 *     Performs pitch pre-filter on excitation in speech decoder.
 *
 *   INPUTS:
 *
 *     pswExcite[0:39]
 *
 *                     Synthetic residual signal to be filtered, a subframe-
 *                     length vector.
 *
 *     ppsrPVecIntFilt[0:9][0:5] ([tap][phase])
 *
 *                     Interpolation filter coefficients.
 *
 *     ppsrSqtrP0[0:2][0:31] ([voicing level-1][gain code])
 *
 *                     Sqrt(P0) look-up table, used to determine pitch
 *                     pre-filtering coefficient.
 *
 *     swRxGsp0
 *
 *                     Coded value from gain quantizer, used to look up
 *                     sqrt(P0).
 *
 *     swRxLag
 *
 *                     Full-resolution lag value (fractional lag *
 *                     oversampling factor), used to index pitch pre-filter
 *                     state.
 *
 *     swUvCode
 *
 *                     Coded voicing level, used to distinguish between
 *                     voiced and unvoiced conditions, and to look up
 *                     sqrt(P0).
 *
 *     swSemiBeta
 *
 *                     The gain applied to the adaptive codebook excitation
 *                     (long-term predictor excitation) limited to a maximum
 *                     of 1.0, used to determine the pitch pre-filter
 *                     coefficient.
 *
 *     snsSqrtRs
 *
 *                     The estimate of the energy in the residual, used only
 *                     for scaling.
 *
 *   OUTPUTS:
 *
 *     pswExciteOut[0:39]
 *
 *                     The output pitch pre-filtered excitation.
 *
 *     pswPPreState[0:44]
 *
 *                     Contains the state of the pitch pre-filter
 *
 *   RETURN VALUE:
 *
 *     none
 *
 *   DESCRIPTION:
 *
 *     If the voicing mode for the frame is unvoiced, then the pitch pre-
 *     filter state is updated with the input excitation, and the input
 *     excitation is copied to the output.
 *
 *     If voiced: first the energy in the input excitation is calculated.
 *     Then, the coefficient of the pitch pre-filter is obtained:
 *
 *     PpfCoef = POST_EPSILON * min(beta, sqrt(P0)).
 *
 *     Then, the pitch pre-filter is performed:
 *
 *     ex_p(n) = ex(n)  +  PpfCoef * ex_p(n-L)
 *
 *     The ex_p(n-L) sample is interpolated from the surrounding samples,
 *     even for integer values of L.
 *
 *     Note: The coefficients of the interpolating filter are multiplied
 *     by PpfCoef, rather multiplying ex_p(n_L) after interpolation.
 *
 *     Finally, the energy in the output excitation is calculated, and
 *     automatic gain control is applied to the output signal so that
 *     its energy matches the original.
 *
 *     The pitch pre-filter is described in section 4.2.2.
 *
 *   REFERENCES: Sub-clause 4.2.2 of GSM Recomendation 06.20
 *
 *   KEYWORDS: prefilter, pitch, pitchprefilter, excitation, residual
 *
 *************************************************************************/

static void pitchPreFilt(Shortword pswExcite[],
                                Shortword swRxGsp0,
                                Shortword swRxLag, Shortword swUvCode,
                              Shortword swSemiBeta, struct NormSw snsSqrtRs,
                                Shortword pswExciteOut[],
                                Shortword pswPPreState[])
{

/*_________________________________________________________________________
 |                                                                         |
 |                              Local Constants                            |
 |_________________________________________________________________________|
*/

#define  POST_EPSILON  0x2666

/*_________________________________________________________________________
 |                                                                         |
 |                           Local Static Variables                        |
 |_________________________________________________________________________|
*/


/*_________________________________________________________________________
 |                                                                         |
 |                            Automatic Variables                          |
 |_________________________________________________________________________|
*/

  Longword L_1,
         L_OrigEnergy;

  Shortword swScale,
         swSqrtP0,
         swIntLag,
         swRemain,
         swEnergy,
         pswInterpCoefs[P_INT_MACS];

  short int i,
         j;

  struct NormSw snsOrigEnergy;

  Shortword *pswPPreCurr = &pswPPreState[LTP_LEN];

/*_________________________________________________________________________
 |                                                                         |
 |                              Executable Code                            |
 |_________________________________________________________________________|
*/

  /* Initialization */
  /*----------------*/

  swEnergy = 0;

  /* Check voicing level */
  /*---------------------*/

  if (swUvCode == 0)
  {

    /* Unvoiced: perform one subframe of delay on state, copy input to */
    /* state, copy input to output (if not same)                       */
    /*-----------------------------------------------------------------*/

    for (i = 0; i < LTP_LEN - S_LEN; i++)
      pswPPreState[i] = pswPPreState[i + S_LEN];

    for (i = 0; i < S_LEN; i++)
      pswPPreState[i + LTP_LEN - S_LEN] = pswExcite[i];

    if (pswExciteOut != pswExcite)
    {

      for (i = 0; i < S_LEN; i++)
        pswExciteOut[i] = pswExcite[i];
    }
  }
  else
  {

    /* Voiced: calculate energy in input, filter, calculate energy in */
    /* output, scale                                                  */
    /*----------------------------------------------------------------*/

    /* Get energy in input excitation vector */
    /*---------------------------------------*/

    swEnergy = add(negate(shl(snsSqrtRs.sh, 1)), 3);

    if (swEnergy > 0)
    {

      /* High-energy residual: scale input vector during energy */
      /* calculation.  The shift count + 1 of the energy of the */
      /* residual estimate is used as an estimate of the shift  */
      /* count needed for the excitation energy                 */
      /*--------------------------------------------------------*/


      snsOrigEnergy.sh = g_corr1s(pswExcite, swEnergy, &L_OrigEnergy);
      snsOrigEnergy.man = round(L_OrigEnergy);

    }
    else
    {

      /* set shift count to zero for AGC later */
      /*---------------------------------------*/

      swEnergy = 0;

      /* Lower-energy residual: no overflow protection needed */
      /*------------------------------------------------------*/

      L_OrigEnergy = 0;
      for (i = 0; i < S_LEN; i++)
      {

        L_OrigEnergy = L_mac(L_OrigEnergy, pswExcite[i], pswExcite[i]);
      }

      snsOrigEnergy.sh = norm_l(L_OrigEnergy);
      snsOrigEnergy.man = round(L_shl(L_OrigEnergy, snsOrigEnergy.sh));
    }

    /* Determine pitch pre-filter coefficient, and scale the appropriate */
    /* phase of the interpolating filter by it                           */
    /*-------------------------------------------------------------------*/

    swSqrtP0 = ppsrSqrtP0[swUvCode - 1][swRxGsp0];

    if (sub(swSqrtP0, swSemiBeta) > 0)
      swScale = swSemiBeta;
    else
      swScale = swSqrtP0;

    swScale = mult_r(POST_EPSILON, swScale);

    get_ipjj(swRxLag, &swIntLag, &swRemain);

    for (i = 0; i < P_INT_MACS; i++)
      pswInterpCoefs[i] = mult_r(ppsrPVecIntFilt[i][swRemain], swScale);

    /* Perform filter */
    /*----------------*/

    for (i = 0; i < S_LEN; i++)
    {

      L_1 = L_deposit_h(pswExcite[i]);

      for (j = 0; j < P_INT_MACS - 1; j++)
      {

        L_1 = L_mac(L_1, pswPPreCurr[i - swIntLag - P_INT_MACS / 2 + j],
                    pswInterpCoefs[j]);
      }

      pswPPreCurr[i] = mac_r(L_1,
                             pswPPreCurr[i - swIntLag + P_INT_MACS / 2 - 1],
                             pswInterpCoefs[P_INT_MACS - 1]);
    }

    /* Get energy in filtered vector, determine automatic-gain-control */
    /* scale factor                                                    */
    /*-----------------------------------------------------------------*/

    swScale = agcGain(pswPPreCurr, snsOrigEnergy, swEnergy);

    /* Scale filtered vector by AGC, put out.  NOTE: AGC scale returned */
    /* by routine above is divided by two, hence the shift below        */
    /*------------------------------------------------------------------*/

    for (i = 0; i < S_LEN; i++)
    {

      L_1 = L_mult(pswPPreCurr[i], swScale);
      L_1 = L_shl(L_1, 1);
      pswExciteOut[i] = round(L_1);
    }

    /* Update pitch pre-filter state */
    /*-------------------------------*/

    for (i = 0; i < LTP_LEN; i++)
      pswPPreState[i] = pswPPreState[i + S_LEN];
  }
}

/***************************************************************************
 *
 *   FUNCTION NAME: r0BasedEnergyShft
 *
 *   PURPOSE:
 *
 *     Given an R0 voicing level, find the number of shifts to be
 *     performed on the energy to ensure that the subframe energy does
 *     not overflow.  example if energy can maximally take the value
 *     4.0, then 2 shifts are required.
 *
 *   INPUTS:
 *
 *     swR0Index
 *                     R0 codeword (0-0x1f)
 *
 *   OUTPUTS:
 *
 *     none
 *
 *   RETURN VALUE:
 *
 *     swShiftDownSignal
 *
 *                    number of right shifts to apply to energy (0..6)
 *
 *   DESCRIPTION:
 *
 *     Based on the R0, the average frame energy, we can get an
 *     upper bound on the energy any one subframe can take on.
 *     Using this upper bound we can calculate what right shift is
 *     needed to ensure an unsaturated output out of a subframe
 *     energy calculation (g_corr).
 *
 *   REFERENCES: Sub-clause 4.1.9 and 4.2.1 of GSM Recomendation 06.20
 *
 *   KEYWORDS: spectral postfilter
 *
 *************************************************************************/

Shortword r0BasedEnergyShft(Shortword swR0Index)
{

/*_________________________________________________________________________
 |                                                                         |
 |                            Automatic Variables                          |
 |_________________________________________________________________________|
*/

  Shortword swShiftDownSignal;

/*_________________________________________________________________________
 |                                                                         |
 |                              Executable Code                            |
 |_________________________________________________________________________|
*/

  if (sub(swR0Index, 26) <= 0)
  {
    if (sub(swR0Index, 23) <= 0)
    {
      if (sub(swR0Index, 21) <= 0)
        swShiftDownSignal = 0;         /* r0  [0,  21] */
      else
        swShiftDownSignal = 1;         /* r0  [22, 23] */
    }
    else
    {
      if (sub(swR0Index, 24) <= 0)
        swShiftDownSignal = 2;         /* r0  [23, 24] */
      else
        swShiftDownSignal = 3;         /* r0  [24, 26] */
    }
  }
  else
  {                                    /* r0 index > 26 */
    if (sub(swR0Index, 28) <= 0)
    {
      swShiftDownSignal = 4;           /* r0  [26, 28] */
    }
    else
    {
      if (sub(swR0Index, 29) <= 0)
        swShiftDownSignal = 5;         /* r0  [28, 29] */
      else
        swShiftDownSignal = 6;         /* r0  [29, 31] */
    }
  }
  if (sub(swR0Index, 18) > 0)
    swShiftDownSignal = add(swShiftDownSignal, 2);

  return (swShiftDownSignal);
}

/***************************************************************************
 *
 *   FUNCTION NAME: rcToADp
 *
 *   PURPOSE:
 *
 *     This subroutine computes a vector of direct form LPC filter
 *     coefficients, given an input vector of reflection coefficients.
 *     Double precision is used internally, but 16 bit direct form
 *     filter coefficients are returned.
 *
 *   INPUTS:
 *
 *     NP
 *                     order of the LPC filter (global constant)
 *
 *     swAscale
 *                     The multiplier which scales down the direct form
 *                     filter coefficients.
 *
 *     pswRc[0:NP-1]
 *                     The input vector of reflection coefficients.
 *
 *   OUTPUTS:
 *
 *     pswA[0:NP-1]
 *                     Array containing the scaled down direct form LPC
 *                     filter coefficients.
 *
 *   RETURN VALUE:
 *
 *     siLimit
 *                     1 if limiting occured in computation, 0 otherwise.
 *
 *   DESCRIPTION:
 *
 *     This function performs the conversion from reflection coefficients
 *     to direct form LPC filter coefficients. The direct form coefficients
 *     are scaled by multiplication by swAscale. NP, the filter order is 10.
 *     The a's and rc's each have NP elements in them. Double precision
 *     calculations are used internally.
 *
 *        The equations are:
 *        for i = 0 to NP-1{
 *
 *          a(i)(i) = rc(i)              (scaling by swAscale occurs here)
 *
 *          for j = 0 to i-1
 *             a(i)(j) = a(i-1)(j) + rc(i)*a(i-1)(i-j-1)
 *       }
 *
 *     See page 443, of
 *     "Digital Processing of Speech Signals" by L.R. Rabiner and R.W.
 *     Schafer; Prentice-Hall; Englewood Cliffs, NJ (USA).  1978.
 *
 *   REFERENCES: Sub-clause 4.1.7 and 4.2.3 of GSM Recomendation 06.20
 *
 *  KEYWORDS: reflectioncoefficients, parcors, conversion, rctoadp, ks, as
 *  KEYWORDS: parcorcoefficients, lpc, flat, vectorquantization
 *
 *************************************************************************/

short  rcToADp(Shortword swAscale, Shortword pswRc[],
                      Shortword pswA[])
{

/*_________________________________________________________________________
 |                                                                         |
 |                            Automatic Variables                          |
 |_________________________________________________________________________|
*/

  Longword pL_ASpace[NP],
         pL_tmpSpace[NP],
         L_temp,
        *pL_A,
        *pL_tmp,
        *pL_swap;

  short int i,
         j,                            /* loop counters */
         siLimit;

/*_________________________________________________________________________
 |                                                                         |
 |                              Executable Code                            |
 |_________________________________________________________________________|
*/

  /* Initialize starting addresses for temporary buffers */
  /*-----------------------------------------------------*/

  pL_A = pL_ASpace;
  pL_tmp = pL_tmpSpace;

  /* Initialize the flag for checking if limiting has occured */
  /*----------------------------------------------------------*/

  siLimit = 0;

  /* Compute direct form filter coefficients, pswA[0],...,pswA[9] */
  /*-------------------------------------------------------------------*/

  for (i = 0; i < NP; i++)
  {

    pL_tmp[i] = L_mult(swAscale, pswRc[i]);
    for (j = 0; j <= i - 1; j++)
    {
      L_temp = L_mpy_ls(pL_A[i - j - 1], pswRc[i]);
      pL_tmp[j] = L_add(L_temp, pL_A[j]);
      siLimit |= isLwLimit(pL_tmp[j]);
    }
    if (i != NP - 1)
    {
      /* Swap swA and swTmp buffers */

      pL_swap = pL_tmp;
      pL_tmp = pL_A;
      pL_A = pL_swap;
    }
  }

  for (i = 0; i < NP; i++)
  {
    pswA[i] = round(pL_tmp[i]);
    siLimit |= isSwLimit(pswA[i]);
  }
  return (siLimit);
}

/***************************************************************************
 *
 *   FUNCTION NAME: rcToCorrDpL
 *
 *   PURPOSE:
 *
 *     This subroutine computes an autocorrelation vector, given a vector
 *     of reflection coefficients as an input. Double precision calculations
 *     are used internally, and a double precision (Longword)
 *     autocorrelation sequence is returned.
 *
 *   INPUTS:
 *
 *     NP
 *                     LPC filter order passed in as a global constant.
 *
 *     swAshift
 *                     Number of right shifts to be applied to the
 *                     direct form filter coefficients being computed
 *                     as an intermediate step to generating the
 *                     autocorrelation sequence.
 *
 *     swAscale
 *                     A multiplicative scale factor corresponding to
 *                     swAshift; i.e. swAscale = 2 ^(-swAshift).
 *
 *     pswRc[0:NP-1]
 *                     An input vector of reflection coefficients.
 *
 *   OUTPUTS:
 *
 *     pL_R[0:NP]
 *                     An output Longword array containing the
 *                     autocorrelation vector where
 *                     pL_R[0] = 0x7fffffff; (i.e., ~1.0).
 *
 *   RETURN VALUE:
 *
 *     none
 *
 *   DESCRIPTION:
 *
 *     The algorithm used for computing the correlation sequence is
 *     described on page 232 of the book "Linear Prediction of Speech",
 *     by J.D.  Markel and A.H. Gray, Jr.; Springer-Verlag, Berlin,
 *     Heidelberg, New York, 1976.
 *
 *   REFERENCES: Sub_Clause 4.1.4 and 4.2.1  of GSM Recomendation 06.20
 *
 *   KEYWORDS: normalized autocorrelation, reflection coefficients
 *   KEYWORDS: conversion
 *
 **************************************************************************/

void   rcToCorrDpL(Shortword swAshift, Shortword swAscale,
                          Shortword pswRc[], Longword pL_R[])
{

/*_________________________________________________________________________
 |                                                                         |
 |                            Automatic Variables                          |
 |_________________________________________________________________________|
*/

  Longword pL_ASpace[NP],
         pL_tmpSpace[NP],
         L_temp,
         L_sum,
        *pL_A,
        *pL_tmp,
        *pL_swap;

  short int i,
         j;                            /* loop control variables */

/*_________________________________________________________________________
 |                                                                         |
 |                              Executable Code                            |
 |_________________________________________________________________________|
*/

  /* Set R[0] = 0x7fffffff, (i.e., R[0] = 1.0) */
  /*-------------------------------------------*/

  pL_R[0] = LW_MAX;

  /* Assign an address onto each of the two temporary buffers */
  /*----------------------------------------------------------*/

  pL_A = pL_ASpace;
  pL_tmp = pL_tmpSpace;

  /* Compute correlations R[1],...,R[10] */
  /*------------------------------------*/

  for (i = 0; i < NP; i++)
  {

    /* Compute, as an intermediate step, the filter coefficients for */
    /* for an i-th order direct form filter (pL_tmp[j],j=0,i)        */
    /*---------------------------------------------------------------*/

    pL_tmp[i] = L_mult(swAscale, pswRc[i]);
    for (j = 0; j <= i - 1; j++)
    {
      L_temp = L_mpy_ls(pL_A[i - j - 1], pswRc[i]);
      pL_tmp[j] = L_add(L_temp, pL_A[j]);
    }

    /* Swap pL_A and pL_tmp buffers */
    /*------------------------------*/

    pL_swap = pL_A;
    pL_A = pL_tmp;
    pL_tmp = pL_swap;

    /* Given the direct form filter coefficients for an i-th order filter  */
    /* and the autocorrelation vector computed up to and including stage i */
    /* compute the autocorrelation coefficient R[i+1]                      */
    /*---------------------------------------------------------------------*/

    L_temp = L_mpy_ll(pL_A[0], pL_R[i]);
    L_sum = L_negate(L_temp);

    for (j = 1; j <= i; j++)
    {
      L_temp = L_mpy_ll(pL_A[j], pL_R[i - j]);
      L_sum = L_sub(L_sum, L_temp);
    }
    pL_R[i + 1] = L_shl(L_sum, swAshift);

  }
}

/***************************************************************************
 *
 *   FUNCTION NAME: res_eng
 *
 *   PURPOSE:
 *
 *     Calculates square root of subframe residual energy estimate:
 *
 *                     sqrt( R(0)(1-k1**2)...(1-k10**2) )
 *
 *   INPUTS:
 *
 *     pswReflecCoefIn[0:9]
 *
 *                     Array of reflection coeffcients.
 *
 *     swRq
 *
 *                     Subframe energy = sqrt(frame_energy * S_LEN/2**S_SH)
 *                     (quantized).
 *
 *   OUTPUTS:
 *
 *     psnsSqrtRsOut
 *
 *                     (Pointer to) the output residual energy estimate.
 *
 *   RETURN VALUE:
 *
 *     The shift count of the normalized residual energy estimate, as int.
 *
 *   DESCRIPTION:
 *
 *     First, the canonic product of the (1-ki**2) terms is calculated
 *     (normalizations are done to maintain precision).  Also, a factor of
 *     2**S_SH is applied to the product to offset this same factor in the
 *     quantized square root of the subframe energy.
 *
 *     Then the product is square-rooted, and multiplied by the quantized
 *     square root of the subframe energy.  This combined product is put
 *     out as a normalized fraction and shift count (mantissa and exponent).
 *
 *   REFERENCES: Sub-clause 4.1.7 and 4.2.1 of GSM Recomendation 06.20
 *
 *   KEYWORDS: residualenergy, res_eng, rs
 *
 *************************************************************************/

void   res_eng(Shortword pswReflecCoefIn[], Shortword swRq,
                      struct NormSw *psnsSqrtRsOut)
{
/*_________________________________________________________________________
 |                                                                         |
 |                              Local Constants                            |
 |_________________________________________________________________________|
*/

#define  S_SH          6               /* ceiling(log2(S_LEN)) */
#define  MINUS_S_SH    -S_SH


/*_________________________________________________________________________
 |                                                                         |
 |                            Automatic Variables                          |
 |_________________________________________________________________________|
*/

  Longword L_Product,
         L_Shift,
         L_SqrtResEng;

  Shortword swPartialProduct,
         swPartialProductShift,
         swTerm,
         swShift,
         swSqrtPP,
         swSqrtPPShift,
         swSqrtResEng,
         swSqrtResEngShift;

  short int i;

/*_________________________________________________________________________
 |                                                                         |
 |                              Executable Code                            |
 |_________________________________________________________________________|
*/

  /* Form canonic product, maintain precision and shift count */
  /*----------------------------------------------------------*/

  /* (Start off with unity product (actually -1), and shift offset) */
  /*----------------------------------------------------------------*/
  swPartialProduct = SW_MIN;
  swPartialProductShift = MINUS_S_SH;

  for (i = 0; i < NP; i++)
  {

    /* Get next (-1 + k**2) term, form partial canonic product */
    /*---------------------------------------------------------*/


    swTerm = mac_r(LW_MIN, pswReflecCoefIn[i], pswReflecCoefIn[i]);

    L_Product = L_mult(swTerm, swPartialProduct);

    /* Normalize partial product, round */
    /*----------------------------------*/

    swShift = norm_s(extract_h(L_Product));
    swPartialProduct = round(L_shl(L_Product, swShift));
    swPartialProductShift = add(swPartialProductShift, swShift);
  }

  /* Correct sign of product, take square root */
  /*-------------------------------------------*/

  swPartialProduct = abs_s(swPartialProduct);

  swSqrtPP = sqroot(L_deposit_h(swPartialProduct));

  L_Shift = L_shr(L_deposit_h(swPartialProductShift), 1);

  swSqrtPPShift = extract_h(L_Shift);

  if (extract_l(L_Shift) != 0)
  {

    /* Odd exponent: shr above needs to be compensated by multiplying */
    /* mantissa by sqrt(0.5)                                          */
    /*----------------------------------------------------------------*/

    swSqrtPP = mult_r(swSqrtPP, SQRT_ONEHALF);
  }

  /* Form final product, the residual energy estimate, and do final */
  /* normalization                                                  */
  /*----------------------------------------------------------------*/

  L_SqrtResEng = L_mult(swRq, swSqrtPP);

  swShift = norm_l(L_SqrtResEng);
  swSqrtResEng = round(L_shl(L_SqrtResEng, swShift));
  swSqrtResEngShift = add(swSqrtPPShift, swShift);

  /* Return */
  /*--------*/
  psnsSqrtRsOut->man = swSqrtResEng;
  psnsSqrtRsOut->sh = swSqrtResEngShift;

  return;
}

/***************************************************************************
 *
 *   FUNCTION NAME: rs_rr
 *
 *   PURPOSE:
 *
 *     Calculates sqrt(RS/R(x,x)) using floating point format,
 *     where RS is the approximate residual energy in a given
 *     subframe and R(x,x) is the power in each long term
 *     predictor vector or in each codevector.
 *     Used in the joint optimization of the gain and the long
 *     term predictor coefficient.
 *
 *   INPUTS:
 *
 *     pswExcitation[0:39] - excitation signal array
 *     snsSqrtRs - structure sqrt(RS) normalized with mantissa and shift
 *
 *   OUTPUTS:
 *
 *     snsSqrtRsRr - structure sqrt(RS/R(x,x)) with mantissa and shift
 *
 *   RETURN VALUE:
 *
 *     None
 *
 *   DESCRIPTION:
 *
 *     Implemented as sqrt(RS)/sqrt(R(x,x)) where both sqrts
 *     are stored normalized (0.5<=x<1.0) and the associated
 *     shift.  See section 4.1.11.1 for details
 *
 *   REFERENCES: Sub-clause 4.1.11.1 and 4.2.1 of GSM
 *               Recomendation 06.20
 *
 *   KEYWORDS: rs_rr, sqroot
 *
 *************************************************************************/

void   rs_rr(Shortword pswExcitation[], struct NormSw snsSqrtRs,
                    struct NormSw *snsSqrtRsRr)
{

/*_________________________________________________________________________
 |                                                                         |
 |                            Automatic Variables                          |
 |_________________________________________________________________________|
*/
  Longword L_Temp;
  Shortword swTemp,
         swTemp2,
         swEnergy,
         swNormShift,
         swShift;

/*_________________________________________________________________________
 |                                                                         |
 |                              Executable Code                            |
 |_________________________________________________________________________|
*/

  swEnergy = sub(shl(snsSqrtRs.sh, 1), 3);      /* shift*2 + margin ==
                                                 * energy. */


  if (swEnergy < 0)
  {

    /* High-energy residual: scale input vector during energy */
    /* calculation. The shift count of the energy of the      */
    /* residual estimate is used as an estimate of the shift  */
    /* count needed for the excitation energy                 */
    /*--------------------------------------------------------*/

    swNormShift = g_corr1s(pswExcitation, negate(swEnergy), &L_Temp);

  }
  else
  {

    /* Lower-energy residual: no overflow protection needed */
    /*------------------------------------------------------*/

    swNormShift = g_corr1(pswExcitation, &L_Temp);
  }

  /* Compute single precision square root of energy sqrt(R(x,x)) */
  /* ----------------------------------------------------------- */
  swTemp = sqroot(L_Temp);

  /* If odd no. of shifts compensate by sqrt(0.5) */
  /* -------------------------------------------- */
  if (swNormShift & 1)
  {

    /* Decrement no. of shifts in accordance with sqrt(0.5) */
    /* ---------------------------------------------------- */
    swNormShift = sub(swNormShift, 1);

    /* sqrt(R(x,x) = sqrt(R(x,x)) * sqrt(0.5) */
    /* -------------------------------------- */
    L_Temp = L_mult(0x5a82, swTemp);
  }
  else
  {
    L_Temp = L_deposit_h(swTemp);
  }

  /* Normalize again and update shifts */
  /* --------------------------------- */
  swShift = norm_l(L_Temp);
  swNormShift = add(shr(swNormShift, 1), swShift);
  L_Temp = L_shl(L_Temp, swShift);

  /* Shift sqrt(RS) to make sure less than divisor */
  /* --------------------------------------------- */
  swTemp = shr(snsSqrtRs.man, 1);

  /* Divide sqrt(RS)/sqrt(R(x,x)) */
  /* ---------------------------- */
  swTemp2 = divide_s(swTemp, round(L_Temp));

  /* Calculate shift for division, compensate for shift before division */
  /* ------------------------------------------------------------------ */
  swNormShift = sub(snsSqrtRs.sh, swNormShift);
  swNormShift = sub(swNormShift, 1);

  /* Normalize and get no. of shifts */
  /* ------------------------------- */
  swShift = norm_s(swTemp2);
  snsSqrtRsRr->sh = add(swNormShift, swShift);
  snsSqrtRsRr->man = shl(swTemp2, swShift);

}

/***************************************************************************
 *
 *   FUNCTION NAME: rs_rrNs
 *
 *   PURPOSE:
 *
 *     Calculates sqrt(RS/R(x,x)) using floating point format,
 *     where RS is the approximate residual energy in a given
 *     subframe and R(x,x) is the power in each long term
 *     predictor vector or in each codevector.
 *
 *     Used in the joint optimization of the gain and the long
 *     term predictor coefficient.
 *
 *   INPUTS:
 *
 *     pswExcitation[0:39] - excitation signal array
 *     snsSqrtRs - structure sqrt(RS) normalized with mantissa and shift
 *
 *   OUTPUTS:
 *
 *     snsSqrtRsRr - structure sqrt(RS/R(x,x)) with mantissa and shift
 *
 *   RETURN VALUE:
 *
 *     None
 *
 *   DESCRIPTION:
 *
 *     Implemented as sqrt(RS)/sqrt(R(x,x)) where both sqrts
 *     are stored normalized (0.5<=x<1.0) and the associated
 *     shift.
 *
 *   REFERENCES: Sub-clause 4.1.11.1 and 4.2.1 of GSM
 *               Recomendation 06.20
 *
 *   KEYWORDS: rs_rr, sqroot
 *
 *************************************************************************/

void   rs_rrNs(Shortword pswExcitation[], struct NormSw snsSqrtRs,
                      struct NormSw *snsSqrtRsRr)
{

/*_________________________________________________________________________
 |                                                                         |
 |                            Automatic Variables                          |
 |_________________________________________________________________________|
*/
  Longword L_Temp;
  Shortword swTemp,
         swTemp2,
         swNormShift,
         swShift;

/*_________________________________________________________________________
 |                                                                         |
 |                              Executable Code                            |
 |_________________________________________________________________________|
*/

  /* Lower-energy residual: no overflow protection needed */
  /*------------------------------------------------------*/

  swNormShift = g_corr1(pswExcitation, &L_Temp);


  /* Compute single precision square root of energy sqrt(R(x,x)) */
  /* ----------------------------------------------------------- */
  swTemp = sqroot(L_Temp);

  /* If odd no. of shifts compensate by sqrt(0.5) */
  /* -------------------------------------------- */
  if (swNormShift & 1)
  {

    /* Decrement no. of shifts in accordance with sqrt(0.5) */
    /* ---------------------------------------------------- */
    swNormShift = sub(swNormShift, 1);

    /* sqrt(R(x,x) = sqrt(R(x,x)) * sqrt(0.5) */
    /* -------------------------------------- */
    L_Temp = L_mult(0x5a82, swTemp);
  }
  else
  {
    L_Temp = L_deposit_h(swTemp);
  }

  /* Normalize again and update shifts */
  /* --------------------------------- */

  swShift = norm_l(L_Temp);
  swNormShift = add(shr(swNormShift, 1), swShift);
  L_Temp = L_shl(L_Temp, swShift);

  /* Shift sqrt(RS) to make sure less than divisor */
  /* --------------------------------------------- */
  swTemp = shr(snsSqrtRs.man, 1);

  /* Divide sqrt(RS)/sqrt(R(x,x)) */
  /* ---------------------------- */
  swTemp2 = divide_s(swTemp, round(L_Temp));

  /* Calculate shift for division, compensate for shift before division */
  /* ------------------------------------------------------------------ */
  swNormShift = sub(snsSqrtRs.sh, swNormShift);
  swNormShift = sub(swNormShift, 1);

  /* Normalize and get no. of shifts */
  /* ------------------------------- */
  swShift = norm_s(swTemp2);
  snsSqrtRsRr->sh = add(swNormShift, swShift);
  snsSqrtRsRr->man = shl(swTemp2, swShift);

}


/***************************************************************************
 *
 *   FUNCTION NAME: scaleExcite
 *
 *   PURPOSE:
 *
 *     Scale an arbitrary excitation vector (codevector or
 *     pitch vector)
 *
 *   INPUTS:
 *
 *     pswVect[0:39] - the unscaled vector.
 *     iGsp0Scale - an positive offset to compensate for the fact
 *                  that GSP0 table is scaled down.
 *     swErrTerm - rather than a gain being passed in, (beta, gamma)
 *                 it is calculated from this error term - either
 *                 Gsp0[][][0] error term A or Gsp0[][][1] error
 *                 term B. Beta is calculated from error term A,
 *                 gamma from error term B.
 *     snsRS - the RS_xx appropriate to pswVect.
 *
 *   OUTPUTS:
 *
 *      pswScldVect[0:39] - the output, scaled excitation vector.
 *
 *   RETURN VALUE:
 *
 *     swGain - One of two things.  Either a clamped value of 0x7fff if the
 *              gain's shift was > 0 or the rounded vector gain otherwise.
 *
 *   DESCRIPTION:
 *
 *     If gain > 1.0 then
 *         (do not shift gain up yet)
 *         partially scale vector element THEN shift and round save
 *      else
 *         shift gain correctly
 *         scale vector element and round save
 *         update state array
 *
 *   REFERENCES: Sub-clause 4.1.10.2 and 4.2.1 of GSM
 *               Recomendation 06.20
 *
 *   KEYWORDS: excite_vl, sc_ex, excitevl, scaleexcite, codevector, p_vec,
 *   KEYWORDS: x_vec, pitchvector, gain, gsp0
 *
 *************************************************************************/

Shortword scaleExcite(Shortword pswVect[],
                             Shortword swErrTerm, struct NormSw snsRS,
                             Shortword pswScldVect[])
{

/*_________________________________________________________________________
 |                                                                         |
 |                            Automatic Variables                          |
 |_________________________________________________________________________|
*/
  Longword L_GainUs,
         L_scaled,
         L_Round_off;
  Shortword swGain,
         swGainUs,
         swGainShift,
         i,
         swGainUsShft;

/*_________________________________________________________________________
 |                                                                         |
 |                              Executable Code                            |
 |_________________________________________________________________________|
*/


  L_GainUs = L_mult(swErrTerm, snsRS.man);
  swGainUsShft = norm_s(extract_h(L_GainUs));
  L_GainUs = L_shl(L_GainUs, swGainUsShft);

  swGainShift = add(swGainUsShft, snsRS.sh);
  swGainShift = sub(swGainShift, GSP0_SCALE);


  /* gain > 1.0 */
  /* ---------- */

  if (swGainShift < 0)
  {
    swGainUs = round(L_GainUs);

    L_Round_off = L_shl((long) 32768, swGainShift);

    for (i = 0; i < S_LEN; i++)
    {
      L_scaled = L_mac(L_Round_off, swGainUs, pswVect[i]);
      L_scaled = L_shr(L_scaled, swGainShift);
      pswScldVect[i] = extract_h(L_scaled);
    }

    if (swGainShift == 0)
      swGain = swGainUs;
    else
      swGain = 0x7fff;
  }

  /* gain < 1.0 */
  /* ---------- */

  else
  {

    /* shift down or not at all */
    /* ------------------------ */
    if (swGainShift > 0)
      L_GainUs = L_shr(L_GainUs, swGainShift);

    /* the rounded actual vector gain */
    /* ------------------------------ */
    swGain = round(L_GainUs);

    /* now scale the vector (with rounding) */
    /* ------------------------------------ */

    for (i = 0; i < S_LEN; i++)
    {
      L_scaled = L_mac((long) 32768, swGain, pswVect[i]);
      pswScldVect[i] = extract_h(L_scaled);
    }
  }
  return (swGain);
}

/***************************************************************************
 *
 *   FUNCTION NAME: spectralPostFilter
 *
 *   PURPOSE:
 *
 *     Perform spectral post filter on the output of the
 *     synthesis filter.
 *
 *
 *   INPUT:
 *
 *      S_LEN         a global constant
 *
 *      pswSPFIn[0:S_LEN-1]
 *
 *                    input to the routine. Unmodified
 *                    pswSPFIn[0] is the oldest point (first to be filtered),
 *                    pswSPFIn[iLen-1] is the last pointer filtered,
 *                    the newest.
 *
 *      pswNumCoef[0:NP-1],pswDenomCoef[0:NP-1]
 *
 *                    numerator and denominator
 *                    direct form coeffs used by postfilter.
 *                    Exactly like lpc coefficients in format.  Shifted down
 *                    by iAShift to ensure that they are < 1.0.
 *
 *      gpswPostFiltStateNum[0:NP-1], gpswPostFiltStateDenom[0:NP-1]
 *
 *                    array of the filter state.
 *                    Same format as coefficients: *praState = state of
 *                    filter for delay n = -1 praState[NP] = state of
 *                    filter for delay n = -NP These numbers are not
 *                    shifted at all.  These states are static to this
 *                    file.
 *
 *   OUTPUT:
 *
 *      gpswPostFiltStateNum[0:NP-1], gpswPostFiltStateDenom[0:NP-1]
 *
 *                      See above for description.  These are updated.
 *
 *      pswSPFOut[0:S_LEN-1]
 *
 *                    same format as pswSPFIn,
 *                    *pswSPFOut is oldest point. The filtered output.
 *                    Note this routine can handle pswSPFOut = pswSPFIn.
 *                    output can be the same as input.  i.e. in place
 *                    calculation.
 *
 *   RETURN:
 *
 *      none
 *
 *   DESCRIPTION:
 *
 *      find energy in input,
 *      perform the numerator fir
 *      perform the denominator iir
 *      perform the post emphasis
 *      find energy in signal,
 *      perform the agc using energy in and energy in signam after
 *      post emphasis signal
 *
 *      The spectral postfilter is described in section 4.2.4.
 *
 *   REFERENCES: Sub-clause 4.2.4 of GSM Recomendation 06.20
 *
 *   KEYWORDS: postfilter, emphasis, postemphasis, brightness,
 *   KEYWORDS: numerator, deminator, filtering, lpc,
 *
 *************************************************************************/

static void spectralPostFilter(Shortword pswSPFIn[],
                                      Shortword pswNumCoef[],
                            Shortword pswDenomCoef[], Shortword pswSPFOut[])
{
/*_________________________________________________________________________
 |                                                                         |
 |                              Local Constants                            |
 |_________________________________________________________________________|
*/

#define  AGC_COEF       (Shortword)0x19a        /* (1.0 - POST_AGC_COEF)
                                                 * 1.0-.9875 */
#define  POST_EMPHASIS  (Shortword)0x199a       /* 0.2 */

/*_________________________________________________________________________
 |                                                                         |
 |                            Automatic Variables                          |
 |_________________________________________________________________________|
*/

  short int i;

  Longword L_OrigEnergy,
         L_runningGain,
         L_Output;

  Shortword swAgcGain,
         swRunningGain,
         swTemp;

  struct NormSw snsOrigEnergy;

/*_________________________________________________________________________
 |                                                                         |
 |                              Executable Code                            |
 |_________________________________________________________________________|
*/

  /* calculate the energy in the input and save it */
  /*-----------------------------------------------*/

  snsOrigEnergy.sh = g_corr1s(pswSPFIn, swEngyRShift, &L_OrigEnergy);
  snsOrigEnergy.man = round(L_OrigEnergy);

  /* numerator of the postfilter */
  /*-----------------------------*/

  lpcFir(pswSPFIn, pswNumCoef, gpswPostFiltStateNum, pswSPFOut);

  /* denominator of the postfilter */
  /*-------------------------------*/

  lpcIir(pswSPFOut, pswDenomCoef, gpswPostFiltStateDenom, pswSPFOut);

  /* postemphasis section of postfilter */
  /*------------------------------------*/

  for (i = 0; i < S_LEN; i++)
  {
    swTemp = msu_r(L_deposit_h(pswSPFOut[i]), swPostEmphasisState,
                   POST_EMPHASIS);
    swPostEmphasisState = pswSPFOut[i];
    pswSPFOut[i] = swTemp;
  }

  swAgcGain = agcGain(pswSPFOut, snsOrigEnergy, swEngyRShift);

  /* scale the speech vector */
  /*-----------------------------*/

  swRunningGain = gswPostFiltAgcGain;
  L_runningGain = L_deposit_h(gswPostFiltAgcGain);
  for (i = 0; i < S_LEN; i++)
  {
    L_runningGain = L_msu(L_runningGain, swRunningGain, AGC_COEF);
    L_runningGain = L_mac(L_runningGain, swAgcGain, AGC_COEF);
    swRunningGain = extract_h(L_runningGain);

    /* now scale input with gain */

    L_Output = L_mult(swRunningGain, pswSPFOut[i]);
    pswSPFOut[i] = extract_h(L_shl(L_Output, 2));
  }
  gswPostFiltAgcGain = swRunningGain;

}

/***************************************************************************
 *
 *   FUNCTION NAME: speechDecoder
 *
 *   PURPOSE:
 *     The purpose of this function is to call all speech decoder
 *     subroutines.  This is the top level routine for the speech
 *     decoder.
 *
 *   INPUTS:
 *
 *     pswParameters[0:21]
 *
 *        pointer to this frame's parameters.  See below for input
 *        data format.
 *
 *   OUTPUTS:
 *
 *     pswDecodedSpeechFrame[0:159]
 *
 *        this frame's decoded 16 bit linear pcm frame
 *
 *   RETURN VALUE:
 *
 *     none
 *
 *   DESCRIPTION:
 *
 *     The sequence of events in the decoder, and therefore this routine
 *     follow a simple plan.  First, the frame based parameters are
 *     decoded.  Second, on a subframe basis, the subframe based
 *     parameters are decoded and the excitation is generated.  Third,
 *     on a subframe basis, the combined and scaled excitation is
 *     passed through the synthesis filter, and then the pitch and
 *     spectral postfilters.
 *
 *     The in-line comments for the routine speechDecoder, are very
 *     detailed.  Here in a more consolidated form, are the main
 *     points.
 *
 *     The R0 parameter is decoded using the lookup table
 *     psrR0DecTbl[].  The LPC codewords are looked up using lookupVQ().
 *     The decoded parameters are reflection coefficients
 *     (pswFrmKs[]).
 *
 *     The decoder does not use reflection coefficients directly.
 *     Instead it converts them to direct form coeficients.  This is
 *     done using rcToADp().  If this conversion results in invalid
 *     results, the previous frames parameters are used.
 *
 *     The direct form coeficients are used to derive the spectal
 *     postfilter's numerator and denominator coeficients.  The
 *     denominators coefficients are widened, and the numerators
 *     coefficients are a spectrally smoothed version of the
 *     denominator.  The smoothing is done with a_sst().
 *
 *     The frame based LPC coefficients are either used directly as the
 *     subframe coefficients, or are derived through interpolation.
 *     The subframe based coeffiecients are calculated in getSfrmLpc().
 *
 *     Based on voicing mode, the decoder will construct and scale the
 *     excitation in one of two ways.  For the voiced mode the lag is
 *     decoded using lagDecode().  The fractional pitch LTP lookup is
 *     done by the function fp_ex().  In both voiced and unvoiced
 *     mode, the VSELP codewords are decoded into excitation vectors
 *     using b_con() and v_con().
 *
 *     rs_rr(), rs_rrNs(), and scaleExcite() are used to calculate
 *     the gamma's, codevector gains, as well as beta, the LTP vector
 *     gain.  Description of this can be found in section 4.1.11.  Once
 *     the vectors have been scaled and combined, the excitation is
 *     stored in the LTP history.
 *
 *     The excitation, pswExcite[], passes through the pitch pre-filter
 *     (pitchPreFilt()).  Then the harmonically enhanced excitation
 *     passes through the synthesis filter, lpcIir(), and finally the
 *     reconstructed speech passes through the spectral post-filter
 *     (spectalPostFilter()).  The final output speech is passed back in
 *     the array pswDecodedSpeechFrame[].
 *
 *   INPUT DATA FORMAT:
 *
 *      The format/content of the input parameters is the so called
 *      bit alloc format.
 *
 *     voiced mode bit alloc format:
 *     -----------------------------
 *     index    number of bits  parameter name
 *     0        5               R0
 *     1        11              k1Tok3
 *     2        9               k4Tok6
 *     3        8               k7Tok10
 *     4        1               softInterpolation
 *     5        2               voicingDecision
 *     6        8               frameLag
 *     7        9               code_1
 *     8        5               gsp0_1
 *     9        4               deltaLag_2
 *     10       9               code_2
 *     11       5               gsp0_2
 *     12       4               deltaLag_3
 *     13       9               code_3
 *     14       5               gsp0_3
 *     15       4               deltaLag_4
 *     16       9               code_4
 *     17       5               gsp0_4
 *
 *     18       1               BFI
 *     19       1               UFI
 *     20       2               SID
 *     21       1               TAF
 *
 *
 *     unvoiced mode bit alloc format:
 *     -------------------------------
 *
 *     index    number of bits  parameter name
 *     0        5               R0
 *     1        11              k1Tok3
 *     2        9               k4Tok6
 *     3        8               k7Tok10
 *     4        1               softInterpolation
 *     5        2               voicingDecision
 *     6        7               code1_1
 *     7        7               code2_1
 *     8        5               gsp0_1
 *     9        7               code1_2
 *     10       7               code2_2
 *     11       5               gsp0_2
 *     12       7               code1_3
 *     13       7               code2_3
 *     14       5               gsp0_3
 *     15       7               code1_4
 *     16       7               code2_4
 *     17       5               gsp0_4
 *
 *     18       1               BFI
 *     19       1               UFI
 *     20       2               SID
 *     21       1               TAF
 *
 *
 *   REFERENCES: Sub_Clause 4.2 of GSM Recomendation 06.20
 *
 *   KEYWORDS: synthesis, speechdecoder, decoding
 *   KEYWORDS: codewords, lag, codevectors, gsp0
 *
 *************************************************************************/

void   speechDecoder(Shortword pswParameters[],
                            Shortword pswDecodedSpeechFrame[])
{

/*_________________________________________________________________________
 |                                                                         |
 |                           Local Static Variables                        |
 |_________________________________________________________________________|
*/

  static Shortword
        *pswLtpStateOut = &pswLtpStateBaseDec[LTP_LEN],
         pswSythAsSpace[NP * N_SUB],
         pswPFNumAsSpace[NP * N_SUB],
         pswPFDenomAsSpace[NP * N_SUB],
        *ppswSynthAs[N_SUB] = {
    &pswSythAsSpace[0],
    &pswSythAsSpace[10],
    &pswSythAsSpace[20],
    &pswSythAsSpace[30],
  },

        *ppswPFNumAs[N_SUB] = {
    &pswPFNumAsSpace[0],
    &pswPFNumAsSpace[10],
    &pswPFNumAsSpace[20],
    &pswPFNumAsSpace[30],
  },
        *ppswPFDenomAs[N_SUB] = {
    &pswPFDenomAsSpace[0],
    &pswPFDenomAsSpace[10],
    &pswPFDenomAsSpace[20],
    &pswPFDenomAsSpace[30],
  };

  static ShortwordRom
         psrSPFDenomWidenCf[NP] = {
    0x6000, 0x4800, 0x3600, 0x2880, 0x1E60,
    0x16C8, 0x1116, 0x0CD0, 0x099C, 0x0735,
  };


  static Longword L_RxPNSeed;          /* DTX mode */
  static Shortword swRxDtxGsIndex;     /* DTX mode */


/*_________________________________________________________________________
 |                                                                         |
 |                            Automatic Variables                          |
 |_________________________________________________________________________|
*/

  short int i,
         j,
         siLagCode,
         siGsp0Code,
         psiVselpCw[2],
         siVselpCw,
         siNumBits,
         siCodeBook;

  Shortword pswFrmKs[NP],
         pswFrmAs[NP],
         pswFrmPFNum[NP],
         pswFrmPFDenom[NP],
         pswPVec[S_LEN],
         ppswVselpEx[2][S_LEN],
         pswExcite[S_LEN],
         pswPPFExcit[S_LEN],
         pswSynthFiltOut[S_LEN],
         swR0Index,
         swLag,
         swSemiBeta,
         pswBitArray[MAXBITS];

  struct NormSw psnsSqrtRs[N_SUB],
         snsRs00,
         snsRs11,
         snsRs22;


  Shortword swMutePermit,
         swLevelMean,
         swLevelMax,                   /* error concealment */
         swMuteFlag;                   /* error concealment */


  Shortword swTAF,
         swSID,
         swBfiDtx;                     /* DTX mode */
  Shortword swFrameType;               /* DTX mode */

  Longword L_RxDTXGs;                  /* DTX mode */

/*_________________________________________________________________________
 |                                                                         |
 |                              Executable Code                            |
 |_________________________________________________________________________|
*/

  /* -------------------------------------------------------------------- */
  /* do bad frame handling (error concealment) and comfort noise          */
  /* insertion                                                            */
  /* -------------------------------------------------------------------- */


  /* This flag indicates whether muting is performed in the actual frame */
  /* ------------------------------------------------------------------- */
  swMuteFlag = 0;


  /* This flag indicates whether muting is allowed in the actual frame */
  /* ----------------------------------------------------------------- */
  swMutePermit = 0;


  /* frame classification */
  /* -------------------- */

  swSID = pswParameters[20];
  swTAF = pswParameters[21];

  swBfiDtx = pswParameters[18] | pswParameters[19];     /* BFI | UFI */

  if ((swSID == 2) && (swBfiDtx == 0))
    swFrameType = VALIDSID;
  else if ((swSID == 0) && (swBfiDtx == 0))
    swFrameType = GOODSPEECH;
  else if ((swSID == 0) && (swBfiDtx != 0))
    swFrameType = UNUSABLE;
  else
    swFrameType = INVALIDSID;


  /* update of decoder state */
  /* ----------------------- */

  if (swDecoMode == SPEECH)
  {
    /* speech decoding mode */
    /* -------------------- */

    if (swFrameType == VALIDSID)
      swDecoMode = CNIFIRSTSID;
    else if (swFrameType == INVALIDSID)
      swDecoMode = CNIFIRSTSID;
    else if (swFrameType == UNUSABLE)
      swDecoMode = SPEECH;
    else if (swFrameType == GOODSPEECH)
      swDecoMode = SPEECH;
  }
  else
  {
    /* comfort noise insertion mode */
    /* ---------------------------- */

    if (swFrameType == VALIDSID)
      swDecoMode = CNICONT;
    else if (swFrameType == INVALIDSID)
      swDecoMode = CNICONT;
    else if (swFrameType == UNUSABLE)
      swDecoMode = CNIBFI;
    else if (swFrameType == GOODSPEECH)
      swDecoMode = SPEECH;
  }


  if (swDecoMode == SPEECH)
  {
    /* speech decoding mode */
    /* -------------------- */

    /* Perform parameter concealment, depending on BFI (pswParameters[18]) */
    /* or UFI (pswParameters[19])                                          */
    /* ------------------------------------------------------------------- */
    para_conceal_speech_decoder(&pswParameters[18],
                                pswParameters, &swMutePermit);


    /* copy the frame rate parameters */
    /* ------------------------------ */

    swR0Index = pswParameters[0];      /* R0 Index */
    pswVq[0] = pswParameters[1];       /* LPC1 */
    pswVq[1] = pswParameters[2];       /* LPC2 */
    pswVq[2] = pswParameters[3];       /* LPC3 */
    swSi = pswParameters[4];           /* INT_LPC */
    swVoicingMode = pswParameters[5];  /* MODE */


    /* lookup R0 and VQ parameters */
    /* --------------------------- */

    swR0Dec = psrR0DecTbl[swR0Index * 2];       /* R0 */
    lookupVq(pswVq, pswFrmKs);


    /* save this frames GS values */
    /* -------------------------- */

    for (i = 0; i < N_SUB; i++)
    {
      pL_RxGsHist[swRxGsHistPtr] =
              ppLr_gsTable[swVoicingMode][pswParameters[(i * 3) + 8]];
      swRxGsHistPtr++;
      if (swRxGsHistPtr > ((OVERHANG - 1) * N_SUB) - 1)
        swRxGsHistPtr = 0;
    }


    /* DTX variables */
    /* ------------- */

    swDtxBfiCnt = 0;
    swDtxMuting = 0;
    swRxDTXState = CNINTPER - 1;

  }
  else
  {
    /* comfort noise insertion mode */
    /*----------------------------- */

    /* copy the frame rate parameters */
    /* ------------------------------ */

    swR0Index = pswParameters[0];      /* R0 Index */
    pswVq[0] = pswParameters[1];       /* LPC1 */
    pswVq[1] = pswParameters[2];       /* LPC2 */
    pswVq[2] = pswParameters[3];       /* LPC3 */
    swSi = 1;                          /* INT_LPC */
    swVoicingMode = 0;                 /* MODE */


    /* bad frame handling in comfort noise insertion mode */
    /* -------------------------------------------------- */

    if (swDecoMode == CNIFIRSTSID)     /* first SID frame */
    {
      swDtxBfiCnt = 0;
      swDtxMuting = 0;
      swRxDTXState = CNINTPER - 1;

      if (swFrameType == VALIDSID)     /* valid SID frame */
      {
        swR0NewCN = psrR0DecTbl[swR0Index * 2];
        lookupVq(pswVq, pswFrmKs);
      }
      else if (swFrameType == INVALIDSID)       /* invalid SID frame */
      {
        swR0NewCN = psrR0DecTbl[swOldR0IndexDec * 2];
        swR0Index = swOldR0IndexDec;
        for (i = 0; i < NP; i++)
          pswFrmKs[i] = pswOldFrmKsDec[i];
      }

    }
    else if (swDecoMode == CNICONT)    /* SID frame detected, but */
    {                                  /* not the first SID       */
      swDtxBfiCnt = 0;
      swDtxMuting = 0;

      if (swFrameType == VALIDSID)     /* valid SID frame */
      {
        swRxDTXState = 0;
        swR0NewCN = psrR0DecTbl[swR0Index * 2];
        lookupVq(pswVq, pswFrmKs);
      }
      else if (swFrameType == INVALIDSID)       /* invalid SID frame */
      {
        if (swRxDTXState < (CNINTPER - 1))
          swRxDTXState = add(swRxDTXState, 1);
        swR0Index = swOldR0IndexDec;
      }

    }
    else if (swDecoMode == CNIBFI)     /* bad frame received in */
    {                                  /* CNI mode              */
      if (swRxDTXState < (CNINTPER - 1))
        swRxDTXState = add(swRxDTXState, 1);
      swR0Index = swOldR0IndexDec;

      if (swDtxMuting == 1)
      {
        swOldR0IndexDec = sub(swOldR0IndexDec, 2);      /* attenuate
                                                         * by 4 dB */
        if (swOldR0IndexDec < 0)
          swOldR0IndexDec = 0;

        swR0Index = swOldR0IndexDec;

        swR0NewCN = psrR0DecTbl[swOldR0IndexDec * 2];   /* R0 */

      }

      swDtxBfiCnt = add(swDtxBfiCnt, 1);
      if ((swTAF == 1) && (swDtxBfiCnt >= (2 * CNINTPER + 1)))  /* 25 */
        swDtxMuting = 1;

    }


    if (swDecoMode == CNIFIRSTSID)
    {

      /* the first SID frame is received */
      /* ------------------------------- */

      /* initialize the decoders pn-generator */
      /* ------------------------------------ */

      L_RxPNSeed = PN_INIT_SEED;


      /* using the stored rx history, generate averaged GS */
      /* ------------------------------------------------- */

      avgGsHistQntz(pL_RxGsHist, &L_RxDTXGs);
      swRxDtxGsIndex = gsQuant(L_RxDTXGs, 0);

    }


    /* Replace the "transmitted" subframe parameters with */
    /* synthetic ones                                     */
    /* -------------------------------------------------- */

    for (i = 0; i < 4; i++)
    {
      /* initialize the GSP0 parameter */
      pswParameters[(i * 3) + 8] = swRxDtxGsIndex;

      /* CODE1 */
      pswParameters[(i * 3) + 6] = getPnBits(7, &L_RxPNSeed);
      /* CODE2 */
      pswParameters[(i * 3) + 7] = getPnBits(7, &L_RxPNSeed);
    }


    /* Interpolation of CN parameters */
    /* ------------------------------ */

    rxInterpR0Lpc(pswOldFrmKsDec, pswFrmKs, swRxDTXState,
                  swDecoMode, swFrameType);

  }


  /* ------------------- */
  /* do frame processing */
  /* ------------------- */

  /* generate the direct form coefs */
  /* ------------------------------ */

  if (!rcToADp(ASCALE, pswFrmKs, pswFrmAs))
  {

    /* widen direct form coefficients using the widening coefs */
    /* ------------------------------------------------------- */

    for (i = 0; i < NP; i++)
    {
      pswFrmPFDenom[i] = mult_r(pswFrmAs[i], psrSPFDenomWidenCf[i]);
    }

    a_sst(ASHIFT, ASCALE, pswFrmPFDenom, pswFrmPFNum);
  }
  else
  {


    for (i = 0; i < NP; i++)
    {
      pswFrmKs[i] = pswOldFrmKsDec[i];
      pswFrmAs[i] = pswOldFrmAsDec[i];
      pswFrmPFDenom[i] = pswOldFrmPFDenom[i];
      pswFrmPFNum[i] = pswOldFrmPFNum[i];
    }
  }

  /* interpolate, or otherwise get sfrm reflection coefs */
  /* --------------------------------------------------- */

  getSfrmLpc(swSi, swOldR0Dec, swR0Dec, pswOldFrmKsDec, pswOldFrmAsDec,
             pswOldFrmPFNum, pswOldFrmPFDenom, pswFrmKs, pswFrmAs,
             pswFrmPFNum, pswFrmPFDenom, psnsSqrtRs, ppswSynthAs,
             ppswPFNumAs, ppswPFDenomAs);

  /* calculate shift for spectral postfilter */
  /* --------------------------------------- */

  swEngyRShift = r0BasedEnergyShft(swR0Index);


  /* ----------------------- */
  /* do sub-frame processing */
  /* ----------------------- */

  for (giSfrmCnt = 0; giSfrmCnt < 4; giSfrmCnt++)
  {

    /* copy this sub-frame's parameters */
    /* -------------------------------- */

    if (sub(swVoicingMode, 0) == 0)
    {                                  /* unvoiced */
      psiVselpCw[0] = pswParameters[(giSfrmCnt * 3) + 6];       /* CODE_1 */
      psiVselpCw[1] = pswParameters[(giSfrmCnt * 3) + 7];       /* CODE_2 */
      siGsp0Code = pswParameters[(giSfrmCnt * 3) + 8];  /* GSP0 */
    }
    else
    {                                  /* voiced */
      siLagCode = pswParameters[(giSfrmCnt * 3) + 6];   /* LAG */
      psiVselpCw[0] = pswParameters[(giSfrmCnt * 3) + 7];       /* CODE */
      siGsp0Code = pswParameters[(giSfrmCnt * 3) + 8];  /* GSP0 */
    }

    /* for voiced mode, reconstruct the pitch vector */
    /* --------------------------------------------- */

    if (swVoicingMode)
    {

      /* convert delta lag to lag and convert to fractional lag */
      /* ------------------------------------------------------ */

      swLag = lagDecode(siLagCode);

      /* state followed by out */
      /* --------------------- */

      fp_ex(swLag, pswLtpStateOut);

      /* extract a piece of pswLtpStateOut into newly named vector pswPVec */
      /* ----------------------------------------------------------------- */

      for (i = 0; i < S_LEN; i++)
      {
        pswPVec[i] = pswLtpStateOut[i];
      }
    }

    /* for unvoiced, do not reconstruct a pitch vector */
    /* ----------------------------------------------- */

    else
    {
      swLag = 0;                       /* indicates invalid lag
                                        * and unvoiced */
    }

    /* now work on the VSELP codebook excitation output */
    /* x_vec, x_a_vec here named ppswVselpEx[0] and [1] */
    /* ------------------------------------------------ */

    if (swVoicingMode)
    {                                  /* voiced */

      siNumBits = C_BITS_V;
      siVselpCw = psiVselpCw[0];

      b_con(siVselpCw, siNumBits, pswBitArray);

      v_con(pppsrVcdCodeVec[0][0], ppswVselpEx[0], pswBitArray, siNumBits);
    }

    else
    {                                  /* unvoiced */

      siNumBits = C_BITS_UV;

      for (siCodeBook = 0; siCodeBook < 2; siCodeBook++)
      {

        siVselpCw = psiVselpCw[siCodeBook];

        b_con(siVselpCw, siNumBits, (Shortword *) pswBitArray);

        v_con(pppsrUvCodeVec[siCodeBook][0], ppswVselpEx[siCodeBook],
              pswBitArray, siNumBits);
      }
    }

    /* all excitation vectors have been created: ppswVselpEx and pswPVec  */
    /* if voiced compute rs00 and rs11; if unvoiced cmpute rs11 and rs22  */
    /* ------------------------------------------------------------------ */

    if (swLag)
    {
      rs_rr(pswPVec, psnsSqrtRs[giSfrmCnt], &snsRs00);
    }

    rs_rrNs(ppswVselpEx[0], psnsSqrtRs[giSfrmCnt], &snsRs11);

    if (!swVoicingMode)
    {
      rs_rrNs(ppswVselpEx[1], psnsSqrtRs[giSfrmCnt], &snsRs22);
    }

    /* now implement synthesis - combine the excitations */
    /* ------------------------------------------------- */

    if (swVoicingMode)
    {                                  /* voiced */

      /* scale pitch and codebook excitations and get beta */
      /* ------------------------------------------------- */
      swSemiBeta = scaleExcite(pswPVec,
                               pppsrGsp0[swVoicingMode][siGsp0Code][0],
                               snsRs00, pswPVec);
      scaleExcite(ppswVselpEx[0],
                  pppsrGsp0[swVoicingMode][siGsp0Code][1],
                  snsRs11, ppswVselpEx[0]);

      /* combine the two scaled excitations */
      /* ---------------------------------- */
      for (i = 0; i < S_LEN; i++)
      {
        pswExcite[i] = add(pswPVec[i], ppswVselpEx[0][i]);
      }
    }
    else
    {                                  /* unvoiced */

      /* scale codebook excitations and set beta to 0 as not applicable */
      /* -------------------------------------------------------------- */
      swSemiBeta = 0;
      scaleExcite(ppswVselpEx[0],
                  pppsrGsp0[swVoicingMode][siGsp0Code][0],
                  snsRs11, ppswVselpEx[0]);
      scaleExcite(ppswVselpEx[1],
                  pppsrGsp0[swVoicingMode][siGsp0Code][1],
                  snsRs22, ppswVselpEx[1]);

      /* combine the two scaled excitations */
      /* ---------------------------------- */
      for (i = 0; i < S_LEN; i++)
      {
        pswExcite[i] = add(ppswVselpEx[1][i], ppswVselpEx[0][i]);
      }
    }

    /* now update the pitch state using the combined/scaled excitation */
    /* --------------------------------------------------------------- */

    for (i = 0; i < LTP_LEN; i++)
    {
      pswLtpStateBaseDec[i] = pswLtpStateBaseDec[i + S_LEN];
    }

    /* add the current sub-frames data to the state */
    /* -------------------------------------------- */

    for (i = -S_LEN, j = 0; j < S_LEN; i++, j++)
    {
      pswLtpStateOut[i] = pswExcite[j];/* add new frame at t = -S_LEN */
    }

    /* given the excitation perform pitch prefiltering */
    /* ----------------------------------------------- */

    pitchPreFilt(pswExcite, siGsp0Code, swLag,
                 swVoicingMode, swSemiBeta, psnsSqrtRs[giSfrmCnt],
                 pswPPFExcit, pswPPreState);


    /* Concealment on subframe signal level: */
    /* ------------------------------------- */
    level_estimator(0, &swLevelMean, &swLevelMax,
                    &pswDecodedSpeechFrame[giSfrmCnt * S_LEN]);

    signal_conceal_sub(pswPPFExcit, ppswSynthAs[giSfrmCnt], pswSynthFiltState,
                    &pswLtpStateOut[-S_LEN], &pswPPreState[LTP_LEN - S_LEN],
                       swLevelMean, swLevelMax,
                       pswParameters[19], swMuteFlagOld,
                       &swMuteFlag, swMutePermit);


    /* synthesize the speech through the synthesis filter */
    /* -------------------------------------------------- */

    lpcIir(pswPPFExcit, ppswSynthAs[giSfrmCnt], pswSynthFiltState,
           pswSynthFiltOut);

    /* pass reconstructed speech through adaptive spectral postfilter */
    /* -------------------------------------------------------------- */

    spectralPostFilter(pswSynthFiltOut, ppswPFNumAs[giSfrmCnt],
                       ppswPFDenomAs[giSfrmCnt],
                       &pswDecodedSpeechFrame[giSfrmCnt * S_LEN]);

    level_estimator(1, &swLevelMean, &swLevelMax,
                    &pswDecodedSpeechFrame[giSfrmCnt * S_LEN]);

  }

  /* Save muting information for next frame */
  /* -------------------------------------- */
  swMuteFlagOld = swMuteFlag;

  /* end of frame processing - save this frame's frame energy,  */
  /* reflection coefs, direct form coefs, and post filter coefs */
  /* ---------------------------------------------------------- */

  swOldR0Dec = swR0Dec;
  swOldR0IndexDec = swR0Index;         /* DTX mode */

  for (i = 0; i < NP; i++)
  {
    pswOldFrmKsDec[i] = pswFrmKs[i];
    pswOldFrmAsDec[i] = pswFrmAs[i];
    pswOldFrmPFNum[i] = pswFrmPFNum[i];
    pswOldFrmPFDenom[i] = pswFrmPFDenom[i];
  }
}


/***************************************************************************
 *
 *   FUNCTION NAME: sqroot
 *
 *   PURPOSE:
 *
 *     The purpose of this function is to perform a single precision square
 *     root function on a Longword
 *
 *   INPUTS:
 *
 *     L_SqrtIn
 *                     input to square root function
 *
 *   OUTPUTS:
 *
 *     none
 *
 *   RETURN VALUE:
 *
 *     swSqrtOut
 *                     output to square root function
 *
 *   DESCRIPTION:
 *
 *      Input assumed to be normalized
 *
 *      The algorithm is based around a six term Taylor expansion :
 *
 *        y^0.5 = (1+x)^0.5
 *             ~= 1 + (x/2) - 0.5*((x/2)^2) + 0.5*((x/2)^3)
 *                - 0.625*((x/2)^4) + 0.875*((x/2)^5)
 *
 *      Max error less than 0.08 % for normalized input ( 0.5 <= x < 1 )
 *
 *   REFERENCES: Sub-clause 4.1.4.1, 4.1.7, 4.1.11.1, 4.2.1,
 *               4.2.2, 4.2.3 and 4.2.4 of GSM Recomendation 06.20
 *
 *   KEYWORDS: sqrt, squareroot, sqrt016
 *
 *************************************************************************/

Shortword sqroot(Longword L_SqrtIn)
{

/*_________________________________________________________________________
 |                                                                         |
 |                              Local Constants                            |
 |_________________________________________________________________________|
*/

#define    PLUS_HALF          0x40000000L       /* 0.5 */
#define    MINUS_ONE          0x80000000L       /* -1 */
#define    TERM5_MULTIPLER    0x5000   /* 0.625 */
#define    TERM6_MULTIPLER    0x7000   /* 0.875 */

/*_________________________________________________________________________
 |                                                                         |
 |                            Automatic Variables                          |
 |_________________________________________________________________________|
*/

  Longword L_Temp0,
         L_Temp1;

  Shortword swTemp,
         swTemp2,
         swTemp3,
         swTemp4,
         swSqrtOut;

/*_________________________________________________________________________
 |                                                                         |
 |                              Executable Code                            |
 |_________________________________________________________________________|
*/

  /* determine 2nd term x/2 = (y-1)/2 */
  /* -------------------------------- */

  L_Temp1 = L_shr(L_SqrtIn, 1);        /* L_Temp1 = y/2 */
  L_Temp1 = L_sub(L_Temp1, PLUS_HALF); /* L_Temp1 = (y-1)/2 */
  swTemp = extract_h(L_Temp1);         /* swTemp = x/2 */

  /* add contribution of 2nd term */
  /* ---------------------------- */

  L_Temp1 = L_sub(L_Temp1, MINUS_ONE); /* L_Temp1 = 1 + x/2 */

  /* determine 3rd term */
  /* ------------------ */

  L_Temp0 = L_msu(0L, swTemp, swTemp); /* L_Temp0 = -(x/2)^2 */
  swTemp2 = extract_h(L_Temp0);        /* swTemp2 = -(x/2)^2 */
  L_Temp0 = L_shr(L_Temp0, 1);         /* L_Temp0 = -0.5*(x/2)^2 */

  /* add contribution of 3rd term */
  /* ---------------------------- */

  L_Temp0 = L_add(L_Temp1, L_Temp0);   /* L_Temp0 = 1 + x/2 - 0.5*(x/2)^2 */

  /* determine 4rd term */
  /* ------------------ */

  L_Temp1 = L_msu(0L, swTemp, swTemp2);/* L_Temp1 = (x/2)^3 */
  swTemp3 = extract_h(L_Temp1);        /* swTemp3 = (x/2)^3 */
  L_Temp1 = L_shr(L_Temp1, 1);         /* L_Temp1 = 0.5*(x/2)^3 */

  /* add contribution of 4rd term */
  /* ---------------------------- */

  /* L_Temp1 = 1 + x/2 - 0.5*(x/2)^2 + 0.5*(x/2)^3 */

  L_Temp1 = L_add(L_Temp0, L_Temp1);

  /* determine partial 5th term */
  /* -------------------------- */

  L_Temp0 = L_mult(swTemp, swTemp3);   /* L_Temp0 = (x/2)^4 */
  swTemp4 = round(L_Temp0);            /* swTemp4 = (x/2)^4 */

  /* determine partial 6th term */
  /* -------------------------- */

  L_Temp0 = L_msu(0L, swTemp2, swTemp3);        /* L_Temp0 = (x/2)^5 */
  swTemp2 = round(L_Temp0);            /* swTemp2 = (x/2)^5 */

  /* determine 5th term and add its contribution */
  /* ------------------------------------------- */

  /* L_Temp0 = -0.625*(x/2)^4 */

  L_Temp0 = L_msu(0L, TERM5_MULTIPLER, swTemp4);

  /* L_Temp1 = 1 + x/2 - 0.5*(x/2)^2 + 0.5*(x/2)^3 - 0.625*(x/2)^4 */

  L_Temp1 = L_add(L_Temp0, L_Temp1);

  /* determine 6th term and add its contribution */
  /* ------------------------------------------- */

  /* swSqrtOut = 1 + x/2 - 0.5*(x/2)^2 + 0.5*(x/2)^3 */
  /* - 0.625*(x/2)^4 + 0.875*(x/2)^5     */

  swSqrtOut = mac_r(L_Temp1, TERM6_MULTIPLER, swTemp2);

  /* return output */
  /* ------------- */

  return (swSqrtOut);
}

/***************************************************************************
 *
 *   FUNCTION NAME: v_con
 *
 *   PURPOSE:
 *
 *     This subroutine constructs a codebook excitation
 *     vector from basis vectors
 *
 *   INPUTS:
 *
 *     pswBVects[0:siNumBVctrs*S_LEN-1]
 *
 *                     Array containing a set of basis vectors.
 *
 *     pswBitArray[0:siNumBVctrs-1]
 *
 *                     Bit array dictating the polarity of the
 *                     basis vectors in the output vector.
 *                     Each element of the bit array is either -0.5 or +0.5
 *
 *     siNumBVctrs
 *                     Number of bits in codeword
 *
 *   OUTPUTS:
 *
 *     pswOutVect[0:39]
 *
 *                     Array containing the contructed output vector
 *
 *   RETURN VALUE:
 *
 *     none
 *
 *   DESCRIPTION:
 *
 *     The array pswBitArray is used to multiply each of the siNumBVctrs
 *     basis vectors.  The input pswBitArray[] is an array whose
 *     elements are +/-0.5.  These multiply the VSELP basis vectors and
 *     when summed produce a VSELP codevector.  b_con() is the function
 *     used to translate a VSELP codeword into pswBitArray[].
 *
 *
 *   REFERENCES: Sub-clause 4.1.10 and 4.2.1 of GSM Recomendation 06.20
 *
 *   KEYWORDS: v_con, codeword, reconstruct, basis vector, excitation
 *
 *************************************************************************/

void   v_con(Shortword pswBVects[], Shortword pswOutVect[],
                    Shortword pswBitArray[], short int siNumBVctrs)
{

/*_________________________________________________________________________
 |                                                                         |
 |                            Automatic Variables                          |
 |_________________________________________________________________________|
*/

  Longword L_Temp;

  short int siSampleCnt,
         siCVectCnt;

/*_________________________________________________________________________
 |                                                                         |
 |                              Executable Code                            |
 |_________________________________________________________________________|
*/

  /* Sample loop  */
  /*--------------*/
  for (siSampleCnt = 0; siSampleCnt < S_LEN; siSampleCnt++)
  {

    /* First element of output vector */
    L_Temp = L_mult(pswBitArray[0], pswBVects[0 * S_LEN + siSampleCnt]);

    /* Construct output vector */
    /*-------------------------*/
    for (siCVectCnt = 1; siCVectCnt < siNumBVctrs; siCVectCnt++)
    {
      L_Temp = L_mac(L_Temp, pswBitArray[siCVectCnt],
                     pswBVects[siCVectCnt * S_LEN + siSampleCnt]);
    }

    /* store the output vector sample */
    /*--------------------------------*/
    L_Temp = L_shl(L_Temp, 1);
    pswOutVect[siSampleCnt] = extract_h(L_Temp);
  }
}
