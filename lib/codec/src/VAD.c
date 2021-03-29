/****************************************************************************
 *
 *     TITLE:     Half-Rate GSM Voice Activity Detector (VAD) Modules
 *
 *     VERSION:   1.2
 *
 *     REFERENCE: Recommendation GSM 06.42
 *
 ***************************************************************************/

/*_________________________________________________________________________
 |                                                                         |
 |                              Include Files                              |
 |_________________________________________________________________________|
*/

#include "typedefs.h"
#include "mathhalf.h"
#include "mathdp31.h"
#include "vad.h"


/*_________________________________________________________________________
 |                                                                         |
 |                              Local Defines                              |
 |_________________________________________________________________________|
*/

/*** Floating point representations of constants pth, plev and margin ***/

#define M_PTH    26250
#define E_PTH    18
#define M_PLEV   17500
#define E_PLEV   20
#define M_MARGIN 27343
#define E_MARGIN 27

/*_________________________________________________________________________
 |                                                                         |
 |                            Static Variables                             |
 |_________________________________________________________________________|
*/

static Shortword
       pswRvad[9],
       swNormRvad,
       swPt_sacf,
       swPt_sav0,
       swE_thvad,
       swM_thvad,
       swAdaptCount,
       swBurstCount,
       swHangCount,
       swOldLagCount,
       swVeryOldLagCount,
       swOldLag;

static Longword
       pL_sacf[27],
       pL_sav0[36],
       L_lastdm;

/****************************************************************************
 *
 *     FUNCTION:  vad_reset
 *
 *     VERSION:   1.2
 *
 *     PURPOSE:   Resets VAD static variables to their initial value.
 *
 ***************************************************************************/

void   vad_reset(void)

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

  pswRvad[0] = 24576;
  swNormRvad = 7;
  swPt_sacf = 0;
  swPt_sav0 = 0;
  L_lastdm = 0;
  swE_thvad = 21;
  swM_thvad = 21875;
  swAdaptCount = 0;
  swBurstCount = 0;
  swHangCount = -1;
  swOldLagCount = 0;
  swVeryOldLagCount = 0;
  swOldLag = 21;

  for (i = 1; i < 9; i++)
    pswRvad[i] = 0;
  for (i = 0; i < 27; i++)
    pL_sacf[i] = 0;
  for (i = 0; i < 36; i++)
    pL_sav0[i] = 0;

}

/****************************************************************************
 *
 *     FUNCTION:  vad_algorithm
 *
 *     VERSION:   1.2
 *
 *     PURPOSE:   Returns a decision as to whether the current frame being
 *                processed by the speech encoder contains speech or not.
 *
 *     INPUTS:    pL_acf[0..8]  autocorrelation of input signal frame
 *                swScaleAcf    L_acf scaling factor
 *                pswRc[0..3]   speech encoder reflection coefficients
 *                swPtch        flag to indicate a periodic signal component
 *
 *     OUTPUTS:   pswVadFlag    vad decision
 *
 ***************************************************************************/

void   vad_algorithm(Longword pL_acf[9],
                            Shortword swScaleAcf,
                            Shortword pswRc[4],
                            Shortword swPtch,
                            Shortword *pswVadFlag)
{

/*_________________________________________________________________________
 |                                                                         |
 |                            Automatic Variables                          |
 |_________________________________________________________________________|
*/

  Longword
         pL_av0[9],
         pL_av1[9];

  Shortword
         swM_acf0,
         swE_acf0,
         pswRav1[9],
         swNormRav1,
         swM_pvad,
         swE_pvad,
         swStat,
         swTone,
         swVvad;


/*_________________________________________________________________________
 |                                                                         |
 |                              Executable Code                            |
 |_________________________________________________________________________|
*/

  energy_computation
          (
           pL_acf, swScaleAcf,
           pswRvad, swNormRvad,
           &swM_pvad, &swE_pvad,
           &swM_acf0, &swE_acf0
          );

  average_acf
          (
           pL_acf, swScaleAcf,
           pL_av0, pL_av1
          );

  predictor_values
          (
           pL_av1,
           pswRav1,
           &swNormRav1
          );

  spectral_comparison
          (
           pswRav1, swNormRav1,
           pL_av0,
           &swStat
          );

  tone_detection
          (
           pswRc,
           &swTone
          );

  threshold_adaptation
          (
           swStat, swPtch, swTone,
           pswRav1, swNormRav1,
           swM_pvad, swE_pvad,
           swM_acf0, swE_acf0,
           pswRvad, &swNormRvad,
           &swM_thvad, &swE_thvad
          );

  vad_decision
          (
           swM_pvad, swE_pvad,
           swM_thvad, swE_thvad,
           &swVvad
          );

  vad_hangover
          (
           swVvad,
           pswVadFlag
          );

}

/****************************************************************************
 *
 *     FUNCTION:  energy_computation
 *
 *     VERSION:   1.2
 *
 *     PURPOSE:   Computes the input and residual energies of the adaptive
 *                filter in a floating point representation.
 *
 *     INPUTS:    pL_acf[0..8]   autocorrelation of input signal frame
 *                swScaleAcf     L_acf scaling factor
 *                pswRvad[0..8]  autocorrelated adaptive filter coefficients
 *                swNormRvad     rvad scaling factor
 *
 *     OUTPUTS:   pswM_pvad      mantissa of filtered signal energy
 *                pswE_pvad      exponent of filtered signal energy
 *                pswM_acf0      mantissa of signal frame energy
 *                pswE_acf0      exponent of signal frame energy
 *
 ***************************************************************************/

void   energy_computation(Longword pL_acf[],
                                 Shortword swScaleAcf,
                                 Shortword pswRvad[],
                                 Shortword swNormRvad,
                                 Shortword *pswM_pvad,
                                 Shortword *pswE_pvad,
                                 Shortword *pswM_acf0,
                                 Shortword *pswE_acf0)
{

/*_________________________________________________________________________
 |                                                                         |
 |                            Automatic Variables                          |
 |_________________________________________________________________________|
*/

  Longword
         L_temp;

  Shortword
         pswSacf[9],
         swNormAcf,
         swNormProd,
         swShift;

  int
         i;


/*_________________________________________________________________________
 |                                                                         |
 |                              Executable Code                            |
 |_________________________________________________________________________|
*/

  /*** Test if acf[0] is zero ***/

  if (pL_acf[0] == 0)
  {
    *pswE_pvad = -0x8000;
    *pswM_pvad = 0;
    *pswE_acf0 = -0x8000;
    *pswM_acf0 = 0;
    return;
  }


  /*** Re-normalisation of L_acf[0..8] ***/

  swNormAcf = norm_l(pL_acf[0]);
  swShift = sub(swNormAcf, 3);

  for (i = 0; i <= 8; i++)
    pswSacf[i] = extract_h(L_shl(pL_acf[i], swShift));


  /*** Computation of e_acf0 and m_acf0 ***/

  *pswE_acf0 = add(32, shl(swScaleAcf, 1));
  *pswE_acf0 = sub(*pswE_acf0, swNormAcf);
  *pswM_acf0 = shl(pswSacf[0], 3);


  /*** Computation of e_pvad and m_pvad ***/

  *pswE_pvad = add(*pswE_acf0, 14);
  *pswE_pvad = sub(*pswE_pvad, swNormRvad);

  L_temp = 0;

  for (i = 1; i <= 8; i++)
    L_temp = L_mac(L_temp, pswSacf[i], pswRvad[i]);

  L_temp = L_add(L_temp, L_shr(L_mult(pswSacf[0], pswRvad[0]), 1));

  if (L_temp <= 0)
    L_temp = 1;

  swNormProd = norm_l(L_temp);
  *pswE_pvad = sub(*pswE_pvad, swNormProd);
  *pswM_pvad = extract_h(L_shl(L_temp, swNormProd));

}

/****************************************************************************
 *
 *     FUNCTION:  average_acf
 *
 *     VERSION:   1.2
 *
 *     PURPOSE:   Computes the arrays L_av0 [0..8] and L_av1 [0..8].
 *
 *     INPUTS:    pL_acf[0..8]  autocorrelation of input signal frame
 *                swScaleAcf    L_acf scaling factor
 *
 *     OUTPUTS:   pL_av0[0..8]  ACF averaged over last four frames
 *                pL_av1[0..8]  ACF averaged over previous four frames
 *
 ***************************************************************************/

void   average_acf(Longword pL_acf[],
                          Shortword swScaleAcf,
                          Longword pL_av0[],
                          Longword pL_av1[])
{

/*_________________________________________________________________________
 |                                                                         |
 |                            Automatic Variables                          |
 |_________________________________________________________________________|
*/

  Longword L_temp;

  Shortword swScale;

  int    i;

/*_________________________________________________________________________
 |                                                                         |
 |                              Executable Code                            |
 |_________________________________________________________________________|
*/

  /*** computation of the scaleing factor ***/

  swScale = sub(10, shl(swScaleAcf, 1));


  /*** Computation of the arrays L_av0 and L_av1 ***/

  for (i = 0; i <= 8; i++)
  {
    L_temp = L_shr(pL_acf[i], swScale);
    pL_av0[i] = L_add(pL_sacf[i], L_temp);
    pL_av0[i] = L_add(pL_sacf[i + 9], pL_av0[i]);
    pL_av0[i] = L_add(pL_sacf[i + 18], pL_av0[i]);
    pL_sacf[swPt_sacf + i] = L_temp;
    pL_av1[i] = pL_sav0[swPt_sav0 + i];
    pL_sav0[swPt_sav0 + i] = pL_av0[i];
  }


  /*** Update the array pointers ***/

  if (swPt_sacf == 18)
    swPt_sacf = 0;
  else
    swPt_sacf = add(swPt_sacf, 9);

  if (swPt_sav0 == 27)
    swPt_sav0 = 0;
  else
    swPt_sav0 = add(swPt_sav0, 9);

}

/****************************************************************************
 *
 *     FUNCTION:  predictor_values
 *
 *     VERSION:   1.2
 *
 *     PURPOSE:   Computes the array rav [0..8] needed for the spectral
 *                comparison and the threshold adaptation.
 *
 *     INPUTS:    pL_av1 [0..8]  ACF averaged over previous four frames
 *
 *     OUTPUTS:   pswRav1 [0..8] ACF obtained from L_av1
 *                pswNormRav1    r_av1 scaling factor
 *
 ***************************************************************************/

void   predictor_values(Longword pL_av1[],
                               Shortword pswRav1[],
                               Shortword *pswNormRav1)
{

/*_________________________________________________________________________
 |                                                                         |
 |                            Automatic Variables                          |
 |_________________________________________________________________________|
*/

  Shortword
         pswVpar[8],
         pswAav1[9];

/*_________________________________________________________________________
 |                                                                         |
 |                              Executable Code                            |
 |_________________________________________________________________________|
*/

  schur_recursion(pL_av1, pswVpar);
  step_up(8, pswVpar, pswAav1);
  compute_rav1(pswAav1, pswRav1, pswNormRav1);

}

/****************************************************************************
 *
 *     FUNCTION:  schur_recursion
 *
 *     VERSION:   1.2
 *
 *     PURPOSE:   Uses the Schur recursion to compute adaptive filter
 *                reflection coefficients from an autorrelation function.
 *
 *     INPUTS:    pL_av1[0..8]   autocorrelation function
 *
 *     OUTPUTS:   pswVpar[0..7]  reflection coefficients
 *
 ***************************************************************************/

void   schur_recursion(Longword pL_av1[],
                              Shortword pswVpar[])
{

/*_________________________________________________________________________
 |                                                                         |
 |                            Automatic Variables                          |
 |_________________________________________________________________________|
*/

  Shortword
         pswAcf[9],
         pswPp[9],
         pswKk[9],
         swTemp;

  int    i,
         k,
         m,
         n;


/*_________________________________________________________________________
 |                                                                         |
 |                              Executable Code                            |
 |_________________________________________________________________________|
*/

  /*** Schur recursion with 16-bit arithmetic ***/

  if (pL_av1[0] == 0)
  {
    for (i = 0; i < 8; i++)
      pswVpar[i] = 0;
    return;
  }

  swTemp = norm_l(pL_av1[0]);

  for (k = 0; k <= 8; k++)
    pswAcf[k] = extract_h(L_shl(pL_av1[k], swTemp));


  /*** Initialise array pp[..] and kk[..] for the recursion: ***/

  for (i = 1; i <= 7; i++)
    pswKk[9 - i] = pswAcf[i];

  for (i = 0; i <= 8; i++)
    pswPp[i] = pswAcf[i];


  /*** Compute Parcor coefficients: ***/

  for (n = 0; n < 8; n++)
  {
    if (pswPp[0] < abs_s(pswPp[1]))
    {
      for (i = n; i < 8; i++)
        pswVpar[i] = 0;
      return;
    }
    pswVpar[n] = divide_s(abs_s(pswPp[1]), pswPp[0]);

    if (pswPp[1] > 0)
      pswVpar[n] = negate(pswVpar[n]);
    if (n == 7)
      return;


    /*** Schur recursion: ***/

    pswPp[0] = add(pswPp[0], mult_r(pswPp[1], pswVpar[n]));

    for (m = 1; m <= (7 - n); m++)
    {
      pswPp[m] = add(pswPp[1 + m], mult_r(pswKk[9 - m], pswVpar[n]));
      pswKk[9 - m] = add(pswKk[9 - m], mult_r(pswPp[1 + m], pswVpar[n]));
    }
  }

}

/****************************************************************************
 *
 *     FUNCTION:  step_up
 *
 *     VERSION:   1.2
 *
 *     PURPOSE:   Computes the transversal filter coefficients from the
 *                reflection coefficients.
 *
 *     INPUTS:    swNp             filter order (2..8)
 *                pswVpar[0..np-1] reflection coefficients
 *
 *     OUTPUTS:   pswAav1[0..np]   transversal filter coefficients
 *
 ***************************************************************************/

void   step_up(Shortword swNp,
                      Shortword pswVpar[],
                      Shortword pswAav1[])
{

/*_________________________________________________________________________
 |                                                                         |
 |                            Automatic Variables                          |
 |_________________________________________________________________________|
*/

  Longword
         pL_coef[9],
         pL_work[9];

  Shortword
         swTemp;

  int
         i,
         m;


/*_________________________________________________________________________
 |                                                                         |
 |                              Executable Code                            |
 |_________________________________________________________________________|
*/

  /*** Initialisation of the step-up recursion ***/

  pL_coef[0] = L_shl(0x4000L, 15);
  pL_coef[1] = L_shl(L_deposit_l(pswVpar[0]), 14);


  /*** Loop on the LPC analysis order: ***/

  for (m = 2; m <= swNp; m++)
  {
    for (i = 1; i < m; i++)
    {
      swTemp = extract_h(pL_coef[m - i]);
      pL_work[i] = L_mac(pL_coef[i], pswVpar[m - 1], swTemp);
    }
    for (i = 1; i < m; i++)
      pL_coef[i] = pL_work[i];
    pL_coef[m] = L_shl(L_deposit_l(pswVpar[m - 1]), 14);
  }


  /*** Keep the aav1[0..swNp] in 15 bits for the following subclause ***/

  for (i = 0; i <= swNp; i++)
    pswAav1[i] = extract_h(L_shr(pL_coef[i], 3));

}

/****************************************************************************
 *
 *     FUNCTION:  compute_rav1
 *
 *     VERSION:   1.2
 *
 *     PURPOSE:   Computes the autocorrelation function of the adaptive
 *                filter coefficients.
 *
 *     INPUTS:    pswAav1[0..8]  adaptive filter coefficients
 *
 *     OUTPUTS:   pswRav1[0..8]  ACF of aav1
 *                pswNormRav1    r_av1 scaling factor
 *
 ***************************************************************************/

void   compute_rav1(Shortword pswAav1[],
                           Shortword pswRav1[],
                           Shortword *pswNormRav1)
{

/*_________________________________________________________________________
 |                                                                         |
 |                            Automatic Variables                          |
 |_________________________________________________________________________|
*/

  Longword
         pL_work[9];

  int
         i,
         k;


/*_________________________________________________________________________
 |                                                                         |
 |                              Executable Code                            |
 |_________________________________________________________________________|
*/

  /*** Computation of the rav1[0..8] ***/

  for (i = 0; i <= 8; i++)
  {
    pL_work[i] = 0;

    for (k = 0; k <= 8 - i; k++)
      pL_work[i] = L_mac(pL_work[i], pswAav1[k], pswAav1[k + i]);
  }

  if (pL_work[0] == 0)
    *pswNormRav1 = 0;
  else
    *pswNormRav1 = norm_l(pL_work[0]);

  for (i = 0; i <= 8; i++)
    pswRav1[i] = extract_h(L_shl(pL_work[i], *pswNormRav1));

}

/****************************************************************************
 *
 *     FUNCTION:  spectral_comparison
 *
 *     VERSION:   1.2
 *
 *     PURPOSE:   Computes the stat flag needed for the threshold
 *                adaptation decision.
 *
 *     INPUTS:    pswRav1[0..8]   ACF obtained from L_av1
 *                swNormRav1      rav1 scaling factor
 *                pL_av0[0..8]    ACF averaged over last four frames
 *
 *     OUTPUTS:   pswStat         flag to indicate spectral stationarity
 *
 ***************************************************************************/

void   spectral_comparison(Shortword pswRav1[],
                                  Shortword swNormRav1,
                                  Longword pL_av0[],
                                  Shortword *pswStat)
{

/*_________________________________________________________________________
 |                                                                         |
 |                            Automatic Variables                          |
 |_________________________________________________________________________|
*/

  Longword
         L_dm,
         L_sump,
         L_temp;

  Shortword
         pswSav0[9],
         swShift,
         swDivShift,
         swTemp;

  int
         i;


/*_________________________________________________________________________
 |                                                                         |
 |                              Executable Code                            |
 |_________________________________________________________________________|
*/

  /*** Re-normalise L_av0[0..8] ***/

  if (pL_av0[0] == 0)
  {
    for (i = 0; i <= 8; i++)
      pswSav0[i] = 4095;
  }
  else
  {
    swShift = sub(norm_l(pL_av0[0]), 3);
    for (i = 0; i <= 8; i++)
      pswSav0[i] = extract_h(L_shl(pL_av0[i], swShift));
  }

  /*** compute partial sum of dm ***/

  L_sump = 0;

  for (i = 1; i <= 8; i++)
    L_sump = L_mac(L_sump, pswRav1[i], pswSav0[i]);

  /*** compute the division of the partial sum by sav0[0] ***/

  if (L_sump < 0)
    L_temp = L_negate(L_sump);
  else
    L_temp = L_sump;

  if (L_temp == 0)
  {
    L_dm = 0;
    swShift = 0;
  }
  else
  {
    pswSav0[0] = shl(pswSav0[0], 3);
    swShift = norm_l(L_temp);
    swTemp = extract_h(L_shl(L_temp, swShift));

    if (pswSav0[0] >= swTemp)
    {
      swDivShift = 0;
      swTemp = divide_s(swTemp, pswSav0[0]);
    }
    else
    {
      swDivShift = 1;
      swTemp = sub(swTemp, pswSav0[0]);
      swTemp = divide_s(swTemp, pswSav0[0]);
    }

    if (swDivShift == 1)
      L_dm = 0x8000L;
    else
      L_dm = 0;

    L_dm = L_shl((L_add(L_dm, L_deposit_l(swTemp))), 1);

    if (L_sump < 0)
      L_dm = L_sub(0L, L_dm);
  }


  /*** Re-normalisation and final computation of L_dm ***/

  L_dm = L_shl(L_dm, 14);
  L_dm = L_shr(L_dm, swShift);
  L_dm = L_add(L_dm, L_shl(L_deposit_l(pswRav1[0]), 11));
  L_dm = L_shr(L_dm, swNormRav1);


  /*** Compute the difference and save L_dm ***/

  L_temp = L_sub(L_dm, L_lastdm);
  L_lastdm = L_dm;

  if (L_temp < 0L)
    L_temp = L_negate(L_temp);


  /*** Evaluation of the state flag  ***/

  L_temp = L_sub(L_temp, 4456L);

  if (L_temp < 0)
    *pswStat = 1;
  else
    *pswStat = 0;

}

/****************************************************************************
 *
 *     FUNCTION:  tone_detection
 *
 *     VERSION:   1.2
 *
 *     PURPOSE:   Computes the tone flag needed for the threshold
 *                adaptation decision.
 *
 *     INPUTS:    pswRc[0..3] reflection coefficients calculated in the
 *                            speech encoder short term predictor
 *
 *     OUTPUTS:   pswTone     flag to indicate a periodic signal component
 *
 ***************************************************************************/

void   tone_detection(Shortword pswRc[4],
                             Shortword *pswTone)
{

/*_________________________________________________________________________
 |                                                                         |
 |                            Automatic Variables                          |
 |_________________________________________________________________________|
*/

  Longword
         L_num,
         L_den,
         L_temp;

  Shortword
         swTemp,
         swPredErr,
         pswA[3];

  int
         i;

/*_________________________________________________________________________
 |                                                                         |
 |                              Executable Code                            |
 |_________________________________________________________________________|
*/

  *pswTone = 0;


  /*** Calculate filter coefficients ***/

  step_up(2, pswRc, pswA);


  /*** Calculate ( a[1] * a[1] ) ***/

  swTemp = shl(pswA[1], 3);
  L_den = L_mult(swTemp, swTemp);


  /*** Calculate ( 4*a[2] - a[1]*a[1] ) ***/

  L_temp = L_shl(L_deposit_h(pswA[2]), 3);
  L_num = L_sub(L_temp, L_den);


  /*** Check if pole frequency is less than 385 Hz ***/

  if (L_num <= 0)
    return;

  if (pswA[1] < 0)
  {
    swTemp = extract_h(L_den);
    L_temp = L_msu(L_num, swTemp, 3189);

    if (L_temp < 0)
      return;
  }


  /*** Calculate normalised prediction error ***/

  swPredErr = 0x7fff;

  for (i = 0; i < 4; i++)
  {
    swTemp = mult(pswRc[i], pswRc[i]);
    swTemp = sub(32767, swTemp);
    swPredErr = mult(swPredErr, swTemp);
  }


  /*** Test if prediction error is smaller than threshold ***/

  swTemp = sub(swPredErr, 1464);

  if (swTemp < 0)
    *pswTone = 1;

}

/****************************************************************************
 *
 *     FUNCTION:  threshold_adaptation
 *
 *     VERSION:   1.2
 *
 *     PURPOSE:   Evaluates the secondary VAD decision.  If speech is not
 *                present then the noise model rvad and adaptive threshold
 *                thvad are updated.
 *
 *     INPUTS:    swStat        flag to indicate spectral stationarity
 *                swPtch        flag to indicate a periodic signal component
 *                swTone        flag to indicate a tone signal component
 *                pswRav1[0..8] ACF obtained from l_av1
 *                swNormRav1    r_av1 scaling factor
 *                swM_pvad      mantissa of filtered signal energy
 *                swE_pvad      exponent of filtered signal energy
 *                swM_acf0      mantissa of signal frame energy
 *                swE_acf0      exponent of signal frame energy
 *
 *     OUTPUTS:   pswRvad[0..8] autocorrelated adaptive filter coefficients
 *                pswNormRvad   rvad scaling factor
 *                pswM_thvad    mantissa of decision threshold
 *                pswE_thvad    exponent of decision threshold
 *
 ***************************************************************************/

void   threshold_adaptation(Shortword swStat,
                                   Shortword swPtch,
                                   Shortword swTone,
                                   Shortword pswRav1[],
                                   Shortword swNormRav1,
                                   Shortword swM_pvad,
                                   Shortword swE_pvad,
                                   Shortword swM_acf0,
                                   Shortword swE_acf0,
                                   Shortword pswRvad[],
                                   Shortword *pswNormRvad,
                                   Shortword *pswM_thvad,
                                   Shortword *pswE_thvad)
{

/*_________________________________________________________________________
 |                                                                         |
 |                            Automatic Variables                          |
 |_________________________________________________________________________|
*/

  Longword
         L_temp;

  Shortword
         swTemp,
         swComp,
         swComp2,
         swM_temp,
         swE_temp;

  int
         i;


/*_________________________________________________________________________
 |                                                                         |
 |                              Executable Code                            |
 |_________________________________________________________________________|
*/

  swComp = 0;

  /*** Test if acf0 < pth; if yes set thvad to plev ***/

  if (swE_acf0 < E_PTH)
    swComp = 1;
  if ((swE_acf0 == E_PTH) && (swM_acf0 < M_PTH))
    swComp = 1;

  if (swComp == 1)
  {
    *pswE_thvad = E_PLEV;
    *pswM_thvad = M_PLEV;

    return;
  }


  /*** Test if an adaption is required ***/

  if (swPtch == 1)
    swComp = 1;
  if (swStat == 0)
    swComp = 1;
  if (swTone == 1)
    swComp = 1;

  if (swComp == 1)
  {
    swAdaptCount = 0;
    return;
  }


  /*** Increment adaptcount ***/

  swAdaptCount = add(swAdaptCount, 1);
  if (swAdaptCount <= 8)
    return;


  /*** computation of thvad-(thvad/dec) ***/

  *pswM_thvad = sub(*pswM_thvad, shr(*pswM_thvad, 5));

  if (*pswM_thvad < 0x4000)
  {
    *pswM_thvad = shl(*pswM_thvad, 1);
    *pswE_thvad = sub(*pswE_thvad, 1);
  }


  /*** computation of pvad*fac ***/

  L_temp = L_mult(swM_pvad, 20889);
  L_temp = L_shr(L_temp, 15);
  swE_temp = add(swE_pvad, 1);

  if (L_temp > 0x7fffL)
  {
    L_temp = L_shr(L_temp, 1);
    swE_temp = add(swE_temp, 1);
  }
  swM_temp = extract_l(L_temp);


  /*** test if thvad < pavd*fac ***/

  if (*pswE_thvad < swE_temp)
    swComp = 1;

  if ((*pswE_thvad == swE_temp) && (*pswM_thvad < swM_temp))
    swComp = 1;


  /*** compute minimum(thvad+(thvad/inc), pvad*fac) when comp = 1 ***/

  if (swComp == 1)
  {

    /*** compute thvad + (thvad/inc) ***/

    L_temp = L_add(L_deposit_l(*pswM_thvad),L_deposit_l(shr(*pswM_thvad, 4)));

    if (L_temp > 0x7fffL)
    {
      *pswM_thvad = extract_l(L_shr(L_temp, 1));
      *pswE_thvad = add(*pswE_thvad, 1);
    }
    else
      *pswM_thvad = extract_l(L_temp);

    swComp2 = 0;

    if (swE_temp < *pswE_thvad)
      swComp2 = 1;

    if ((swE_temp == *pswE_thvad) && (swM_temp < *pswM_thvad))
      swComp2 = 1;

    if (swComp2 == 1)
    {
      *pswE_thvad = swE_temp;
      *pswM_thvad = swM_temp;
    }
  }


  /*** compute pvad + margin ***/

  if (swE_pvad == E_MARGIN)
  {
    L_temp = L_add(L_deposit_l(swM_pvad), L_deposit_l(M_MARGIN));
    swM_temp = extract_l(L_shr(L_temp, 1));
    swE_temp = add(swE_pvad, 1);
  }
  else
  {
    if (swE_pvad > E_MARGIN)
    {
      swTemp = sub(swE_pvad, E_MARGIN);
      swTemp = shr(M_MARGIN, swTemp);
      L_temp = L_add(L_deposit_l(swM_pvad), L_deposit_l(swTemp));

      if (L_temp > 0x7fffL)
      {
        swE_temp = add(swE_pvad, 1);
        swM_temp = extract_l(L_shr(L_temp, 1));
      }
      else
      {
        swE_temp = swE_pvad;
        swM_temp = extract_l(L_temp);
      }
    }
    else
    {
      swTemp = sub(E_MARGIN, swE_pvad);
      swTemp = shr(swM_pvad, swTemp);
      L_temp = L_add(L_deposit_l(M_MARGIN), L_deposit_l(swTemp));

      if (L_temp > 0x7fffL)
      {
        swE_temp = add(E_MARGIN, 1);
        swM_temp = extract_l(L_shr(L_temp, 1));
      }
      else
      {
        swE_temp = E_MARGIN;
        swM_temp = extract_l(L_temp);
      }
    }
  }

  /*** Test if thvad > pvad + margin ***/

  swComp = 0;

  if (*pswE_thvad > swE_temp)
    swComp = 1;

  if ((*pswE_thvad == swE_temp) && (*pswM_thvad > swM_temp))
    swComp = 1;

  if (swComp == 1)
  {
    *pswE_thvad = swE_temp;
    *pswM_thvad = swM_temp;
  }

  /*** Normalise and retain rvad[0..8] in memory ***/

  *pswNormRvad = swNormRav1;

  for (i = 0; i <= 8; i++)
    pswRvad[i] = pswRav1[i];

  /*** Set adaptcount to adp + 1 ***/

  swAdaptCount = 9;

}

/****************************************************************************
 *
 *     FUNCTION:  vad_decision
 *
 *     VERSION:   1.2
 *
 *     PURPOSE:   Computes the VAD decision based on the comparison of the
 *                floating point representations of pvad and thvad.
 *
 *     INPUTS:    swM_pvad      mantissa of filtered signal energy
 *                swE_pvad      exponent of filtered signal energy
 *                swM_thvad     mantissa of decision threshold
 *                swE_thvad     exponent of decision threshold
 *
 *     OUTPUTS:   pswVvad       vad decision before hangover is added
 *
 ***************************************************************************/

void   vad_decision(Shortword swM_pvad,
                           Shortword swE_pvad,
                           Shortword swM_thvad,
                           Shortword swE_thvad,
                           Shortword *pswVvad)
{

/*_________________________________________________________________________
 |                                                                         |
 |                              Executable Code                            |
 |_________________________________________________________________________|
*/

  *pswVvad = 0;

  if (swE_pvad > swE_thvad)
    *pswVvad = 1;
  if ((swE_pvad == swE_thvad) && (swM_pvad > swM_thvad))
    *pswVvad = 1;

}

/****************************************************************************
 *
 *     FUNCTION:  vad_hangover
 *
 *     VERSION:   1.2
 *
 *     PURPOSE:   Computes the final VAD decision for the current frame
 *                being processed.
 *
 *     INPUTS:    swVvad        vad decision before hangover is added
 *
 *     OUTPUTS:   pswVadFlag    vad decision after hangover is added
 *
 ***************************************************************************/

void   vad_hangover(Shortword swVvad,
                           Shortword *pswVadFlag)
{

/*_________________________________________________________________________
 |                                                                         |
 |                              Executable Code                            |
 |_________________________________________________________________________|
*/

  if (swVvad == 1)
    swBurstCount = add(swBurstCount, 1);
  else
    swBurstCount = 0;

  if (swBurstCount >= 3)
  {
    swHangCount = 5;
    swBurstCount = 3;
  }

  *pswVadFlag = swVvad;

  if (swHangCount >= 0)
  {
    *pswVadFlag = 1;
    swHangCount = sub(swHangCount, 1);
  }

}

/****************************************************************************
 *
 *     FUNCTION:  periodicity_update
 *
 *     VERSION:   1.2
 *
 *     PURPOSE:   Computes the ptch flag needed for the threshold
 *                adaptation decision for the next frame.
 *
 *     INPUTS:    pswLags[0..3]    speech encoder long term predictor lags
 *
 *     OUTPUTS:   pswPtch          Boolean voiced / unvoiced decision
 *
 ***************************************************************************/

void   periodicity_update(Shortword pswLags[4],
                                 Shortword *pswPtch)
{

/*_________________________________________________________________________
 |                                                                         |
 |                            Automatic Variables                          |
 |_________________________________________________________________________|
*/

  Shortword
         swMinLag,
         swMaxLag,
         swSmallLag,
         swLagCount,
         swTemp;

  int
         i,
         j;

/*_________________________________________________________________________
 |                                                                         |
 |                              Executable Code                            |
 |_________________________________________________________________________|
*/

  /*** Run loop for No. of sub-segments in the frame ***/

  swLagCount = 0;

  for (i = 0; i <= 3; i++)
  {
    /*** Search the maximum and minimum of consecutive lags ***/

    if (swOldLag > pswLags[i])
    {
      swMinLag = pswLags[i];
      swMaxLag = swOldLag;
    }
    else
    {
      swMinLag = swOldLag;
      swMaxLag = pswLags[i];
    }


    /*** Compute smallag (modulo operation not defined) ***/

    swSmallLag = swMaxLag;

    for (j = 0; j <= 2; j++)
    {
      if (swSmallLag >= swMinLag)
        swSmallLag = sub(swSmallLag, swMinLag);
    }


    /***  Minimum of smallag and minlag - smallag ***/

    swTemp = sub(swMinLag, swSmallLag);

    if (swTemp < swSmallLag)
      swSmallLag = swTemp;

    if (swSmallLag < 2)
      swLagCount = add(swLagCount, 1);


    /*** Save the current LTP lag ***/

    swOldLag = pswLags[i];
  }


  /*** Update the veryoldlagcount and oldlagcount ***/

  swVeryOldLagCount = swOldLagCount;
  swOldLagCount = swLagCount;


  /*** Make ptch decision ready for next frame ***/

  swTemp = add(swOldLagCount, swVeryOldLagCount);

  if (swTemp >= 7)
    *pswPtch = 1;
  else
    *pswPtch = 0;

}
