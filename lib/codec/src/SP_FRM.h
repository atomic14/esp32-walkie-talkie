#ifndef __SP_FRM
#define __SP_FRM

#include "typedefs.h"
#include "sp_rom.h"


struct QuantList
{
  /* structure which points to the beginning of a block of candidate vq
   * vectors.  It also stores the residual error for each vector. */
  int    iNum;                         /* total number in list */
  int    iRCIndex;                     /* an index to the first vector of the
                                        * block */
  Shortword pswPredErr[PREQ1_NUM_OF_ROWS];  /* PREQ1 is the biggest block */
};

/*_________________________________________________________________________
 |                                                                         |
 |                            Function Prototypes                          |
 |_________________________________________________________________________|
*/

void   iir_d(Shortword pswCoeff[], Shortword pswIn[],
                    Shortword pswXstate[],
                    Shortword pswYstate[],
                    int npts, int shifts,
                    Shortword swPreFirDownSh,
                    Shortword swFinalUpShift);


  void   filt4_2nd(Shortword pswCoeff[],
                          Shortword pswIn[],
                          Shortword pswXstate[],
                          Shortword pswYstate[],
                          int npts,
                          int shifts);

  void   initPBarVBarL(Longword pL_PBarFull[],
                              Shortword pswPBar[],
                              Shortword pswVBar[]);

  void   initPBarFullVBarFullL(Longword pL_CorrelSeq[],
                                      Longword pL_PBarFull[],
                                      Longword pL_VBarFull[]);

  Shortword aflatRecursion(Shortword pswQntRc[],
                                  Shortword pswPBar[],
                                  Shortword pswVBar[],
                                  Shortword *ppswPAddrs[],
                                  Shortword *ppswVAddrs[],
                                  Shortword swSegmentOrder);

  void   aflatNewBarRecursionL(Shortword pswQntRc[],
                                      int iSegment,
                                      Longword pL_PBar[],
                                      Longword pL_VBar[],
                                      Shortword pswPBar[],
                                      Shortword pswVBar[]);


  void   setupPreQ(int iSeg, int iVector);

  void   setupQuant(int iSeg, int iVector);

  void   getNextVec(Shortword pswRc[]);

  void   aflat(Shortword pswSpeechToLPC[],
                      int piR0Index[],
                      Shortword pswFinalRc[],
                      int piVQCodewds[],
                      Shortword swPtch,
                      Shortword *pswVadFlag,
                      Shortword *pswSP);


  Shortword fnExp2(Longword L_Input);

  Shortword fnLog2(Longword L_Input);

  void   weightSpeechFrame(Shortword pswSpeechFrm[],
                                  Shortword pswWNumSpace[],
                                  Shortword pswWDenomSpace[],
                                  Shortword pswWSpeechBuffBase[]);

  void   getSfrmLpcTx(Shortword swPrevR0, Shortword swNewR0,
                             Shortword pswPrevFrmKs[],
                             Shortword pswPrevFrmAs[],
                             Shortword pswPrevFrmSNWCoef[],
                             Shortword pswNewFrmKs[],
                             Shortword pswNewFrmAs[],
                             Shortword pswNewFrmSNWCoef[],
                             Shortword pswHPFSpeech[],
                             short *pswSoftInterp,
                             struct NormSw *psnsSqrtRs,
                             Shortword ppswSynthAs[][NP],
                             Shortword ppswSNWCoefAs[][NP]);

  short int fnBest_CG(Shortword pswCframe[],
                             Shortword pswGframe[],
                             Shortword *pswCmaxSqr,
                             Shortword *pswGmax,
                             short int siNumPairs);

  short  compResidEnergy(Shortword pswSpeech[],
                                Shortword ppswInterpCoef[][NP],
                                Shortword pswPreviousCoef[],
                                Shortword pswCurrentCoef[],
                                struct NormSw psnsSqrtRs[]);

  Shortword r0Quant(Longword L_UnqntzdR0);

  Shortword cov32(Shortword pswIn[],
                         Longword pppL_B[NP][NP][2],
                         Longword pppL_F[NP][NP][2],
                         Longword pppL_C[NP][NP][2],
                         Longword *pL_R0,
                         Longword pL_VadAcf[],
                         Shortword *pswVadScalAuto);

  Longword flat(Shortword pswSpeechIn[],
                       Shortword pswRc[],
                       int *piR0Inx,
                       Longword pL_VadAcf[],
                       Shortword *pswVadScalAuto);



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
                                  Shortword swSP);

  Shortword getCCThreshold(Shortword swRp0,
                                  Shortword swCC,
                                  Shortword swG);

  void   pitchLags(Shortword swBestIntLag,
                          Shortword pswIntCs[],
                          Shortword pswIntGs[],
                          Shortword swCCThreshold,
                          Shortword pswLPeaksSorted[],
                          Shortword pswCPeaksSorted[],
                          Shortword pswGPeaksSorted[],
                          Shortword *psiNumSorted,
                          Shortword *pswPitch,
                          Shortword *pswHNWCoef);

  short  CGInterpValid(Shortword swFullResLag,
                              Shortword pswCIn[],
                              Shortword pswGIn[],
                              Shortword pswLOut[],
                              Shortword pswCOut[],
                              Shortword pswGOut[]);

  void   CGInterp(Shortword pswLIn[],
                         short siNum,
                         Shortword pswCIn[],
                         Shortword pswGIn[],
                         short siLoIntLag,
                         Shortword pswCOut[],
                         Shortword pswGOut[]);

  Shortword quantLag(Shortword swRawLag,
                            Shortword *psiCode);

  void   findBestInQuantList(struct QuantList psqlInList,
                                    int iNumVectOut,
                                    struct QuantList psqlBestOutList[]);

  Shortword findPeak(Shortword swSingleResLag,
                            Shortword pswCIn[],
                            Shortword pswGIn[]);

  void   bestDelta(Shortword pswLagList[],
                          Shortword pswCSfrm[],
                          Shortword pswGSfrm[],
                          short int siNumLags,
                          short int siSfrmIndex,
                          Shortword pswLTraj[],
                          Shortword pswCCTraj[],
                          Shortword pswGTraj[]);

  Shortword
         maxCCOverGWithSign(Shortword pswCIn[],
                                   Shortword pswGIn[],
                                   Shortword *pswCCMax,
                                   Shortword *pswGMax,
                                   Shortword swNum);

  void   getNWCoefs(Shortword pswACoefs[],
                           Shortword pswHCoefs[]);

#endif
