#ifndef __VAD
#define __VAD

#include "typedefs.h"


/*_________________________________________________________________________
 |                                                                         |
 |                            Function Prototypes                          |
 |_________________________________________________________________________|
*/

void   vad_reset(void);

void   vad_algorithm
       (
               Longword pL_acf[9],
               Shortword swScaleAcf,
               Shortword pswRc[4],
               Shortword swPtch,
               Shortword *pswVadFlag
);

void   energy_computation
       (
               Longword pL_acf[],
               Shortword swScaleAcf,
               Shortword pswRvad[],
               Shortword swNormRvad,
               Shortword *pswM_pvad,
               Shortword *pswE_pvad,
               Shortword *pswM_acf0,
               Shortword *pswE_acf0
);


void   average_acf
       (
               Longword pL_acf[],
               Shortword swScaleAcf,
               Longword pL_av0[],
               Longword pL_av1[]
);

void   predictor_values
       (
               Longword pL_av1[],
               Shortword pswRav1[],
               Shortword *pswNormRav1
);

void   schur_recursion
       (
               Longword pL_av1[],
               Shortword pswVpar[]
);

void   step_up
       (
               Shortword swNp,
               Shortword pswVpar[],
               Shortword pswAav1[]
);

void   compute_rav1
       (
               Shortword pswAav1[],
               Shortword pswRav1[],
               Shortword *pswNormRav1
);

void   spectral_comparison
       (
               Shortword pswRav1[],
               Shortword swNormRav1,
               Longword pL_av0[],
               Shortword *pswStat
);

void   tone_detection
       (
               Shortword pswRc[4],
               Shortword *pswTone
);


void   threshold_adaptation
       (
               Shortword swStat,
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
               Shortword *pswE_thvad
);

void   vad_decision
       (
               Shortword swM_pvad,
               Shortword swE_pvad,
               Shortword swM_thvad,
               Shortword swE_thvad,
               Shortword *pswVvad
);

void   vad_hangover
       (
               Shortword swVvad,
               Shortword *pswVadFlag
);

void   periodicity_update
       (
               Shortword pswLags[4],
               Shortword *pswPtch
);

#endif
