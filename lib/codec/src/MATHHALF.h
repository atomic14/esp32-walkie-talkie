#ifndef __MATHHALF
#define __MATHHALF

#include "typedefs.h"

/*_________________________________________________________________________
 |                                                                         |
 |                            Function Prototypes                          |
 |_________________________________________________________________________|
*/

/* addition */
/************/

Shortword add(Shortword var1, Shortword var2);  /* 1 ops */
Shortword sub(Shortword var1, Shortword var2);  /* 1 ops */
Longword L_add(Longword L_var1, Longword L_var2);       /* 2 ops */
Longword L_sub(Longword L_var1, Longword L_var2);       /* 2 ops */

/* multiplication */
/******************/

Shortword mult(Shortword var1, Shortword var2); /* 1 ops */
Longword L_mult(Shortword var1, Shortword var2);        /* 1 ops */
Shortword mult_r(Shortword var1, Shortword var2);       /* 2 ops */


/* arithmetic shifts */
/*********************/

Shortword shr(Shortword var1, Shortword var2);  /* 1 ops */
Shortword shl(Shortword var1, Shortword var2);  /* 1 ops */
Longword L_shr(Longword L_var1, Shortword var2);        /* 2 ops */
Longword L_shl(Longword L_var1, Shortword var2);        /* 2 ops */
Shortword shift_r(Shortword var, Shortword var2);       /* 2 ops */
Longword L_shift_r(Longword L_var, Shortword var2);     /* 3 ops */

/* absolute value  */
/*******************/

Shortword abs_s(Shortword var1);       /* 1 ops */
Longword L_abs(Longword var1);         /* 3 ops */


/* multiply accumulate  */
/************************/

Longword L_mac(Longword L_var3,
                      Shortword var1, Shortword var2);  /* 1 op */
Shortword mac_r(Longword L_var3,
                       Shortword var1, Shortword var2); /* 2 op */
Longword L_msu(Longword L_var3,
                      Shortword var1, Shortword var2);  /* 1 op */
Shortword msu_r(Longword L_var3,
                       Shortword var1, Shortword var2); /* 2 op */

/* negation  */
/*************/

Shortword negate(Shortword var1);      /* 1 ops */
Longword L_negate(Longword L_var1);    /* 2 ops */


/* Accumulator manipulation */
/****************************/

Longword L_deposit_l(Shortword var1);  /* 1 ops */
Longword L_deposit_h(Shortword var1);  /* 1 ops */
Shortword extract_l(Longword L_var1);  /* 1 ops */
Shortword extract_h(Longword L_var1);  /* 1 ops */

/* Round */
/*********/

Shortword round(Longword L_var1);      /* 1 ops */

/* Normalization */
/*****************/

Shortword norm_l(Longword L_var1);     /* 30 ops */
Shortword norm_s(Shortword var1);      /* 15 ops */

/* Division */
/************/
Shortword divide_s(Shortword var1, Shortword var2);     /* 18 ops */

/* Non-saturating instructions */
/*******************************/
Longword L_add_c(Longword L_Var1, Longword L_Var2);     /* 2 ops */
Longword L_sub_c(Longword L_Var1, Longword L_Var2);     /* 2 ops */
Longword L_sat(Longword L_var1);       /* 4 ops */
Longword L_macNs(Longword L_var3,
                        Shortword var1, Shortword var2);        /* 1 ops */
Longword L_msuNs(Longword L_var3,
                        Shortword var1, Shortword var2);        /* 1 ops */

#endif
