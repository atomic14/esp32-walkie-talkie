/***************************************************************************
 *
 *   File Name:  globdefs.c
 *
 *   Purpose:  Defines global frame count, subframe count,
 *             and DTX on/off switch variables.
 *
 **************************************************************************/
/*
  Storage of global variables (definition)
  to reference (declaration) these you must include "typeDefs.h"
  */

int    giFrmCnt;                       /* Frame count: 0,1,2,3,4,... */
int    giSfrmCnt = 0;                  /* Subframe count: 0,1,2,3 */


/*
  Global flag which indicates: DTX mode switched on/off
  1 - DTX mode switched on (speech encoder)
  0 - DTX mode switched off (speech encoder)
  */

int    giDTXon = 1;
