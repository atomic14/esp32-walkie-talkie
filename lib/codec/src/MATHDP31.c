/***************************************************************************
 *
 *   File Name:  mathdp31.c
 *
 *   Purpose:  Contains functions increased-precision arithmetic operations.
 *
 *      Below is a listing of all the functions in this file.  There
 *      is no interdependence among the functions.
 *
 *      L_mpy_ls()
 *      L_mpy_ll()
 *      isLwLimit()
 *      isSwLimit()
 *
 ***************************************************************************/
/*_________________________________________________________________________
 |                                                                         |
 |                              Include Files                              |
 |_________________________________________________________________________|
*/

#include "mathhalf.h"
#include "typedefs.h"

/****************************************************************************
 *
 *     FUNCTION NAME: isLwLimit
 *
 *     PURPOSE:
 *
 *        Check to see if the input Longword is at the
 *        upper or lower limit of its range.  i.e.
 *        0x7fff ffff or -0x8000 0000
 *
 *        Ostensibly this is a check for an overflow.
 *        This does not truly mean an overflow occurred,
 *        it means the value reached is the
 *        maximum/minimum value representable.  It may
 *        have come about due to an overflow.
 *
 *     INPUTS:
 *
 *       L_In               A Longword input variable
 *
 *
 *     OUTPUTS:             none
 *
 *     RETURN VALUE:        1 if input == 0x7fff ffff or -0x8000 0000
 *                          0 otherwise
 *
 *     KEYWORDS: saturation, limit
 *
 ***************************************************************************/

short  isLwLimit(Longword L_In)
{

  Longword L_ls;
  short  siOut;

  if (L_In != 0)
  {
    L_ls = L_shl(L_In, 1);
    if (L_sub(L_In, L_ls) == 0)
      siOut = 1;
    else
      siOut = 0;
  }
  else
  {
    siOut = 0;
  }
  return (siOut);
}

/****************************************************************************
 *
 *     FUNCTION NAME: isSwLimit
 *
 *     PURPOSE:
 *
 *        Check to see if the input Shortword is at the
 *        upper or lower limit of its range.  i.e.
 *        0x7fff or -0x8000
 *
 *        Ostensibly this is a check for an overflow.
 *        This does not truly mean an overflow occurred,
 *        it means the value reached is the
 *        maximum/minimum value representable.  It may
 *        have come about due to an overflow.
 *
 *     INPUTS:
 *
 *       swIn               A Shortword input variable
 *
 *
 *     OUTPUTS:             none
 *
 *     RETURN VALUE:        1 if input == 0x7fff or -0x8000
 *                          0 otherwise
 *
 *     KEYWORDS: saturation, limit
 *
 ***************************************************************************/

short  isSwLimit(Shortword swIn)
{

  Shortword swls;
  short  siOut;

  if (swIn != 0)
  {
    swls = shl(swIn, 1);
    if (sub(swIn, swls) == 0)          /* logical compare outputs 1/0 */
      siOut = 1;
    else
      siOut = 0;
  }
  else
  {
    siOut = 0;
  }
  return (siOut);

}

/****************************************************************************
 *
 *     FUNCTION NAME: L_mpy_ll
 *
 *     PURPOSE:    Multiply a 32 bit number (L_var1) and a 32 bit number
 *                 (L_var2), and return a 32 bit result.
 *
 *     INPUTS:
 *
 *       L_var1             A Longword input variable
 *
 *       L_var2             A Longword input variable
 *
 *     OUTPUTS:             none
 *
 *     IMPLEMENTATION:
 *
 *        Performs a 31x31 bit multiply, Complexity=24 Ops.
 *
 *        Let x1x0, or y1y0, be the two constituent halves
 *        of a 32 bit number.  This function performs the
 *        following:
 *
 *        low = ((x0 >> 1)*(y0 >> 1)) >> 16     (low * low)
 *        mid1 = [(x1 * (y0 >> 1)) >> 1 ]       (high * low)
 *        mid2 = [(y1 * (x0 >> 1)) >> 1]        (high * low)
 *        mid =  (mid1 + low + mid2) >> 14      (sum so far)
 *        output = (y1*x1) + mid                (high * high)
 *
 *
 *     RETURN VALUE:        A Longword value
 *
 *     KEYWORDS: mult,mpy,multiplication
 *
 ***************************************************************************/

Longword L_mpy_ll(Longword L_var1, Longword L_var2)
{
  Shortword swLow1,
         swLow2,
         swHigh1,
         swHigh2;
  Longword L_varOut,
         L_low,
         L_mid1,
         L_mid2,
         L_mid;


  swLow1 = shr(extract_l(L_var1), 1);
  swLow1 = SW_MAX & swLow1;

  swLow2 = shr(extract_l(L_var2), 1);
  swLow2 = SW_MAX & swLow2;
  swHigh1 = extract_h(L_var1);
  swHigh2 = extract_h(L_var2);

  L_low = L_mult(swLow1, swLow2);
  L_low = L_shr(L_low, 16);

  L_mid1 = L_mult(swHigh1, swLow2);
  L_mid1 = L_shr(L_mid1, 1);
  L_mid = L_add(L_mid1, L_low);

  L_mid2 = L_mult(swHigh2, swLow1);
  L_mid2 = L_shr(L_mid2, 1);
  L_mid = L_add(L_mid, L_mid2);

  L_mid = L_shr(L_mid, 14);
  L_varOut = L_mac(L_mid, swHigh1, swHigh2);

  return (L_varOut);
}

/****************************************************************************
 *
 *     FUNCTION NAME: L_mpy_ls
 *
 *     PURPOSE:    Multiply a 32 bit number (L_var2) and a 16 bit
 *                 number (var1) returning a 32 bit result. L_var2
 *                 is truncated to 31 bits prior to executing the
 *                 multiply.
 *
 *     INPUTS:
 *
 *       L_var2             A Longword input variable
 *
 *       var1               A Shortword input variable
 *
 *     OUTPUTS:             none
 *
 *     RETURN VALUE:        A Longword value
 *
 *     KEYWORDS: mult,mpy,multiplication
 *
 ***************************************************************************/

Longword L_mpy_ls(Longword L_var2, Shortword var1)
{
  Longword L_varOut;
  Shortword swtemp;

  swtemp = shr(extract_l(L_var2), 1);
  swtemp = (short) 32767 & (short) swtemp;

  L_varOut = L_mult(var1, swtemp);
  L_varOut = L_shr(L_varOut, 15);
  L_varOut = L_mac(L_varOut, var1, extract_h(L_var2));
  return (L_varOut);
}
