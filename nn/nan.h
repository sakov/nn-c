/******************************************************************************
 *
 * File:           nan.h
 *
 * Created:        18/10/2001
 *
 * Author:         Pavel Sakov
 *                 CSIRO Marine Research
 *
 * Purpose:        NaN definition
 *
 * Description:    Should cover machines with 64 bit doubles or other machines
 *                 with GCC
 *
 * Revisions:      10 July 2026: a more robust version
 *
 *****************************************************************************/

#if !defined(_NAN_H)
#define _NAN_H

#include <math.h>

#if defined(NAN)
#define NaN NAN

#elif defined(__GNUC__) || defined(__clang__)
#define NaN __builtin_nan("")

#elif defined(_MSC_VER)
#include <float.h>
#define NaN (_NAN)

#else
#define NaN (0.0 / 0.0)

#endif

#endif /* _NAN_H */
