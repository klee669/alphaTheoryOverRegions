//
//  MPFI.cpp
//  Homotopy Continuation
//
//  Created by burr2 on 12/3/18.
//  Copyright © 2018 burr2. All rights reserved.
//

#include "MPFI.hpp"

namespace homotopy
{
    ///////////////////////
    // Declaring static members
    ///////////////////////
    MPFR MPFI::temp1;
    MPFR MPFI::temp2;

    // Taken from mpfi.h.  Changed to use static mpfr's instead of redeclaring.
    int my_mpfi_mul (MPFI& a, const MPFI& b, const MPFI& c)
    {
        int inexact_left, inexact_right;
        int inexact = 0;
        
        /* Handling the NaN cases */
        if ( (mpfr_nan_p(&(*b.interval).left) || (mpfr_nan_p(&(*b.interval).right))) ||
            (mpfr_nan_p(&(*c.interval).left) || (mpfr_nan_p(&(*c.interval).right))))
        {
            mpfr_set_nan (&(*a.interval).left);
            mpfr_set_nan (&(*a.interval).right);
            do {mpfr_set_nanflag(); return 0;} while (0);
        }
        
        /* Handling the case where one operand is 0, in order */
        /* to avoid problems with 0 * an infinite interval    */
        if (((mpfr_nan_p(&(*b.interval).left) || (mpfr_nan_p(&(*b.interval).right))) ?
             0 : ((mpfr_sgn(&((*b.interval).right))==0) && (mpfr_sgn(&((*b.interval).left))==0)))) {
            return (mpfi_set (a.interval, b.interval));
        }
        if (((mpfr_nan_p(&(*c.interval).left) || (mpfr_nan_p(&(*c.interval).right))) ?
             0 : ((mpfr_sgn(&((*c.interval).right))==0) && (mpfr_sgn(&((*c.interval).left))==0)))) {
            return (mpfi_set (a.interval, c.interval));
        }
        
        if (mpfr_sgn (&((*b.interval).left)) >= 0) {
            if (mpfr_sgn (&((*c.interval).left)) >=0) {
                /* b nonnegative and c nonnegative */
                inexact_left = mpfr_mul (&((*a.interval).left), &((*b.interval).left), &((*c.interval).left), MPFR_RNDD);
                inexact_right = mpfr_mul (&((*a.interval).right), &((*b.interval).right), &((*c.interval).right), MPFR_RNDU);
            }
            else {
                set_precision(a.temp1,mpfr_get_prec (&((*a.interval).left)));
                if (mpfr_sgn (&((*c.interval).right)) <= 0) {
                    /* b nonnegative and c non-positive */
                    inexact_left = mpfr_mul (a.temp1.point, &((*b.interval).right), &((*c.interval).left), MPFR_RNDD);
                    inexact_right = mpfr_mul (&((*a.interval).right), &((*b.interval).left), &((*c.interval).right), MPFR_RNDU);
                }
                else {
                    /* b nonnegative and c overlapping 0 */
                    inexact_left = mpfr_mul (a.temp1.point, &((*b.interval).right), &((*c.interval).left), MPFR_RNDD);
                    inexact_right = mpfr_mul (&((*a.interval).right), &((*b.interval).right), &((*c.interval).right), MPFR_RNDU);
                }
                mpfr_set (&((*a.interval).left), a.temp1.point, MPFR_RNDD); /* exact */
            }
        }
        else {
            if (mpfr_sgn (&((*b.interval).right)) <= 0) {
                /* b non-positive */
                set_precision(a.temp1,mpfr_get_prec (&((*a.interval).left)));
                if (mpfr_sgn (&((*c.interval).left)) >= 0) {
                    /* b non-positive and c nonnegative */
                    inexact_left = mpfr_mul (a.temp1.point, &((*b.interval).left), &((*c.interval).right), MPFR_RNDD);
                    inexact_right = mpfr_mul (&((*a.interval).right), &((*b.interval).right), &((*c.interval).left), MPFR_RNDU);
                }
                else {
                    if (mpfr_sgn (&((*c.interval).right)) <= 0) {
                        /* b non-positive and c non-positive */
                        inexact_left = mpfr_mul (a.temp1.point, &((*b.interval).right), &((*c.interval).right), MPFR_RNDD);
                        inexact_right = mpfr_mul (&((*a.interval).right), &((*b.interval).left), &((*c.interval).left), MPFR_RNDU);
                    }
                    else {
                        /* b non-positive and c overlapping 0 */
                        inexact_left = mpfr_mul (a.temp1.point, &((*b.interval).left), &((*c.interval).right), MPFR_RNDD);
                        inexact_right = mpfr_mul (&((*a.interval).right), &((*b.interval).left), &((*c.interval).left), MPFR_RNDU);
                    }
                }
                mpfr_set (&((*a.interval).left), a.temp1.point, MPFR_RNDD); /* exact */
            }
            else {
                /* b contains 0 */
                if (mpfr_sgn (&((*c.interval).left)) >= 0) {
                    /* b overlapping 0 and c nonnegative  */
                    set_precision(a.temp1,mpfr_get_prec (&((*a.interval).left)));
                    
                    inexact_left = mpfr_mul (a.temp1.point, &((*b.interval).left), &((*c.interval).right), MPFR_RNDD);
                    inexact_right = mpfr_mul (&((*a.interval).right), &((*b.interval).right), &((*c.interval).right), MPFR_RNDU);
                    
                    mpfr_set (&((*a.interval).left), a.temp1.point, MPFR_RNDD);
                }
                else {
                    if (mpfr_sgn (&((*c.interval).right)) <= 0) {
                        /* b overlapping 0 and c non-positive */
                        set_precision(a.temp1,mpfr_get_prec (&((*a.interval).left)));
                        
                        inexact_left = mpfr_mul (a.temp1.point, &((*b.interval).right), &((*c.interval).left), MPFR_RNDD);
                        inexact_right = mpfr_mul (&((*a.interval).right), &((*b.interval).left), &((*c.interval).left), MPFR_RNDU);
                        
                        mpfr_set (&((*a.interval).left), a.temp1.point, MPFR_RNDD);
                    }
                    else {
                        /* b overlapping 0 and c overlapping 0
                         Beware the case where the result is one of the operands! */
                        int inexact_tmp;
                        
                        set_precision(a.temp1,mpfr_get_prec (&((*a.interval).left)));
                        set_precision(a.temp2,mpfr_get_prec (&((*a.interval).left)));
                        
                        inexact_tmp = mpfr_mul (a.temp1.point, &((*b.interval).left), &((*c.interval).right), MPFR_RNDD);
                        inexact_left = mpfr_mul (a.temp2.point, &((*b.interval).right), &((*c.interval).left), MPFR_RNDD);
                        if (mpfr_cmp (a.temp1.point, a.temp2.point) < 0) {
                            mpfr_swap (a.temp2.point, a.temp1.point); /* same precision */
                            inexact_left = inexact_tmp;
                        }
                        
                        set_precision(a.temp1,mpfr_get_prec (&((*a.interval).left)));
                        inexact_tmp = mpfr_mul (a.temp1.point, &((*b.interval).left), &((*c.interval).left), MPFR_RNDU);
                        inexact_right = mpfr_mul (&((*a.interval).right), &((*b.interval).right), &((*c.interval).right), MPFR_RNDU);
                        if (mpfr_cmp (a.temp1.point, &((*a.interval).right)) > 0) {
                            mpfr_set (&((*a.interval).right),a.temp1.point, MPFR_RNDU); /* exact */
                            inexact_right = inexact_tmp;
                        }
                        mpfr_set (&((*a.interval).left), a.temp2.point, MPFR_RNDD); /* exact */
                    }
                }
            }
        }
        
        /* no need to check to sign of an endpoint equal to 0, it should be OK */
        if (inexact_left)
            inexact += 1;
        if (inexact_right)
            inexact += 2;
        
        return inexact;
    }
}
