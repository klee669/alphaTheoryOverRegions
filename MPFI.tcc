//
//  MPFI.tcc
//  Homotopy Continuation
//
//  Created by burr2 on 8/10/18.
//  Copyright © 2018 burr2. All rights reserved.
//

#ifndef MPFI_tcc
#define MPFI_tcc

#include "MPFI.hpp"
#include "NumberUtils.hpp"
#include "MPFI.cpp"

namespace homotopy
{
    ///////////////////////
    // Constructors
    ///////////////////////
    
    // Default constructor.  Sets everything equal to 0 using default precision.
    MPFI::MPFI()
    {
        mpfi_init(interval);
        mpfi_set_si(interval,0);
        initialized = true;
    }

    // Constructor with precision.  Sets everything equal to 0 using input precision.
    MPFI::MPFI(mpfr_prec_t prec)
    {
        mpfi_init2(interval,prec);
        mpfi_set_si(interval,0);
        initialized = true;
    }

    // Copy constructor.  Sets upper and lower equal to values in input.
    // Sets precision level equal to that of the input.
    MPFI::MPFI(const MPFI &mpfi)
    {
        mpfi_init2(interval,mpfi_get_prec(mpfi.interval));
        mpfi_set(interval,mpfi.interval);
        initialized = true;
    }
    
    // Copy constructor with precision.  Copies value of one point to another.
    // Sets precision level equal to the user specified precision.
    MPFI::MPFI(const MPFI &mpfi, mpfr_prec_t prec)
    {
        mpfi_init2(interval,prec);
        mpfi_set(interval,mpfi.interval);
        initialized = true;
    }
    
    // Move constructor.  Swaps the data between two inputs.
    MPFI::MPFI(MPFI &&mpfi)
    {
        mpfi_swap(interval,mpfi.interval);
        mpfi.initialized = false;
        initialized = true;
    }

    // Parameter constructors.  Sets upper and lower values to input value.
    // MPFI is not templated, so we need individual functions.
    // MPFI provides more functionality, but this should be enough for our purposes.
    
    // Using default precision.
    MPFI::MPFI(int value)
    {
        mpfi_init(interval);
        mpfi_set_si(interval,value);
        initialized = true;
    }

    MPFI::MPFI(unsigned int value)
    {
        mpfi_init(interval);
        mpfi_set_ui(interval,value);
        initialized = true;
    }

    MPFI::MPFI(float value)
    {
        mpfi_init(interval);
        mpfi_set_d(interval,value);
        initialized = true;
    }

    MPFI::MPFI(double value)
    {
        mpfi_init(interval);
        mpfi_set_d(interval,value);
        initialized = true;
    }
    
    MPFI::MPFI(const MPFR &mpfr)
    {
        mpfi_init2(interval,mpfr_get_prec(mpfr.point));
        mpfi_set_fr(interval,mpfr.point);
        initialized = true;
    }

    // Using user defined precision.
    MPFI::MPFI(int value, mpfr_prec_t prec)
    {
        mpfi_init2(interval, prec);
        mpfi_set_si(interval,value);
        initialized = true;
    }
    
    MPFI::MPFI(unsigned int value, mpfr_prec_t prec)
    {
        mpfi_init2(interval, prec);
        mpfi_set_ui(interval,value);
        initialized = true;
    }
    
    MPFI::MPFI(float value, mpfr_prec_t prec)
    {
        mpfi_init2(interval, prec);
        mpfi_set_d(interval,value);
        initialized = true;
    }
    
    MPFI::MPFI(double value, mpfr_prec_t prec)
    {
        mpfi_init2(interval, prec);
        mpfi_set_d(interval,value);
        initialized = true;
    }
 
    MPFI::MPFI(const MPFR &mpfr, mpfr_prec_t prec)
    {
        mpfi_init2(interval,prec);
        mpfi_set_fr(interval,mpfr.point);
        initialized = true;
    }

    // Assumes that the type can be cast to double.  Using default precision.
    template <typename Type>
    MPFI::MPFI(const Type &value)
    {
        double input = (double)value;
        mpfi_init(interval);
        mpfi_set_d(interval,input);
        initialized = true;
    }
    
    // Assumes that the type can be cast to double.  Using user defined precision.
    template <typename Type>
    MPFI::MPFI(const Type &value, mpfr_prec_t prec)
    {
        double input = (double)value;
        mpfi_init2(interval, prec);
        mpfi_set_d(interval,input);
        initialized = true;
    }
    
    // Parameter constructors.  Sets upper and lower values to input values.
    // MPFI is not templated, so we need individual functions.
    // MPFI provides more functionality, but this should be enough for our purposes.
    // The smaller of the two inputs is always the left (order does not matter).
    
    // Using default precision.
    MPFI::MPFI(int value1, int value2)
    {
        mpfi_init(interval);
        mpfi_interv_si(interval,value1,value2);
        initialized = true;
    }
    
    MPFI::MPFI(unsigned int value1, unsigned int value2)
    {
        mpfi_init(interval);
        mpfi_interv_ui(interval,value1,value2);
        initialized = true;
    }
    
    MPFI::MPFI(float value1, float value2)
    {
        mpfi_init(interval);
        mpfi_interv_d(interval,value1,value2);
        initialized = true;
    }
    
    MPFI::MPFI(double value1, double value2)
    {
        mpfi_init(interval);
        mpfi_interv_d(interval,value1,value2);
        initialized = true;
    }

    // Using maximum precision of the two endpoints.
    MPFI::MPFI(const MPFR & mpfr1,const MPFR & mpfr2)
    {
        mpfi_init2(interval,std::max(mpfr_get_prec(mpfr1.point),mpfr_get_prec(mpfr2.point)));
        mpfi_interv_fr(interval,mpfr1.point,mpfr2.point);
        initialized = true;
    }
    
    // Using user defined precision.
    MPFI::MPFI(int value1, int value2, mpfr_prec_t prec)
    {
        mpfi_init2(interval,prec);
        mpfi_interv_si(interval,value1,value2);
        initialized = true;
    }
    
    MPFI::MPFI(unsigned int value1, unsigned int value2, mpfr_prec_t prec)
    {
        mpfi_init2(interval,prec);
        mpfi_interv_ui(interval,value1,value2);
        initialized = true;
    }
    
    MPFI::MPFI(float value1, float value2, mpfr_prec_t prec)
    {
        mpfi_init2(interval,prec);
        mpfi_interv_d(interval,value1,value2);
        initialized = true;
    }
    
    MPFI::MPFI(double value1, double value2, mpfr_prec_t prec)
    {
        mpfi_init2(interval,prec);
        mpfi_interv_d(interval,value1,value2);
        initialized = true;
    }
    
    MPFI::MPFI(const MPFR & mpfr1,const MPFR & mpfr2, mpfr_prec_t prec)
    {
        mpfi_init2(interval,prec);
        mpfi_interv_fr(interval,mpfr1.point,mpfr2.point);
        initialized = true;
    }
    
    // Assumes that the type can be cast to double.  Using default precision.
    template <typename Type1, typename Type2>
    MPFI::MPFI(const Type1 & value1, const Type2 & value2)
    {
        double input1 = (double) value1, input2 = (double) value2;
        mpfi_init(interval);
        mpfi_interv_d(interval,input1,input2);
        initialized = true;
    }
    
    // Assumes that the type can be cast to double.  Using user defined precision.
    template <typename Type1, typename Type2>
    MPFI::MPFI(const Type1 & value1, const Type2 &value2, mpfr_prec_t prec)
    {
        double input1 = (double) value1, input2 = (double) value2;
        mpfi_init2(interval,prec);
        mpfi_interv_d(interval,input1,input2);
        initialized = true;
    }

    ///////////////////////
    // Destructor
    ///////////////////////
    
    // Destructor.  Frees memory in point.
    MPFI::~MPFI()
    {
        if(initialized)
             mpfi_clear(interval);
    }

    ///////////////////////
    // Operators
    ///////////////////////
    
    // Assignment operator.  Changes precision to match input precision.
    MPFI & MPFI::operator=(const MPFI &mpfi)
    {
        mpfi_set_prec(interval,mpfi_get_prec(mpfi.interval));
        mpfi_set(interval,mpfi.interval);
        return *this;
    }
    
    MPFI & MPFI::operator=(const MPFR &mpfr)
    {
        mpfi_set_prec(interval,mpfr_get_prec(mpfr.point));
        mpfi_set_fr(interval,mpfr.point);
        return *this;
    }
    
    // Assignment operator for Intervals.  Should be done via conversions.
    inline MPFI & MPFI::operator=(const Interval &interval)
    {
        convert(*this,interval);
        return *this;
    }
    
    // Move-assignment operator.  Swaps the data in the input.
    MPFI & MPFI::operator=(MPFI &&mpfi)
    {
        mpfi_swap(interval,mpfi.interval);
        bool backup = initialized;
        initialized = mpfi.initialized;
        mpfi.initialized = backup;
        return *this;
    }
    
    // Assignment operators.  Use current precision levels.
    // Using user defined precision.
    MPFI & MPFI::operator=(int value)
    {
        mpfi_set_si(interval,value);
        return *this;
    }
    
    MPFI & MPFI::operator=(unsigned int value)
    {
        mpfi_set_ui(interval,value);
        return *this;
    }
    
    MPFI & MPFI::operator=(float value)
    {
        mpfi_set_d(interval,value);
        return *this;
    }
    
    MPFI & MPFI::operator=(double value)
    {
        mpfi_set_d(interval,value);
        return *this;
    }
    
    // Assignment operator
    // Assumes that the type can be cast to double.  Using default precision.
    template <typename Type>
    MPFI & MPFI::operator=(const Type &value)
    {
        double input = (double)value;
        mpfi_set_d(interval,input);
        return *this;
    }
    
    // In-place addition operator.  Increases precision to maximum of input precisions.
    MPFI & MPFI::operator+=(const MPFI &mpfi)
    {
        if(mpfi_get_prec(interval)<mpfi_get_prec(mpfi.interval))
            mpfi_round_prec(interval,mpfi_get_prec(mpfi.interval));
        mpfi_add(interval,interval,mpfi.interval);
        return *this;
    }

    // In-place addition operator for Intervals.  Should be done via conversions.
    inline MPFI & MPFI::operator+=(const Interval &value)
    {
        mpfr_add_d(&(*interval).right,&(*interval).right,upper(value),MPFR_RNDU);
        mpfr_add_d(&(*interval).left,&(*interval).left,lower(value),MPFR_RNDD);
        return *this;
    }
    
    MPFI & MPFI::operator+=(const MPFR &mpfr)
    {
        if(mpfi_get_prec(interval)<mpfr_get_prec(mpfr.point))
            mpfi_round_prec(interval,mpfr_get_prec(mpfr.point));
        mpfi_add_fr(interval,interval,mpfr.point);
        return *this;
    }
    
    // In-place addition operators.  Using current precision.
    MPFI & MPFI::operator+=(int value)
    {
        mpfi_add_si(interval,interval,value);
        return *this;
    }
    
    MPFI & MPFI::operator+=(unsigned int value)
    {
        mpfi_add_ui(interval,interval,value);
        return *this;
    }
    
    MPFI & MPFI::operator+=(float value)
    {
        mpfi_add_d(interval,interval,value);
        return *this;
    }
    
    MPFI & MPFI::operator+=(double value)
    {
        mpfi_add_d(interval,interval,value);
        return *this;
    }
    
    // In-place addition operator.
    // Assumes that the type can be cast to double.  Using default precision.
    template <typename Type>
    MPFI & MPFI::operator+=(const Type &value)
    {
        double input = (double)value;
        mpfi_add_d(interval,interval,input);
        return *this;
    }

    // In-place subtraction operator.  Increases precision to maximum of input precisions.
    MPFI & MPFI::operator-=(const MPFI &mpfi)
    {
        if(mpfi_get_prec(interval)<mpfi_get_prec(mpfi.interval))
            mpfi_round_prec(interval,mpfi_get_prec(mpfi.interval));
        mpfi_sub(interval,interval,mpfi.interval);
        return *this;
    }
    
    MPFI & MPFI::operator-=(const MPFR &mpfr)
    {
        if(mpfi_get_prec(interval)<mpfr_get_prec(mpfr.point))
            mpfi_round_prec(interval,mpfr_get_prec(mpfr.point));
        mpfi_sub_fr(interval,interval,mpfr.point);
        return *this;
    }
    
    // In-place subtraction operator for Intervals.  Should be done via conversions.
    inline MPFI & MPFI::operator-=(const Interval &value)
    {
        mpfr_sub_d(&(*interval).right,&(*interval).right,lower(value),MPFR_RNDU);
        mpfr_sub_d(&(*interval).left,&(*interval).left,upper(value),MPFR_RNDD);
        return *this;
    }
    
    // In-place subtraction operators.  Using user defined precision.
    MPFI & MPFI::operator-=(int value)
    {
        mpfi_sub_si(interval,interval,value);
        return *this;
    }
    
    MPFI & MPFI::operator-=(unsigned int value)
    {
        mpfi_sub_ui(interval,interval,value);
        return *this;
    }
    
    MPFI & MPFI::operator-=(float value)
    {
        mpfi_sub_d(interval,interval,value);
        return *this;
    }
    
    MPFI & MPFI::operator-=(double value)
    {
        mpfi_sub_d(interval,interval,value);
        return *this;
    }
    
    // In-place subtraction operator.
    // Assumes that the type can be cast to double.  Using default precision.
    template <typename Type>
    MPFI & MPFI::operator-=(const Type &value)
    {
        double input = (double)value;
        mpfi_sub_d(interval,interval,input);
        return *this;
    }

    // In-place multiplication operator.  Increases precision to maximum of input precisions.
    MPFI & MPFI::operator*=(const MPFI &mpfi)
    {
        if(mpfi_get_prec(interval)<mpfi_get_prec(mpfi.interval))
            mpfi_round_prec(interval,mpfi_get_prec(mpfi.interval));
        //mpfi_mul(interval,interval,mpfi.interval);
        my_mpfi_mul(*this,*this,mpfi);
        return *this;
    }
    
    MPFI & MPFI::operator*=(const MPFR &mpfr)
    {
        if(mpfi_get_prec(interval)<mpfr_get_prec(mpfr.point))
            mpfi_round_prec(interval,mpfr_get_prec(mpfr.point));
        mpfi_mul_fr(interval,interval,mpfr.point);
        return *this;
    }
    
    // In-place multiplication operators.  Using user defined precision.
    MPFI & MPFI::operator*=(int value)
    {
        mpfi_mul_si(interval,interval,value);
        return *this;
    }
    
    MPFI & MPFI::operator*=(unsigned int value)
    {
        mpfi_mul_ui(interval,interval,value);
        return *this;
    }
    
    MPFI & MPFI::operator*=(float value)
    {
        mpfi_mul_d(interval,interval,value);
        return *this;
    }
    
    MPFI & MPFI::operator*=(double value)
    {
        mpfi_mul_d(interval,interval,value);
        return *this;
    }
    
    // In-place multiplication operator.
    // Assumes that the type can be cast to double.  Using default precision.
    template <typename Type>
    MPFI & MPFI::operator*=(const Type &value)
    {
        double input;
        convert(input,value);
        mpfi_mul_d(interval,interval,input);
        return *this;
    }

    // In-place division operator.  Increases precision to maximum of input precisions.
    MPFI & MPFI::operator/=(const MPFI &mpfi)
    {
        if(mpfi_get_prec(interval)<mpfi_get_prec(mpfi.interval))
            mpfi_round_prec(interval,mpfi_get_prec(mpfi.interval));
        mpfi_div(interval,interval,mpfi.interval);
        return *this;
    }
    
    MPFI & MPFI::operator/=(const MPFR &mpfr)
    {
        if(mpfi_get_prec(interval)<mpfr_get_prec(mpfr.point))
            mpfi_round_prec(interval,mpfr_get_prec(mpfr.point));
        mpfi_div_fr(interval,interval,mpfr.point);
        return *this;
    }
    
    // In-place division operators.  Using user defined precision.
    MPFI & MPFI::operator/=(int value)
    {
        mpfi_div_si(interval,interval,value);
        return *this;
    }
    
    MPFI & MPFI::operator/=(unsigned int value)
    {
        mpfi_div_ui(interval,interval,value);
        return *this;
    }
    
    MPFI & MPFI::operator/=(float value)
    {
        mpfi_div_d(interval,interval,value);
        return *this;
    }
    
    MPFI & MPFI::operator/=(double value)
    {
        mpfi_div_d(interval,interval,value);
        return *this;
    }
    
    // In-place division operator.
    // Assumes that the type can be cast to double.  Using default precision.
    template <typename Type>
    MPFI & MPFI::operator/=(const Type &value)
    {
        double input = (double)value;
        mpfi_div_d(interval,interval,input);
        return *this;
    }
    
    ///////////////////////
    // Interval data
    ///////////////////////
    
    // Returns or sets the width of an MPFI
    MPFR width(const MPFI &mpfi)
    {
        MPFR mpfr(mpfi_get_prec(mpfi.interval));
        mpfi_diam_abs(mpfr.point,mpfi.interval);
        return mpfr;
    }
    
    bool width(MPFR &mpfr, const MPFI &mpfi)
    {
        if(mpfi_get_prec(mpfi.interval)>mpfr_get_prec(mpfr.point))
            mpfr_prec_round(mpfr.point,mpfi_get_prec(mpfi.interval),MPFR_RNDN);
        return (mpfi_diam_abs(mpfr.point,mpfi.interval) == 0);
    }
    
    // Returns or assigns the upper endpoint of an MPFI
    MPFR upper(const MPFI &mpfi)
    {
        MPFR mpfr(mpfi_get_prec(mpfi.interval));
        mpfi_get_right(mpfr.point,mpfi.interval);
        return mpfr;
    }
    
    bool upper(MPFR &mpfr, const MPFI &mpfi)
    {
        if(mpfi_get_prec(mpfi.interval)>mpfr_get_prec(mpfr.point))
            mpfr_prec_round(mpfr.point,mpfi_get_prec(mpfi.interval),MPFR_RNDN);
        return (mpfi_get_right(mpfr.point,mpfi.interval) == 0);
    }
    
    // Returns or assigns the lower endpoint of an MPFI
    MPFR lower(const MPFI &mpfi)
    {
        MPFR mpfr(mpfi_get_prec(mpfi.interval));
        mpfi_get_left(mpfr.point,mpfi.interval);
        return mpfr;
    }
    
    bool lower(MPFR &mpfr, const MPFI &mpfi)
    {
        if(mpfi_get_prec(mpfi.interval)>mpfr_get_prec(mpfr.point))
            mpfr_prec_round(mpfr.point,mpfi_get_prec(mpfi.interval),MPFR_RNDN);
        return (mpfi_get_left(mpfr.point,mpfi.interval) == 0);
    }
    
    // Returns or sets the midpoint of an MPFI
    MPFR median(const MPFI &mpfi)
    {
        MPFR mpfr(mpfi_get_prec(mpfi.interval));
        mpfi_mid(mpfr.point,mpfi.interval);
        return mpfr;
    }
    
    bool median(MPFR &mpfr, const MPFI &mpfi)
    {
        if(mpfi_get_prec(mpfi.interval)>mpfr_get_prec(mpfr.point))
            mpfr_prec_round(mpfr.point,mpfi_get_prec(mpfi.interval),MPFR_RNDN);
        return (mpfi_mid(mpfr.point,mpfi.interval) == 0);
    }
    
    // Returns or sets largest absolute value of element of an MPFI
    MPFR abs(const MPFI &mpfi)
    {
        MPFR mpfr(mpfi_get_prec(mpfi.interval));
        mpfi_mag(mpfr.point,mpfi.interval);
        return mpfr;
    }
    
    bool abs(MPFR &mpfr, const MPFI &mpfi)
    {
        if(mpfi_get_prec(mpfi.interval)>mpfr_get_prec(mpfr.point))
            mpfr_prec_round(mpfr.point,mpfi_get_prec(mpfi.interval),MPFR_RNDN);
        return (mpfi_mag(mpfr.point,mpfi.interval) == 0);
    }

    // Gets the square of all the elements of an MPFI
    // The first function save the resulting interval to itself.
    // the second function save the resulting interval to its 1st argument
    inline MPFI square_norm(MPFI mpfi)
    {
        mpfi_sqr(mpfi.interval,mpfi.interval);
        return mpfi;
    }

    inline bool square_norm(MPFI& value, const MPFI & mpfi)
    {
        return mpfi_sqr(value.interval,mpfi.interval);
    }

    ///////////////////////
    // Output operator
    ///////////////////////
    
    // Display function for an MPFI.  Begins by converting MPFI to a double.
    // If higher output precision is needed, then we will need a new function.
    std::ostream& operator<<(std::ostream &out, const MPFI &mpfi)
    {
        if (mpfr_equal_p(&(*mpfi.interval).left,&(*mpfi.interval).right)!=0)
        {
            double value = mpfr_get_d(&(*mpfi.interval).left,MPFR_RNDN);
            return out << value;
        }
        double value1 = mpfr_get_d(&(*mpfi.interval).left,MPFR_RNDD),
               value2 = mpfr_get_d(&(*mpfi.interval).right,MPFR_RNDU);
        return out << "(" << value1 << "," << value2 << ")";
    }


    // Display function for an MPFI. Display the MPFI to a given length.
    inline bool print_more_digits(const MPFI &mpfi, unsigned int length)
    {
        MPFR value(mpfi_get_prec(mpfi.interval));
        mpfi_get_left(value.point, mpfi.interval);
        std::cout << "(";
        mpfr_out_str(stdout, 10, length, value.point, MPFR_RNDD);
        std::cout << ",";
        mpfi_get_right(value.point, mpfi.interval);
        mpfr_out_str(stdout, 10, length, value.point, MPFR_RNDU);
        std::cout << ")";
        return true;
    }
    ///////////////////////
    // Precision operators
    ///////////////////////
    
    // Returns the current precision.
    mpfr_prec_t get_precision(const MPFI & mpfi)
    {
        return mpfi_get_prec(mpfi.interval);
    }
    
    // Sets the precision of the MPFR (without destroying data).
    // Does not set precision below default 64.
    bool set_precision(MPFI & mpfi, mpfr_prec_t prec)
    {
        if ((prec<MPFR_PREC_MIN)||(prec>MPFR_PREC_MAX))
            return false;
        if (prec<64)
            prec = 64;
        mpfi_round_prec(mpfi.interval,prec);
        return true;
    }
    
    // Increases the precision to the input precision.
    // If there is only one input, the precision is doubled
    bool increase_precision(MPFI &mpfi)
    {
        return set_precision(mpfi,mpfi_get_prec(mpfi.interval)*2);
    }
    
    bool increase_precision(MPFI &mpfi, mpfr_prec_t prec)
    {
        if(prec>mpfi_get_prec(mpfi.interval))
            return set_precision(mpfi,prec);
        return false;
    }
    
    // Decreases the precision to the input precision.
    // If there is only one input, the precision is halved
    bool decrease_precision(MPFI &mpfi)
    {
        return set_precision(mpfi,mpfi_get_prec(mpfi.interval)/2);
    }
    
    bool decrease_precision(MPFI &mpfi, mpfr_prec_t prec)
    {
        if(prec<mpfi_get_prec(mpfi.interval))
            return set_precision(mpfi,prec);
        return false;
    }
    
    // Get exponent of MPFI (mantissa is assumed to be between 1/2 and 1)
    // Versions for the largest and smallest
    inline long max_exponent(MPFI &mpfi)
    {
        return (long)std::max(mpfr_get_exp(&(*mpfi.interval).left),
                              mpfr_get_exp(&(*mpfi.interval).right))-1;
    }
    
    inline long min_exponent(MPFI &mpfi)
    {
        return (long)std::min(mpfr_get_exp(&(*mpfi.interval).left),
                              mpfr_get_exp(&(*mpfi.interval).right))-1;
    }

    ///////////////////////
    // Operations
    ///////////////////////
    
    // Returns exponential of an MPFI.  The power must be an unsigned integer.
    // MPFI doesn't have a built-in power function, so we write it explicitly
    // The first input is set in the second version.
    MPFI exp(MPFI mpfi,unsigned int power)
    {
        MPFI base(1);
        while (power != 0)
        {
            if((power % 2)==1)
                base *= mpfi;
            power /= 2;
            mpfi *= mpfi;
        }
        return base;
    }
    
    bool exp(MPFI &mpfi,const MPFI &value, unsigned int power)
    {
        MPFI base = value;
        mpfi = 1;
        while (power != 0)
        {
            if((power % 2)==1)
                mpfi *= base;
            power /= 2;
            base *= base;
        }
        return true;
    }

    // Returns or sets the square root of an MPFI.  Rounds outward.
    // The first input is set in the second version.
    MPFI sqrt(MPFI mpfi)
    {
        mpfi_sqrt(mpfi.interval,mpfi.interval);
        return mpfi;
    }
    
    bool sqrt(MPFI &mpfi,const MPFI &value)
    {
        return (mpfi_sqrt(mpfi.interval,value.interval) == 0);
    }
 
    ///////////////////////
    // Friend Operators
    ///////////////////////

    // Addition operator.  Uses maximum precision of inputs
    MPFI operator+(MPFI mpfi1,const MPFI & mpfi2)
    {
        mpfi1 += mpfi2;
        return mpfi1;
    }
  
    MPFI operator+(MPFI mpfi,const MPFR &mpfr)
    {
        mpfi += mpfr;
        return mpfi;
    }
    
    
    MPFI operator+(const MPFR &mpfr, MPFI mpfi)
    {
        mpfi += mpfr;
        return mpfi;
    }
    
    // Addition operator for use with Intervals
    inline MPFI operator+(MPFI mpfi, const Interval &interval)
    {
        mpfi += interval;
        return mpfi;
    }
    
    inline MPFI operator+(const Interval &interval, MPFI mpfi)
    {
        mpfi += interval;
        return mpfi;
    }
  
    // Addition operator.  Using current precision
    MPFI operator+(MPFI mpfi, int value)
    {
        mpfi += value;
        return mpfi;
    }
    
    MPFI operator+(int value, MPFI mpfi)
    {
        mpfi += value;
        return mpfi;
    }

    MPFI operator+(MPFI mpfi,unsigned int value)
    {
        mpfi += value;
        return mpfi;
    }
    
    MPFI operator+(unsigned int value, MPFI mpfi)
    {
        mpfi += value;
        return mpfi;
    }
    
    MPFI operator+(MPFI mpfi, float value)
    {
        mpfi += value;
        return mpfi;
    }
    
    MPFI operator+(float value, MPFI mpfi)
    {
        mpfi += value;
        return mpfi;
    }
    
    MPFI operator+(MPFI mpfi, double value){
        mpfi += value;
        return mpfi;
    }
    
    MPFI operator+(double value, MPFI mpfi)
    {
        mpfi += value;
        return mpfi;
    }
    
    // Addition operator.  Using current precision
    // Assumes the input can be cast to a double.
    template <typename Type>
    MPFI operator+(MPFI mpfi, const Type &value)
    {
        mpfi += value;
        return mpfi;
    }
    
    template <typename Type>
    MPFI operator+(const Type &value, MPFI mpfi)
    {
        mpfi += value;
        return mpfi;
    }
 
    // Subtraction operator.  Uses maximum precision of inputs
    MPFI operator-(MPFI mpfi1,const MPFI &mpfi2)
    {
        mpfi1 -= mpfi2;
        return mpfi1;
    }
    
    MPFI operator-(MPFI mpfi,const MPFR &mpfr)
    {
        mpfi -= mpfr;
        return mpfi;
    }
    
    MPFI operator-(const MPFR &mpfr, MPFI mpfi)
    {
        mpfi_fr_sub(mpfi.interval,mpfr.point,mpfi.interval);
        return mpfi;
    }
    
    // Subtraction operator.  Using current precision
    MPFI operator-(MPFI mpfi, int value)
    {
        mpfi -= value;
        return mpfi;
    }
    
    MPFI operator-(int value, MPFI mpfi)
    {
        mpfi_si_sub(mpfi.interval,value,mpfi.interval);
        return mpfi;
    }

    MPFI operator-(MPFI mpfi,unsigned int value)
    {
        mpfi -= value;
        return mpfi;
    }
    
    MPFI operator-(unsigned int value, MPFI mpfi)
    {
        mpfi_ui_sub(mpfi.interval,value,mpfi.interval);
        return mpfi;
    }
    
    MPFI operator-(MPFI mpfi, float value)
    {
        mpfi -= value;
        return mpfi;
    }
    
    MPFI operator-(float value, MPFI mpfi)
    {
        mpfi_d_sub(mpfi.interval,value,mpfi.interval);
        return mpfi;
    }
    
    MPFI operator-(MPFI mpfi, double value)
    {
        mpfi -= value;
        return mpfi;
    }
    
    MPFI operator-(double value, MPFI mpfi)
    {
        mpfi_d_sub(mpfi.interval,value,mpfi.interval);
        return mpfi;
    }
    
    // Subtraction operator.  Using current precision
    // Assumes the input can be cast to a double.
    template <typename Type>
    MPFI operator-(MPFI mpfi, const Type &value)
    {
        double input = (double)value;
        mpfi -= input;
        return mpfi;
    }
    
    template <typename Type>
    MPFI operator-(const Type &value, MPFI mpfi)
    {
        double input = (double)value;
        mpfi_neg(mpfi.interval,mpfi.interval);
        mpfi += input;
        return mpfi;
    }

    // Multiplication operator.  Uses maximum precision of inputs
    MPFI operator*(MPFI mpfi1,const MPFI & mpfi2)
    {
        mpfi1 *= mpfi2;
        return mpfi1;
    }
    
    MPFI operator*(MPFI mpfi,const MPFR &mpfr)
    {
        mpfi *= mpfr;
        return mpfi;
    }
    
    MPFI operator*(const MPFR &mpfr, MPFI mpfi)
    {
        mpfi *= mpfr;
        return mpfi;
    }

    // Multiplication operator.  Using current precision
    MPFI operator*(MPFI mpfi, int value)
    {
        mpfi *= value;
        return mpfi;
    }
    
    MPFI operator*(int value, MPFI mpfi)
    {
        mpfi *= value;
        return mpfi;
    }
    
    MPFI operator*(MPFI mpfi,unsigned int value)
    {
        mpfi *= value;
        return mpfi;
    }
    
    MPFI operator*(unsigned int value, MPFI mpfi)
    {
        mpfi *= value;
        return mpfi;
    }
    
    MPFI operator*(MPFI mpfi, float value)
    {
        mpfi *= value;
        return mpfi;
    }
    
    MPFI operator*(float value, MPFI mpfi)
    {
        mpfi *= value;
        return mpfi;
    }
    
    MPFI operator*(MPFI mpfi, double value){
        mpfi *= value;
        return mpfi;
    }
    
    MPFI operator*(double value, MPFI mpfi)
    {
        mpfi *= value;
        return mpfi;
    }
    
    // Multiplication operator.  Using current precision
    // Assumes the input can be cast to a double.
    template <typename Type>
    MPFI operator*(MPFI mpfi, const Type &value)
    {
        mpfi *= value;
        return mpfi;
    }
    
    template <typename Type>
    MPFI operator*(const Type &value, MPFI mpfi)
    {
        mpfi *= value;
        return mpfi;
    }

    // Division operator.  Uses maximum precision of inputs
    MPFI operator/(MPFI mpfi1,const MPFI & mpfi2)
    {
        mpfi1 /= mpfi2;
        return mpfi1;
    }
  
    MPFI operator/(MPFI mpfi,const MPFR &mpfr)
    {
        mpfi /= mpfr;
        return mpfi;
    }
    
    MPFI operator/(const MPFR &mpfr, MPFI mpfi)
    {
        mpfi_fr_div(mpfi.interval,mpfr.point,mpfi.interval);
        return mpfi;
    }
  
    // Division operator.  Using current precision
    MPFI operator/(MPFI mpfi, int value)
    {
        mpfi /= value;
        return mpfi;
    }
    
    MPFI operator/(int value, MPFI mpfi)
    {
        mpfi_si_div(mpfi.interval,value,mpfi.interval);
        return mpfi;
    }
    
    MPFI operator/(MPFI mpfi,unsigned int value)
    {
        mpfi /= value;
        return mpfi;
    }
    
    MPFI operator/(unsigned int value, MPFI mpfi)
    {
        mpfi_ui_div(mpfi.interval,value,mpfi.interval);
        return mpfi;
    }
    
    MPFI operator/(MPFI mpfi, float value)
    {
        mpfi /= value;
        return mpfi;
    }
    
    MPFI operator/(float value, MPFI mpfi)
    {
        mpfi_d_div(mpfi.interval,value,mpfi.interval);
        return mpfi;
    }
    
    MPFI operator/(MPFI mpfi, double value){
        mpfi /= value;
        return mpfi;
    }
    
    MPFI operator/(double value, MPFI mpfi)
    {
        mpfi_d_div(mpfi.interval,value,mpfi.interval);
        return mpfi;
    }
    
    // Division operator.  Using current precision
    // Assumes the input can be cast to a double.
    template <typename Type>
    MPFI operator/(MPFI mpfi, const Type &value)
    {
        mpfi /= value;
        return mpfi;
    }
    
    template <typename Type>
    MPFI operator/(const Type &value, MPFI mpfi)
    {
        double input = (double)value;
        mpfi_d_div(mpfi.interval,input,mpfi.interval);
        return mpfi;
    }

    ///////////////////////
    // Comparison Operators
    ///////////////////////
    
    // Equality operator.
    bool operator==(const MPFI &mpfi1,const MPFI &mpfi2)
    {
        return (mpfr_equal_p(&(*mpfi1.interval).left,&(*mpfi2.interval).left)!=0)&&
                (mpfr_equal_p(&(*mpfi1.interval).right,&(*mpfi2.interval).right)!=0);
    }

    bool operator==(const MPFI &mpfi,const MPFR &mpfr)
    {
        return (mpfr_equal_p(&(*mpfi.interval).left,mpfr.point)!=0)&&
               (mpfr_equal_p(&(*mpfi.interval).right,mpfr.point)!=0);

    }
    
    bool operator==(const MPFR &mpfr, const MPFI &mpfi)
    {
        return (mpfr_equal_p(&(*mpfi.interval).left,mpfr.point)!=0)&&
               (mpfr_equal_p(&(*mpfi.interval).right,mpfr.point)!=0);
    }
    
    // Equality operator.  May have some problems if the MPFR is NaN.
    bool operator==(const MPFI & mpfi,int value)
    {
        return (mpfi.initialized)&&
               (mpfr_cmp_si(&(*mpfi.interval).left,value)==0)&&
               (mpfr_cmp_si(&(*mpfi.interval).right,value)==0);
    }
    
    bool operator==(int value,const MPFI & mpfi)
    {
        return (mpfi.initialized)&&
               (mpfr_cmp_si(&(*mpfi.interval).left,value)==0)&&
               (mpfr_cmp_si(&(*mpfi.interval).right,value)==0);
    }

    bool operator==(const MPFI & mpfi,unsigned int value)
    {
        return (mpfi.initialized)&&
               (mpfr_cmp_ui(&(*mpfi.interval).left,value)==0)&&
               (mpfr_cmp_ui(&(*mpfi.interval).right,value)==0);
    }
    
    bool operator==(unsigned int value, const MPFI &mpfi)
    {
        return (mpfi.initialized)&&
               (mpfr_cmp_ui(&(*mpfi.interval).left,value)==0)&&
               (mpfr_cmp_ui(&(*mpfi.interval).right,value)==0);
    }
    
    bool operator==(const MPFI &mpfi, float value)
    {
        return (mpfi.initialized)&&
               (mpfr_cmp_d(&(*mpfi.interval).left,value)==0)&&
               (mpfr_cmp_d(&(*mpfi.interval).right,value)==0);
    }
    
    bool operator==(float value, const MPFI & mpfi)
    {
        return (mpfi.initialized)&&
               (mpfr_cmp_d(&(*mpfi.interval).left,value)==0)&&
               (mpfr_cmp_d(&(*mpfi.interval).right,value)==0);
    }
    
    bool operator==(const MPFI & mpfi,double value)
    {
        return (mpfi.initialized)&&
               (mpfr_cmp_d(&(*mpfi.interval).left,value)==0)&&
               (mpfr_cmp_d(&(*mpfi.interval).right,value)==0);
    }
    
    bool operator==(double value, const MPFI & mpfi)
    {
        return (mpfi.initialized)&&
               (mpfr_cmp_d(&(*mpfi.interval).left,value)==0)&&
               (mpfr_cmp_d(&(*mpfi.interval).right,value)==0);
    }
    
    // Equality operator.
    // Assumes the input can be cast to a double.
    template <typename Type>
    bool operator==(const MPFI &mpfi,const Type & value)
    {
        double input = (double)value;
        return mpfi == input;
    }
    
    template <typename Type>
    bool operator==(const Type &value,const MPFI &mpfi)
    {
        double input = (double)value;
        return input == mpfi;
    }
    
    // Inquality operator.
    bool operator!=(const MPFI &mpfi1,const MPFI &mpfi2)
    {
        return (mpfi_nan_p(mpfi1.interval)==0)&&(mpfi_nan_p(mpfi2.interval)==0)&&!(mpfi1==mpfi2);
    }

    bool operator!=(const MPFI &mpfi,const MPFR &mpfr)
    {
        return (mpfi.initialized)&&!(mpfi==mpfr);
    }
    
    bool operator!=(const MPFR &mpfr,const MPFI &mpfi)
    {
        return (mpfi.initialized)&&!(mpfi==mpfr);
    }

    // Inequality operator.  May have some problems if the MPFR is NaN.
    bool operator!=(const MPFI & mpfi,int value)
    {
        return (mpfi.initialized)&&!(mpfi==value);
    }
    
    bool operator!=(int value,const MPFI & mpfi)
    {
        return (mpfi.initialized)&&!(mpfi==value);
    }

    bool operator!=(const MPFI & mpfi,unsigned int value)
    {
        return (mpfi.initialized)&&!(mpfi==value);
    }
    
    bool operator!=(unsigned int value, const MPFI &mpfi)
    {
        return (mpfi.initialized)&&!(mpfi==value);
    }
    
    bool operator!=(const MPFI &mpfi, float value)
    {
        return (mpfi.initialized)&&!(mpfi==value);
    }
    
    bool operator!=(float value, const MPFI & mpfi)
    {
        return (mpfi.initialized)&&!(mpfi==value);
    }
    
    bool operator!=(const MPFI & mpfi,double value)
    {
        return (mpfi.initialized)&&!(mpfi==value);
    }
    
    bool operator!=(double value, const MPFI & mpfi)
    {
        return (mpfi.initialized)&&!(mpfi==value);
    }
    
    // Inequality operator.
    // Assumes the input can be cast to a double.
    template <typename Type>
    bool operator!=(const MPFI &mpfi,const Type & value)
    {
        double input = (double) value;
        return input != mpfi;
    }
    
    template <typename Type>
    bool operator!=(const Type &value,const MPFI &mpfi)
    {
        double input = (double) value;
        return input != mpfi;
    }
    
    // Less than operator.
    bool operator<(const MPFI & mpfi1,const MPFI &mpfi2)
    {
        return mpfi_cmp(mpfi1.interval,mpfi2.interval)<0;
    }
 
    bool operator<(const MPFI &mpfi,const MPFR &mpfr)
    {
        return mpfi_cmp_fr(mpfi.interval,mpfr.point)<0;
    }
    
    bool operator<(const MPFR &mpfr,const MPFI &mpfi)
    {
        return (mpfi.initialized)&&(mpfr.initialized)&&(mpfi_cmp_fr(mpfi.interval,mpfr.point)>0);
    }
 
    // Less than operator.
    bool operator<(const MPFI &mpfi,int value)
    {
        return mpfi_cmp_si(mpfi.interval,value)<0;
    }
    
    bool operator<(int value,const MPFI &mpfi)
    {
        return (mpfi.initialized)&&(mpfi_cmp_si(mpfi.interval,value)>0);
    }
    
    bool operator<(const MPFI &mpfi,unsigned int value)
    {
        return mpfi_cmp_ui(mpfi.interval,value)<0;
    }
    
    bool operator<(unsigned int value, const MPFI &mpfi)
    {
        return (mpfi.initialized)&&(mpfi_cmp_ui(mpfi.interval,value)>0);
    }
    
    bool operator<(const MPFI &mpfi,float value)
    {
        return mpfi_cmp_d(mpfi.interval,value)<0;
    }
    
    bool operator<(float value, const MPFI &mpfi)
    {
        return (mpfi.initialized)&&(mpfi_cmp_d(mpfi.interval,value)>0);
    }
    
    bool operator<(const MPFI &mpfi,double value)
    {
        return mpfi_cmp_d(mpfi.interval,value)<0;
    }
    
    bool operator<(double value, const MPFI &mpfi)
    {
        return (mpfi.initialized)&&(mpfi_cmp_d(mpfi.interval,value)>0);
    }
    
    // Less than operator.
    // Assumes the input can be cast to a double.
    template <typename Type>
    bool operator<(const MPFI &mpfi,const Type &value)
    {
        double input = (double)value;
        return mpfi_cmp_d(mpfi.interval,input)<0;
    }
    
    template <typename Type>
    bool operator<(const Type &value,const MPFI &mpfi)
    {
        double input = (double)value;
        return (mpfi.initialized)&&(mpfi_cmp_d(mpfi.interval,input)>0);
    }

    // Less than or equal operator.
    bool operator<=(const MPFI & mpfi1,const MPFI &mpfi2)
    {
        return mpfr_lessequal_p(&(*mpfi1.interval).right,&(*mpfi2.interval).left)!=0;
    }
    
    bool operator<=(const MPFI &mpfi,const MPFR &mpfr)
    {
        return mpfr_lessequal_p(&(*mpfi.interval).right,mpfr.point)!=0;
    }
    
    bool operator<=(const MPFR &mpfr,const MPFI &mpfi)
    {
        return mpfr_greaterequal_p(&(*mpfi.interval).left,mpfr.point)!=0;
    }

    // Less than or equal operator.
    bool operator<=(const MPFI &mpfi,int value)
    {
        return (mpfi.initialized)&&(mpfr_cmp_si(&(*mpfi.interval).right,value)<=0);
    }
    
    bool operator<=(int value,const MPFI &mpfi)
    {
        return (mpfi.initialized)&&(mpfr_cmp_si(&(*mpfi.interval).left,value)>=0);
    }
    
    bool operator<=(const MPFI &mpfi,unsigned int value)
    {
        return (mpfi.initialized)&&(mpfr_cmp_ui(&(*mpfi.interval).right,value)<=0);
    }

    bool operator<=(unsigned int value, const MPFI &mpfi)
    {
        return (mpfi.initialized)&&(mpfr_cmp_ui(&(*mpfi.interval).left,value)>=0);
    }

    bool operator<=(const MPFI &mpfi,float value)
    {
        return (mpfi.initialized)&&(mpfr_cmp_d(&(*mpfi.interval).right,value)<=0);
    }

    bool operator<=(float value, const MPFI &mpfi)
    {
        return (mpfi.initialized)&&(mpfr_cmp_d(&(*mpfi.interval).left,value)>=0);
    }
    
    bool operator<=(const MPFI &mpfi,double value)
    {
        return (mpfi.initialized)&&(mpfr_cmp_d(&(*mpfi.interval).right,value)<=0);
    }

    bool operator<=(double value, const MPFI &mpfi)
    {
        return (mpfi.initialized)&&(mpfr_cmp_d(&(*mpfi.interval).left,value)>=0);
    }

    // Less than or equal operator.
    // Assumes the input can be cast to a double.
    template <typename Type>
    bool operator<=(const MPFI &mpfi,const Type &value)
    {
        double input = (double)value;
        return mpfi <= input;
    }
 
    template <typename Type>
    bool operator<=(const Type &value,const MPFI &mpfi)
    
    {
        double input = (double)value;
        return input <= mpfi;
    }
 
    // Greater than operator.
    bool operator>(const MPFI & mpfi1,const MPFI &mpfi2)
    {
        return mpfi2 < mpfi1;
    }
    
    bool operator>(const MPFI &mpfi,const MPFR &mpfr)
    {
        return mpfr < mpfi;
    }
    
    bool operator>(const MPFR &mpfr,const MPFI &mpfi)
    {
        return mpfi < mpfr;
    }
    
    // Greater than operator.
    bool operator>(const MPFI &mpfi,int value)
    {
        return value < mpfi;
    }
    
    bool operator>(int value,const MPFI &mpfi)
    {
        return mpfi < value;
    }
    
    bool operator>(const MPFI &mpfi,unsigned int value)
    {
        return value < mpfi;
    }
    
    bool operator>(unsigned int value, const MPFI &mpfi)
    {
        return mpfi < value;
    }
    
    bool operator>(const MPFI &mpfi,float value)
    {
        return value < mpfi;
    }
    
    bool operator>(float value, const MPFI &mpfi)
    {
        return mpfi < value;
    }
    
    bool operator>(const MPFI &mpfi,double value)
    {
        return value < mpfi;
    }
    
    bool operator>(double value, const MPFI &mpfi)
    {
        return mpfi < value;
    }
    
    // Greater than operator.
    // Assumes the input can be cast to a double.
    template <typename Type>
    bool operator>(const MPFI &mpfi,const Type &value)
    {
        double input = (double)value;
        return input < mpfi;
    }
    
    template <typename Type>
    bool operator>(const Type &value,const MPFI &mpfi)
    {
        double input = (double)value;
        return mpfi < input;
    }
    
    // Greater than or equal operator.
    bool operator>=(const MPFI & mpfi1,const MPFI &mpfi2)
    {
        return mpfi2 <= mpfi1;
    }
    
    bool operator>=(const MPFI &mpfi,const MPFR &mpfr)
    {
        return mpfr <= mpfi;
    }
    
    bool operator>=(const MPFR &mpfr,const MPFI &mpfi)
    {
        return mpfi <= mpfr;
    }
    
    // Greater than or equal operator.
    bool operator>=(const MPFI &mpfi,int value)
    {
        return value <= mpfi;
    }
    
    bool operator>=(int value,const MPFI &mpfi)
    {
        return mpfi <= value;
    }
    
    bool operator>=(const MPFI &mpfi,unsigned int value)
    {
        return value <= mpfi;
    }
    
    bool operator>=(unsigned int value, const MPFI &mpfi)
    {
        return mpfi <= value;
    }
    
    bool operator>=(const MPFI &mpfi,float value)
    {
        return value <= mpfi;
    }
    
    bool operator>=(float value, const MPFI &mpfi)
    {
        return mpfi <= value;
    }
    
    bool operator>=(const MPFI &mpfi,double value)
    {
        return value <= mpfi;
    }
    
    bool operator>=(double value, const MPFI &mpfi)
    {
        return mpfi <= value;
    }
    
    // Less than or equal operator.
    // Assumes the input can be cast to a double.
    template <typename Type>
    bool operator>=(const MPFI &mpfi,const Type &value)
    {
        double input = (double)value;
        return input <= mpfi;
    }
    
    template <typename Type>
    bool operator>=(const Type &value,const MPFI &mpfi)
    {
        double input = (double)value;
        return mpfi <= input;
    }
    
    ///////////////////////
    // Friend Inclusion Operators
    ///////////////////////

    // Returns true if the first value is contained in the interval
    bool contains(int value, const MPFI& mpfi)
    {
        return (mpfi.initialized) &&
               (mpfr_cmp_si(&(*mpfi.interval).left,value)<=0) &&
               (mpfr_cmp_si(&(*mpfi.interval).right,value)>=0);
    }
    
    bool contains(unsigned int value, const MPFI& mpfi)
    {
        return (mpfi.initialized) &&
               (mpfr_cmp_ui(&(*mpfi.interval).left,value)<=0) &&
               (mpfr_cmp_ui(&(*mpfi.interval).right,value)>=0);
    }
    
    bool contains(float value, const MPFI& mpfi)
    {
        return (mpfi.initialized) &&
               (mpfr_cmp_d(&(*mpfi.interval).left,value)<=0) &&
               (mpfr_cmp_d(&(*mpfi.interval).right,value)>=0);
    }
    
    bool contains(double value, const MPFI& mpfi)
    {
        return (mpfi.initialized) &&
               (mpfr_cmp_d(&(*mpfi.interval).left,value)<=0) &&
               (mpfr_cmp_d(&(*mpfi.interval).right,value)>=0);
    }
    
    bool contains(const MPFR& mpfr, const MPFI& mpfi)
    {
        return (mpfi.initialized) &&
               (mpfr_cmp(&(*mpfi.interval).left,mpfr.point)<=0) &&
               (mpfr_cmp(&(*mpfi.interval).right,mpfr.point)>=0);
    }
    
    // Templated version
    template <typename Type>
    bool contains(const Type& value, const MPFI& mpfi)
    {
        return (mpfi.initialized) &&
               (mpfr_cmp_d(&(*mpfi.interval).left,value)<=0) &&
               (mpfr_cmp_d(&(*mpfi.interval).right,value)>=0);
    }
    
    // Returns true if 0 is in the interval
    bool containsZero(const MPFI& value)
    {
        return mpfi_has_zero(value.interval)>0;
    }
    
    // Returns true if the interval has finite endpoints and
    // false if at least one endpoint is infinite
    bool is_finite(const MPFI &mpfi)
    {
        return (mpfi_bounded_p(mpfi.interval)!=0);
    }

}

#endif /* MPFI_tcc */
