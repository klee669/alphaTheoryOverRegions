//
//  MPFR.tcc
//  Homotopy Continuation
//
//  Created by burr2 on 8/8/18.
//  Copyright © 2018 burr2. All rights reserved.
//

#ifndef MPFR_tcc
#define MPFR_tcc

#include <iostream>
#include <cfenv>
#include <algorithm>
#include <cmath>
#include "MPFR.hpp"

namespace homotopy
{
    ///////////////////////
    // Constructors
    ///////////////////////

    // Default constructor.  Sets value to 0 using default precision.
    MPFR::MPFR()
    {
        mpfr_init(point);
        mpfr_set_si(point,0,MPFR_RNDN);
        initialized = true;
    }
    
    // Constructor with precision.  Sets value to 0 using input precision.
    MPFR::MPFR(mpfr_prec_t prec)
    {
        mpfr_init2(point,prec);
        mpfr_set_si(point,0,MPFR_RNDN);
        initialized = true;
    }
    
    // Copy constructor.  Copies value of one point to another.
    // Sets copy to have the correct precision level
    MPFR::MPFR(const MPFR &mpfr)
    {
        mpfr_init2(point,mpfr_get_prec(mpfr.point));
        mpfr_set(point,mpfr.point,MPFR_RNDN);
        initialized = true;
    }
    
    // Copy constructor with precision.  Copies value of one point to another.
    // Sets precision level equal to the user specified precision.
    MPFR::MPFR(const MPFR & mpfr, mpfr_prec_t prec)
    {
        mpfr_init2(point,prec);
        mpfr_set(point,mpfr.point,MPFR_RNDN);
        initialized = true;
    }
    
    // Move constructor.  Swaps the data between two inputs.
    MPFR::MPFR(MPFR &&mpfr)
    {
        mpfr_swap(point,mpfr.point);
        mpfr.initialized = false;
        initialized = true;
    }
    
    // Parameter constructors.
    // MPFR is not templated, so we need individual functions.
    // MPFR provides more functionality, but this should be enough for our purposes.
    
    // Using default precision.
    MPFR::MPFR(int value)
    {
        mpfr_init(point);
        mpfr_set_si(point,value,MPFR_RNDN);
        initialized = true;
    }
    
    MPFR::MPFR(unsigned int value)
    {
        mpfr_init(point);
        mpfr_set_ui(point,value,MPFR_RNDN);
        initialized = true;
    }
    
    MPFR::MPFR(float value)
    {
        mpfr_init(point);
        mpfr_set_flt(point,value,MPFR_RNDN);
        initialized = true;
    }
    
    MPFR::MPFR(double value)
    {
        mpfr_init(point);
        mpfr_set_d(point,value,MPFR_RNDN);
        initialized = true;
    }
    
    // Using user defined precision.
    MPFR::MPFR(int value, mpfr_prec_t prec)
    {
        mpfr_init2(point, prec);
        mpfr_set_si(point,value,MPFR_RNDN);
        initialized = true;
    }
    
    MPFR::MPFR(unsigned int value, mpfr_prec_t prec)
    {
        mpfr_init2(point, prec);
        mpfr_set_ui(point,value,MPFR_RNDN);
        initialized = true;
    }
    
    MPFR::MPFR(float value, mpfr_prec_t prec)
    {
        mpfr_init2(point, prec);
        mpfr_set_flt(point,value,MPFR_RNDN);
        initialized = true;
    }
    
    MPFR::MPFR(double value, mpfr_prec_t prec)
    {
        mpfr_init2(point, prec);
        mpfr_set_d(point,value,MPFR_RNDN);
        initialized = true;
    }

    MPFR::MPFR(const char *value, mpfr_prec_t prec)
    {
        mpfr_init2(point, prec);
        mpfr_init_set_str(point,value,10,MPFR_RNDA);
        initialized = true;
    }

    MPFR::MPFR(const char *value)
    {
        mpfr_init(point);
        mpfr_init_set_str(point,value,10,MPFR_RNDA);
        initialized = true;
    }


    MPFR::MPFR(std::string value, mpfr_prec_t prec)
    {
        mpfr_init2(point, prec);
        mpfr_init_set_str(point,value.c_str(),10,MPFR_RNDA);
        initialized = true;
    }
    MPFR::MPFR(std::string value)
    {
        mpfr_init(point);
        mpfr_init_set_str(point,value.c_str(),10,MPFR_RNDA);
        initialized = true;
    }



    // Assumes that the type can be cast to double.  Using default precision.
    template <typename Type>
    MPFR::MPFR(const Type &value)
    {
        double input = (double)value;
        mpfr_init(point);
        mpfr_set_d(point,input,MPFR_RNDN);
        initialized = true;
    }
    
    // Assumes that the type can be cast to double.  Using user defined precision.
    template <typename Type>
    MPFR::MPFR(const Type &value, mpfr_prec_t prec)
    {
        double input = (double)value;
        mpfr_init2(point, prec);
        mpfr_set_d(point,input,MPFR_RNDN);
        initialized = true;
    }
    
    ///////////////////////
    // Destructor
    ///////////////////////
 
    // Destructor.  Frees memory in point.
    MPFR::~MPFR()
    {
        if(initialized)
            mpfr_clear(point);
    }

    ///////////////////////
    // Operators
    ///////////////////////
    
    // Assignment operator.  Changes precision to match input precision.
    MPFR & MPFR::operator=(const MPFR &mpfr)
    {
#ifdef TEST_POLY
        std::cout << "In MPFR copy assignment, input precision " << mpfr_get_prec(mpfr.point) << std::endl;
#endif
        mpfr_set_prec(point,mpfr_get_prec(mpfr.point));
        mpfr_set(point,mpfr.point,MPFR_RNDN);
        return *this;
    }
    
    // Move-assignment operator.  Swaps the data in the input.
    MPFR & MPFR::operator=(MPFR &&mpfr)
    {
#ifdef TEST_POLY
        std::cout << "In MPFR copy assignment, input precision " << mpfr_get_prec(mpfr.point) << std::endl;
#endif
        mpfr_swap(point,mpfr.point);
        bool backup = initialized;
        initialized = mpfr.initialized;
        mpfr.initialized = backup;
        return *this;
    }
    
    // Assignment operators.  Use current precision levels.
    // Using current defined precision.
    MPFR & MPFR::operator=(int value)
    {
        mpfr_set_si(point,value,MPFR_RNDN);
        return *this;
    }
    
    MPFR & MPFR::operator=(unsigned int value)
    {
        mpfr_set_ui(point,value,MPFR_RNDN);
        return *this;
    }
    
    MPFR & MPFR::operator=(float value)
    {
        mpfr_set_flt(point,value,MPFR_RNDN);
        return *this;
    }
    
    MPFR & MPFR::operator=(double value)
    {
        mpfr_set_d(point,value,MPFR_RNDN);
        return *this;
    }

    // Assignment operators.  Using mpfr_t.
    MPFR & MPFR::operator=(const mpfr_t & value )
    {
        mpfr_set(point,value,MPFR_RNDN);
        return *this;
    } 

    // Assignment operator
    // Assumes that the type can be cast to double.  Using default precision.
    template <typename Type>
    MPFR & MPFR::operator=(const Type &value)
    {
        double input = (double)value;
        mpfr_set_d(point,input,MPFR_RNDN);
        return *this;
    }
    
    // In-place addition operator.  Increases precision to maximum of input precisions.
    MPFR & MPFR::operator+=(const MPFR &mpfr)
    {
        if(mpfr_get_prec(point)<mpfr_get_prec(mpfr.point))
            mpfr_prec_round(point,mpfr_get_prec(mpfr.point),MPFR_RNDN);
        mpfr_add(point,point,mpfr.point,MPFR_RNDN);
        return *this;
    }
    
    // In-place addition operators.  Using current defined precision.
    MPFR & MPFR::operator+=(int value)
    {
        mpfr_add_si(point,point,value,MPFR_RNDN);
        return *this;
    }
    
    MPFR & MPFR::operator+=(unsigned int value)
    {
        mpfr_add_ui(point,point,value,MPFR_RNDN);
        return *this;
    }
    
    MPFR & MPFR::operator+=(float value)
    {
        mpfr_add_d(point,point,value,MPFR_RNDN);
        return *this;
    }
    
    MPFR & MPFR::operator+=(double value)
    {
        mpfr_add_d(point,point,value,MPFR_RNDN);
        return *this;
    }
    
    // In-place addition operator.
    // Assumes that the type can be cast to double.  Using default precision.
    template <typename Type>
    MPFR & MPFR::operator+=(const Type &value)
    {
        double input = (double)value;
        mpfr_add_d(point,point,input,MPFR_RNDN);
        return *this;
    }
    
    // In-place subtraction operator.  Increases precision to maximum of input precisions.
    MPFR & MPFR::operator-=(const MPFR &mpfr)
    {
        if(mpfr_get_prec(point)<mpfr_get_prec(mpfr.point))
            mpfr_prec_round(point,mpfr_get_prec(mpfr.point),MPFR_RNDN);
        mpfr_sub(point,point,mpfr.point,MPFR_RNDN);
        return *this;
    }
    
    // In-place subtraction operators.  Using current defined precision.
    MPFR & MPFR::operator-=(int value)
    {
        mpfr_sub_si(point,point,value,MPFR_RNDN);
        return *this;
    }
    
    MPFR & MPFR::operator-=(unsigned int value)
    {
        mpfr_sub_ui(point,point,value,MPFR_RNDN);
        return *this;
    }
    
    MPFR & MPFR::operator-=(float value)
    {
        mpfr_sub_d(point,point,value,MPFR_RNDN);
        return *this;
    }
    
    MPFR & MPFR::operator-=(double value)
    {
        mpfr_sub_d(point,point,value,MPFR_RNDN);
        return *this;
    }
    
    // In-place subtraction operator.
    // Assumes that the type can be cast to double.  Using default precision.
    template <typename Type>
    MPFR & MPFR::operator-=(const Type &value)
    {
        double input = (double)value;
        mpfr_sub_d(point,point,input,MPFR_RNDN);
        return *this;
    }
    
    // In-place multiplication operator.  Increases precision to maximum of input precisions.
    MPFR & MPFR::operator*=(const MPFR &mpfr)
    {
        if(mpfr_get_prec(point)<mpfr_get_prec(mpfr.point))
            mpfr_prec_round(point,mpfr_get_prec(mpfr.point),MPFR_RNDN);
        mpfr_mul(point,point,mpfr.point,MPFR_RNDN);
        return *this;
    }
    
    // In-place multiplication operators.  Using current defined precision.
    MPFR & MPFR::operator*=(int value)
    {
        mpfr_mul_si(point,point,value,MPFR_RNDN);
        return *this;
    }
    
    MPFR & MPFR::operator*=(unsigned int value)
    {
        mpfr_mul_ui(point,point,value,MPFR_RNDN);
        return *this;
    }
    
    MPFR & MPFR::operator*=(float value)
    {
        mpfr_mul_d(point,point,value,MPFR_RNDN);
        return *this;
    }
    
    MPFR & MPFR::operator*=(double value)
    {
        mpfr_mul_d(point,point,value,MPFR_RNDN);
        return *this;
    }
    
    // In-place multiplication operator.
    // Assumes that the type can be cast to double.  Using default precision.
    template <typename Type>
    MPFR & MPFR::operator*=(const Type &value)
    {
        double input;
        convert(input,value);
        mpfr_mul_d(point,point,input,MPFR_RNDN);
        return *this;
    }

    // In-place division operator.  Increases precision to maximum of input precisions.
    MPFR & MPFR::operator/=(const MPFR &mpfr)
    {
        if(mpfr_get_prec(point)<mpfr_get_prec(mpfr.point))
            mpfr_prec_round(point,mpfr_get_prec(mpfr.point),MPFR_RNDN);
        mpfr_div(point,point,mpfr.point,MPFR_RNDN);
        return *this;
    }
    
    // In-place division operators.  Using current defined precision.
    MPFR & MPFR::operator/=(int value)
    {
        mpfr_div_si(point,point,value,MPFR_RNDN);
        return *this;
    }
    
    MPFR & MPFR::operator/=(unsigned int value)
    {
        mpfr_div_ui(point,point,value,MPFR_RNDN);
        return *this;
    }
    
    MPFR & MPFR::operator/=(float value)
    {
        mpfr_div_d(point,point,value,MPFR_RNDN);
        return *this;
    }
    
    MPFR & MPFR::operator/=(double value)
    {
        mpfr_div_d(point,point,value,MPFR_RNDN);
        return *this;
    }
    
    // In-place division operator.
    // Assumes that the type can be cast to double.  Using default precision.
    template <typename Type>
    MPFR & MPFR::operator/=(const Type &value)
    {
        double input = (double)value;
        mpfr_div_d(point,point,input,MPFR_RNDN);
        return *this;
    }

    ///////////////////////
    // Output operator
    ///////////////////////
 
    // Display function for an MPFR.  Begins by converting MPFR to a double.
    // If higher output precision is needed, then we will need a new function.
    std::ostream& operator<<(std::ostream &out, const MPFR &mpfr)
    {
        double value = mpfr_get_d(mpfr.point,MPFR_RNDN);
        return out << value;
    }

    // Display an MPFR to a given length.
    inline bool print_more_digits(const MPFR & mpfr, unsigned int length)
    {
        mpfr_out_str(stdout, 10, length, mpfr.point, MPFR_RNDN) ;
        return true;
    }
 
    ///////////////////////
    // Precision operators
    ///////////////////////
    
    // Returns the current precision.
    mpfr_prec_t get_precision(const MPFR & mpfr)
    {
        return mpfr_get_prec(mpfr.point);
    }
    
    // Sets the precision of the MPFR (without destroying data).
    // Does not set precision below default 64.
    bool set_precision(MPFR & mpfr, mpfr_prec_t prec)
    {
        if ((prec<MPFR_PREC_MIN)||(prec>MPFR_PREC_MAX))
            return false;
        if (prec<64)
            prec = 64;
        mpfr_prec_round(mpfr.point,prec,MPFR_RNDN);
        return true;
    }
    
    // Increases the precision to the input precision.
    // If there is only one input, the precision is doubled
    bool increase_precision(MPFR &mpfr)
    {
        return set_precision(mpfr,mpfr_get_prec(mpfr.point)*2);
    }
    
    bool increase_precision(MPFR &mpfr, mpfr_prec_t prec)
    {
        if(prec>mpfr_get_prec(mpfr.point))
            return set_precision(mpfr,prec);
        return false;
    }
    
    // Decreases the precision to the input precision.
    // If there is only one input, the precision is halved
    bool decrease_precision(MPFR &mpfr)
    {
        return set_precision(mpfr,mpfr_get_prec(mpfr.point)/2);
    }
    
    bool decrease_precision(MPFR &mpfr, mpfr_prec_t prec)
    {
        if(prec<mpfr_get_prec(mpfr.point))
            return set_precision(mpfr,prec);
        return false;
    }
    
    // Get exponent of MPFR (mantissa is assumed to be between 1/2 and 1)
    long exponent(MPFR &mpfr)
    {
        return (long)mpfr_get_exp(mpfr.point)-1;
    }





    ///////////////////////
    // Operations
    ///////////////////////
    
    // Return or sets the square of the mpfr
    MPFR square_norm(MPFR mpfr)
    {
        mpfr_sqr(mpfr.point, mpfr.point, MPFR_RNDN);
        return mpfr;
    }
    bool square_norm(MPFR &value, const MPFR& mpfr)
    {
        return mpfr_sqr(value.point, mpfr.point, MPFR_RNDN);
    }
    
    // Returns or sets the absolute value of the mpfr
    MPFR abs(MPFR mpfr)
    {
        mpfr_abs(mpfr.point,mpfr.point,MPFR_RNDN);
        return mpfr;
    }
    
    bool abs(MPFR &mpfr,const MPFR &value)
    {
        return (mpfr_abs(mpfr.point,value.point,MPFR_RNDN) == 0);
    }
    
    // Returns or assigns exponential of an MPFR.  The power must be an unsigned integer.
    MPFR exp(MPFR mpfr,unsigned int power)
    {
        mpfr_pow_ui(mpfr.point,mpfr.point,power,MPFR_RNDN);
        return mpfr;
    }
    
    bool exp(MPFR &mpfr, const MPFR &value, unsigned int power)
    {
        return (mpfr_pow_ui(mpfr.point,value.point,power,MPFR_RNDN) == 0);
    }
    
    // Returns or assigns the square root of an MPFR.  Rounds to nearest.
    MPFR sqrt(MPFR mpfr)
    {
        mpfr_sqrt(mpfr.point,mpfr.point,MPFR_RNDN);
        return mpfr;
    }
    
    bool sqrt(MPFR &mpfr,const MPFR &value)
    {
        return (mpfr_sqrt(mpfr.point,value.point,MPFR_RNDN) == 0);
    }

    
    ///////////////////////
    // Operators
    ///////////////////////
    
    // Addition operator.  Uses maximum precision of inputs
    MPFR operator+(MPFR mpfr1,const MPFR & mpfr2)
    {
        mpfr1 += mpfr2;
        return mpfr1;
    }
    
    // Addition operator.  Using current precision
    MPFR operator+(MPFR mpfr, int value)
    {
        mpfr += value;
        return mpfr;
    }
    
    MPFR operator+(int value, MPFR mpfr)
    {
        mpfr += value;
        return mpfr;
    }
    
    MPFR operator+(MPFR mpfr,unsigned int value)
    {
        mpfr += value;
        return mpfr;
    }
    
    MPFR operator+(unsigned int value, MPFR mpfr)
    {
        mpfr += value;
        return mpfr;
    }
    
    MPFR operator+(MPFR mpfr, float value)
    {
        mpfr += value;
        return mpfr;
    }
    
    MPFR operator+(float value, MPFR mpfr)
    {
        mpfr += value;
        return mpfr;
    }
    
    MPFR operator+(MPFR mpfr, double value){
        mpfr += value;
        return mpfr;
    }
    
    MPFR operator+(double value, MPFR mpfr)
    {
        mpfr += value;
        return mpfr;
    }
    
    // Addition operator.  Using current precision
    // Assumes the input can be cast to a double.
    template <typename Type>
    MPFR operator+(MPFR mpfr, const Type &value)
    {
        mpfr += value;
        return mpfr;
    }
    
    template <typename Type>
    MPFR operator+(const Type &value, MPFR mpfr)
    {
        mpfr += value;
        return mpfr;
    }
    
    // Subtraction operator.  Uses maximum precision of inputs
    MPFR operator-(MPFR mpfr1,const MPFR &mpfr2)
    {
        mpfr1 -= mpfr2;
        return mpfr1;
    }
    
    // Subtraction operator.  Using current precision
    MPFR operator-(MPFR mpfr, int value)
    {
        mpfr -= value;
        return mpfr;
    }
    
    MPFR operator-(int value, MPFR mpfr)
    {
        mpfr_si_sub(mpfr.point,value,mpfr.point,MPFR_RNDN);
        return mpfr;
    }
    
    MPFR operator-(MPFR mpfr,unsigned int value)
    {
        mpfr -= value;
        return mpfr;
    }
    
    MPFR operator-(unsigned int value, MPFR mpfr)
    {
        mpfr_ui_sub(mpfr.point,value,mpfr.point,MPFR_RNDN);
        return mpfr;
    }
    
    MPFR operator-(MPFR mpfr, float value)
    {
        mpfr -= value;
        return mpfr;
    }
    
    MPFR operator-(float value, MPFR mpfr)
    {
        mpfr_d_sub(mpfr.point,value,mpfr.point,MPFR_RNDN);
        return mpfr;
    }
    
    MPFR operator-(MPFR mpfr, double value)
    {
        mpfr -= value;
        return mpfr;
    }
    
    MPFR operator-(double value, MPFR mpfr)
    {
        mpfr_d_sub(mpfr.point,value,mpfr.point,MPFR_RNDN);
        return mpfr;
    }
    
    // Subtraction operator.  Using current precision
    // Assumes the input can be cast to a double.
    template <typename Type>
    MPFR operator-(MPFR mpfr, const Type &value)
    {
        double input = (double)value;
        mpfr -= input;
        return mpfr;
    }
    
    template <typename Type>
    MPFR operator-(const Type &value, MPFR mpfr)
    {
        double input = (double)value;
        mpfr_neg(mpfr.point,mpfr.point,MPFR_RNDN);
        mpfr += input;
        return mpfr;
    }
    
    // Multiplication operator.  Uses maximum precision of inputs
    MPFR operator*(MPFR mpfr1,const MPFR & mpfr2)
    {
        mpfr1 *= mpfr2;
        return mpfr1;
    }
    
    // Multiplication operator.  Using current precision
    MPFR operator*(MPFR mpfr, int value)
    {
        mpfr *= value;
        return mpfr;
    }
    
    MPFR operator*(int value, MPFR mpfr)
    {
        mpfr *= value;
        return mpfr;
    }
    
    MPFR operator*(MPFR mpfr,unsigned int value)
    {
        mpfr *= value;
        return mpfr;
    }
    
    MPFR operator*(unsigned int value, MPFR mpfr)
    {
        mpfr *= value;
        return mpfr;
    }
    
    MPFR operator*(MPFR mpfr, float value)
    {
        mpfr *= value;
        return mpfr;
    }
    
    MPFR operator*(float value, MPFR mpfr)
    {
        mpfr *= value;
        return mpfr;
    }
    
    MPFR operator*(MPFR mpfr, double value){
        mpfr *= value;
        return mpfr;
    }
    
    MPFR operator*(double value, MPFR mpfr)
    {
        mpfr *= value;
        return mpfr;
    }
    
    // Multiplication operator.  Using current precision
    // Assumes the input can be cast to a double.
    template <typename Type>
    MPFR operator*(MPFR mpfr, const Type &value)
    {
        mpfr *= value;
        return mpfr;
    }
    
    template <typename Type>
    MPFR operator*(const Type &value, MPFR mpfr)
    {
        mpfr *= value;
        return mpfr;
    }
    
    // Division operator.  Uses maximum precision of inputs
    MPFR operator/(MPFR mpfr1,const MPFR & mpfr2)
    {
        mpfr1 /= mpfr2;
        return mpfr1;
    }
    
    // Division operator.  Using current precision
    MPFR operator/(MPFR mpfr, int value)
    {
        mpfr /= value;
        return mpfr;
    }
    
    MPFR operator/(int value, MPFR mpfr)
    {
        mpfr_si_div(mpfr.point,value,mpfr.point,MPFR_RNDN);
        return mpfr;
    }
    
    MPFR operator/(MPFR mpfr,unsigned int value)
    {
        mpfr /= value;
        return mpfr;
    }
    
    MPFR operator/(unsigned int value, MPFR mpfr)
    {
        mpfr_ui_div(mpfr.point,value,mpfr.point,MPFR_RNDN);
        return mpfr;
    }
    
    MPFR operator/(MPFR mpfr, float value)
    {
        mpfr /= value;
        return mpfr;
    }
    
    MPFR operator/(float value, MPFR mpfr)
    {
        mpfr_d_div(mpfr.point,value,mpfr.point,MPFR_RNDN);
        return mpfr;
    }
    
    MPFR operator/(MPFR mpfr, double value){
        mpfr /= value;
        return mpfr;
    }
    
    MPFR operator/(double value, MPFR mpfr)
    {
        mpfr_d_div(mpfr.point,value,mpfr.point,MPFR_RNDN);
        return mpfr;
    }
    
    // Division operator.  Using current precision
    // Assumes the input can be cast to a double.
    template <typename Type>
    MPFR operator/(MPFR mpfr, const Type &value)
    {
        mpfr /= value;
        return mpfr;
    }
    
    template <typename Type>
    MPFR operator/(const Type &value, MPFR mpfr)
    {
        double input = (double)value;
        mpfr_d_div(mpfr.point,input,mpfr.point,MPFR_RNDN);
        return mpfr;
    }

    ///////////////////////
    // Comparison Operators
    ///////////////////////
    
    // Equality operator.
    bool operator==(const MPFR &mpfr1,const MPFR &mpfr2)
    {
        return (mpfr_equal_p(mpfr1.point,mpfr2.point)!=0);
    }
    
    // Equality operator.  May have some problems if the MPFR is NaN.
    bool operator==(const MPFR & mpfr,int value)
    {
        return (mpfr.initialized)&&(mpfr_cmp_si(mpfr.point,value)==0);
    }
    
    bool operator==(int value,const MPFR & mpfr)
    {
        return (mpfr.initialized)&&(mpfr_cmp_si(mpfr.point,value)==0);
    }
    
    bool operator==(const MPFR & mpfr,unsigned int value)
    {
        return (mpfr.initialized)&&(mpfr_cmp_ui(mpfr.point,value)==0);
    }
    
    bool operator==(unsigned int value, const MPFR &mpfr)
    {
        return (mpfr.initialized)&&(mpfr_cmp_ui(mpfr.point,value)==0);
    }
    
    bool operator==(const MPFR &mpfr, float value)
    {
        return (mpfr.initialized)&&(mpfr_cmp_d(mpfr.point,value)==0);
    }
    
    bool operator==(float value, const MPFR & mpfr)
    {
        return (mpfr.initialized)&&(mpfr_cmp_d(mpfr.point,value)==0);
    }
    
    bool operator==(const MPFR & mpfr,double value)
    {
        return (mpfr.initialized)&&(mpfr_cmp_d(mpfr.point,value)==0);
    }
    
    bool operator==(double value, const MPFR & mpfr)
    {
        return (mpfr.initialized)&&(mpfr_cmp_d(mpfr.point,value)==0);
    }
    
    // Equality operator.
    // Assumes the input can be cast to a double.
    template <typename Type>
    bool operator==(const MPFR &mpfr,const Type & value)
    {
        double input = (double)value;
        return (mpfr.initialized)&&(mpfr_cmp_d(mpfr.point,input)==0);
    }
    
    template <typename Type>
    bool operator==(const Type &value,const MPFR &mpfr)
    {
        double input = (double)value;
        return (mpfr.initialized)&&(mpfr_cmp_d(mpfr.point,input)==0);
    }
    
    // Inquality operator.
    bool operator!=(const MPFR &mpfr1,const MPFR &mpfr2)
    {
        return (mpfr_equal_p(mpfr1.point,mpfr2.point)==0);
    }
    
    // Inequality operator.  May have some problems if the MPFR is NaN.
    bool operator!=(const MPFR & mpfr,int value)
    {
        return ((mpfr.initialized)&&!(mpfr==value));
    }
    
    bool operator!=(int value,const MPFR & mpfr)
    {
        return ((mpfr.initialized)&&!(mpfr==value));
    }
    
    bool operator!=(const MPFR & mpfr,unsigned int value)
    {
        return ((mpfr.initialized)&&!(mpfr==value));
    }
    
    bool operator!=(unsigned int value, const MPFR &mpfr)
    {
        return ((mpfr.initialized)&&!(mpfr==value));
    }
    
    bool operator!=(const MPFR &mpfr, float value)
    {
        return ((mpfr.initialized)&&!(mpfr==value));
    }
    
    bool operator!=(float value, const MPFR & mpfr)
    {
        return ((mpfr.initialized)&&!(mpfr==value));
    }
    
    bool operator!=(const MPFR & mpfr,double value)
    {
        return ((mpfr.initialized)&&!(mpfr==value));
    }
    
    bool operator!=(double value, const MPFR & mpfr)
    {
        return ((mpfr.initialized)&&!(mpfr==value));
    }
    
    // Inequality operator.
    // Assumes the input can be cast to a double.
    template <typename Type>
    bool operator!=(const MPFR &mpfr,const Type & value)
    {
        return ((mpfr.initialized)&&!(mpfr==value));
    }
    
    template <typename Type>
    bool operator!=(const Type &value,const MPFR &mpfr)
    {
        return ((mpfr.initialized)&&!(mpfr==value));
    }
    
    // Less than operator.
    bool operator<(const MPFR & mpfr1,const MPFR &mpfr2)
    {
        return mpfr_less_p(mpfr1.point,mpfr2.point);
    }
    
    // Less than operator.
    bool operator<(const MPFR &mpfr,int value)
    {
        return (mpfr_cmp_si(mpfr.point,value)<0);
    }
    
    bool operator<(int value,const MPFR &mpfr)
    {
        return (mpfr_cmp_si(mpfr.point,value)>0);
    }
    
    bool operator<(const MPFR &mpfr,unsigned int value)
    {
        return (mpfr_cmp_ui(mpfr.point,value)<0);
    }
    
    bool operator<(unsigned int value, const MPFR &mpfr)
    {
        return (mpfr_cmp_ui(mpfr.point,value)>0);
    }
    
    bool operator<(const MPFR &mpfr,float value)
    {
        return (mpfr_cmp_d(mpfr.point,value)<0);
    }
    
    bool operator<(float value, const MPFR &mpfr)
    {
        return (mpfr_cmp_d(mpfr.point,value)>0);
    }
    
    bool operator<(const MPFR &mpfr,double value)
    {
        return (mpfr_cmp_d(mpfr.point,value)<0);
    }
    
    bool operator<(double value, const MPFR &mpfr)
    {
        return (mpfr_cmp_d(mpfr.point,value)>0);
    }
    
    // Less than operator.
    // Assumes the input can be cast to a double.
    template <typename Type>
    bool operator<(const MPFR &mpfr,const Type &value)
    {
        double input = (double)value;
        return (mpfr_cmp_d(mpfr.point,input)<0);
    }
    
    template <typename Type>
    bool operator<(const Type &value,const MPFR &mpfr)
    {
        double input = (double)value;
        return (mpfr_cmp_d(mpfr.point,input)>0);
    }
    
    // Less than or equal operator.
    bool operator<=(const MPFR & mpfr1,const MPFR &mpfr2)
    {
        return mpfr_lessequal_p(mpfr1.point,mpfr2.point);
    }
    
    // Less than or equal operator.
    bool operator<=(const MPFR &mpfr,int value)
    {
        return (mpfr.initialized)&&(mpfr_cmp_si(mpfr.point,value)<=0);
    }
    
    bool operator<=(int value,const MPFR &mpfr)
    {
        return (mpfr.initialized)&&(mpfr_cmp_si(mpfr.point,value)>=0);
    }
    
    bool operator<=(const MPFR &mpfr,unsigned int value)
    {
        return (mpfr.initialized)&&(mpfr_cmp_ui(mpfr.point,value)<=0);
    }
    
    bool operator<=(unsigned int value, const MPFR &mpfr)
    {
        return (mpfr.initialized)&&(mpfr_cmp_ui(mpfr.point,value)>=0);
    }
    
    bool operator<=(const MPFR &mpfr,float value)
    {
        return (mpfr.initialized)&&(mpfr_cmp_d(mpfr.point,value)<=0);
    }
    
    bool operator<=(float value, const MPFR &mpfr)
    {
        return (mpfr.initialized)&&(mpfr_cmp_d(mpfr.point,value)>=0);
    }
    
    bool operator<=(const MPFR &mpfr,double value)
    {
        return (mpfr.initialized)&&(mpfr_cmp_d(mpfr.point,value)<=0);
    }
    
    bool operator<=(double value, const MPFR &mpfr)
    {
        return (mpfr.initialized)&&(mpfr_cmp_d(mpfr.point,value)>=0);
    }
    
    // Less than or equal operator.
    // Assumes the input can be cast to a double.
    template <typename Type>
    bool operator<=(const MPFR &mpfr,const Type &value)
    {
        double input = (double)value;
        return (mpfr.initialized)&&(mpfr_cmp_d(mpfr.point,input)<=0);
    }
    
    template <typename Type>
    bool operator<=(const Type &value,const MPFR &mpfr)
    {
        double input = (double)value;
        return (mpfr.initialized)&&(mpfr_cmp_d(mpfr.point,input)>=0);
    }
    
// *****************
    
    // Greater than operator.
    bool operator>(const MPFR & mpfr1,const MPFR &mpfr2)
    {
        return mpfr_greater_p(mpfr1.point,mpfr2.point);
    }
    
    // Greater than operator.
    bool operator>(const MPFR &mpfr,int value)
    {
        return (mpfr_cmp_si(mpfr.point,value)>0);
    }
    
    bool operator>(int value,const MPFR &mpfr)
    {
        return (mpfr_cmp_si(mpfr.point,value)<0);
    }
    
    bool operator>(const MPFR &mpfr,unsigned int value)
    {
        return (mpfr_cmp_ui(mpfr.point,value)>0);
    }
    
    bool operator>(unsigned int value, const MPFR &mpfr)
    {
        return (mpfr_cmp_ui(mpfr.point,value)<0);
    }
    
    bool operator>(const MPFR &mpfr,float value)
    {
        return (mpfr_cmp_d(mpfr.point,value)>0);
    }
    
    bool operator>(float value, const MPFR &mpfr)
    {
        return (mpfr_cmp_d(mpfr.point,value)<0);
    }
    
    bool operator>(const MPFR &mpfr,double value)
    {
        return (mpfr_cmp_d(mpfr.point,value)>0);
    }
    
    bool operator>(double value, const MPFR &mpfr)
    {
        return (mpfr_cmp_d(mpfr.point,value)<0);
    }
    
    // Greater than operator.
    // Assumes the input can be cast to a double.
    template <typename Type>
    bool operator>(const MPFR &mpfr,const Type &value)
    {
        double input = (double)value;
        return (mpfr_cmp_d(mpfr.point,input)>0);
    }
    
    template <typename Type>
    bool operator>(const Type &value,const MPFR &mpfr)
    {
        double input = (double)value;
        return (mpfr_cmp_d(mpfr.point,input)<0);
    }
    
    // Less than or equal operator.
    bool operator>=(const MPFR & mpfr1,const MPFR &mpfr2)
    {
        return mpfr_greaterequal_p(mpfr1.point,mpfr2.point);
    }
    
    // Less than or equal operator.
    bool operator>=(const MPFR &mpfr,int value)
    {
        return (mpfr.initialized)&&(mpfr_cmp_si(mpfr.point,value)>=0);
    }
    
    bool operator>=(int value,const MPFR &mpfr)
    {
        return (mpfr.initialized)&&(mpfr_cmp_si(mpfr.point,value)<=0);
    }
    
    bool operator>=(const MPFR &mpfr,unsigned int value)
    {
        return (mpfr.initialized)&&(mpfr_cmp_ui(mpfr.point,value)>=0);
    }
    
    bool operator>=(unsigned int value, const MPFR &mpfr)
    {
        return (mpfr.initialized)&&(mpfr_cmp_ui(mpfr.point,value)<=0);
    }
    
    bool operator>=(const MPFR &mpfr,float value)
    {
        return (mpfr.initialized)&&(mpfr_cmp_d(mpfr.point,value)>=0);
    }
    
    bool operator>=(float value, const MPFR &mpfr)
    {
        return (mpfr.initialized)&&(mpfr_cmp_d(mpfr.point,value)<=0);
    }
    
    bool operator>=(const MPFR &mpfr,double value)
    {
        return (mpfr.initialized)&&(mpfr_cmp_d(mpfr.point,value)>=0);
    }
    
    bool operator>=(double value, const MPFR &mpfr)
    {
        return (mpfr.initialized)&&(mpfr_cmp_d(mpfr.point,value)<=0);
    }
    
    // Greater than or equal operator.
    // Assumes the input can be cast to a double.
    template <typename Type>
    bool operator>=(const MPFR &mpfr,const Type &value)
    {
        double input = (double)value;
        return (mpfr.initialized)&&(mpfr_cmp_d(mpfr.point,input)>=0);
    }
    
    template <typename Type>
    bool operator>=(const Type &value,const MPFR &mpfr)
    {
        double input = (double)value;
        return (mpfr.initialized)&&(mpfr_cmp_d(mpfr.point,input)<=0);
    }

    inline MPFR floor(MPFR mpfr)
    {
        mpfr_floor(mpfr.point, mpfr.point);
        return mpfr;
    }

    inline MPFR log2(MPFR mpfr)
    {
        mpfr_log2(mpfr.point, mpfr.point, MPFR_RNDN);
        return mpfr;
    }
    
    // Returns true if the interval has finite endpoints and
    // false if at least one endpoint is infinite
    bool is_finite(const MPFR &mpfr)
    {
        return (mpfr_number_p(mpfr.point)!=0);
    }
}

#endif /* MPFR_tcc */
