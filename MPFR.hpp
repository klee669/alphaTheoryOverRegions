//
//  MPFR.hpp
//  Homotopy Continuation
//
//  Created by burr2 on 8/8/18.
//  Copyright © 2018 burr2. All rights reserved.
//

#ifndef MPFR_hpp
#define MPFR_hpp

#include <iostream>
#include <mpfr.h>

namespace homotopy
{
    class MPFI;
    class Interval;
    
    class MPFR
    {
        friend class MPFI;
        
    public:
        ///////////////////////
        // Constructors
        ///////////////////////
        
        // Default constructor.  Sets value to 0 using default precision.
        inline MPFR();
        
        // Constructor with precision.  Sets value to 0 using input precision.
        inline MPFR(mpfr_prec_t);
        
        // Copy constructor.  Copies value of one point to another.
        // Sets precision level equal to that of the input.
        inline MPFR(const MPFR &);
        
        // Copy constructor with precision.  Copies value of one point to another.
        // Sets precision level equal to the user specified precision.
        inline MPFR(const MPFR &, mpfr_prec_t);
        
        // Move constructor.  Swaps the data between two inputs.
        inline MPFR(MPFR &&);
        
        // Parameter constructors.
        // MPFR is not templated, so we need individual functions.
        // MPFR provides more functionality, but this should be enough for our purposes.
        
        // Using default precision.
        inline MPFR(int);
        inline MPFR(unsigned int);
        inline MPFR(float);
        inline MPFR(double);
        
        // Using user defined precision.
        inline MPFR(int, mpfr_prec_t);
        inline MPFR(unsigned int, mpfr_prec_t);
        inline MPFR(float, mpfr_prec_t);
        inline MPFR(double, mpfr_prec_t);
        inline MPFR(const char *, mpfr_prec_t);
        inline MPFR(const char *);
      inline MPFR(std::string, mpfr_prec_t);
      inline MPFR(std::string);
        
        // Assumes that the type can be cast to double.  Using default precision.
        template <typename Type>
        inline MPFR(const Type &);
        
        // Assumes that the type can be cast to double.  Using user defined precision.
        template <typename Type>
        inline MPFR(const Type &, mpfr_prec_t);
  
        ///////////////////////
        // Destructor
        ///////////////////////
  
        // Destructor.  Frees memory in point.
        inline ~MPFR();
  
        ///////////////////////
        // Operators
        ///////////////////////
        
        // Assignment operator.  Changes precision to match input precision.
        inline MPFR & operator=(const MPFR &);
        
        // Move-assignment operator.  Swaps the data in the input.
        inline MPFR & operator=(MPFR &&);
        
        // Assignment operators.  Using current defined precision.
        inline MPFR & operator=(int);
        inline MPFR & operator=(unsigned int);
        inline MPFR & operator=(float);
        inline MPFR & operator=(double);

        // Assignment operators.  Using mpfr_t.
        inline MPFR & operator=(const mpfr_t &);
        
        // Assignment operator.
        // Assumes that the type can be cast to double.  Using current precision.
        template <typename Type>
        inline MPFR & operator=(const Type &);
        
        // In-place addition operator.  Increases precision to maximum of input precisions.
        inline MPFR & operator+=(const MPFR &);
        
        // In-place addition operators.  Using current defined precision.
        inline MPFR & operator+=(int);
        inline MPFR & operator+=(unsigned int);
        inline MPFR & operator+=(float);
        inline MPFR & operator+=(double);
        
        // In-place addition operator.
        // Assumes that the type can be cast to double.  Using current precision.
        template <typename Type>
        inline MPFR & operator+=(const Type &);
        
        // In-place subtraction operator.  Increases precision to maximum of input precisions.
        inline MPFR & operator-=(const MPFR &);
        
        // In-place subtraction operators.  Using current defined precision.
        inline MPFR & operator-=(int);
        inline MPFR & operator-=(unsigned int);
        inline MPFR & operator-=(float);
        inline MPFR & operator-=(double);
        
        // In-place subtraction operator.
        // Assumes that the type can be cast to double.  Using current precision.
        template <typename Type>
        inline MPFR & operator-=(const Type &);
        
        // In-place multiplication operator.  Increases precision to maximum of input precisions.
        inline MPFR & operator*=(const MPFR &);
        
        // In-place multiplication operators.  Using current defined precision.
        inline MPFR & operator*=(int);
        inline MPFR & operator*=(unsigned int);
        inline MPFR & operator*=(float);
        inline MPFR & operator*=(double);
        
        // In-place multiplication operator.
        // Assumes that the type can be cast to double.  Using current precision.
        template <typename Type>
        inline MPFR & operator*=(const Type &);
        
        // In-place division operator.  Increases precision to maximum of input precisions.
        inline MPFR & operator/=(const MPFR &);
        
        // In-place division operators.  Using current defined precision.
        inline MPFR & operator/=(int);
        inline MPFR & operator/=(unsigned int);
        inline MPFR & operator/=(float);
        inline MPFR & operator/=(double);
        
        // In-place division operator.
        // Assumes that the type can be cast to double.  Using current precision.
        template <typename Type>
        inline MPFR & operator/=(const Type &);

        ///////////////////////
        // Output operator
        ///////////////////////
        
        // Display function for an MPFR.  Begins by converting MPFR to a double.
        // If higher output precision is needed, then we will need a new function.
        inline friend std::ostream& operator<<(std::ostream &, const MPFR &);

        // Display an MPFR to a given length.
        inline friend bool print_more_digits(const MPFR &, unsigned int );

        // Display an MPFI to a given length.
        inline friend bool print_more_digits(const MPFI &, unsigned int );
 
        ///////////////////////
        // Precision operators
        ///////////////////////
        
        // Returns the current precision.
        inline friend mpfr_prec_t get_precision(const MPFR &);
        
        // Sets the precision of the MPFR (without destroying data).
        // Does not set precision below default (IEEE standard) 53.
        // Also checks that the precision is within allowable amounts.
        inline friend bool set_precision(MPFR &, mpfr_prec_t);
        
        // Increases the precision to the input precision.
        // If there is only one input, the precision is doubled
        inline friend bool increase_precision(MPFR &);
        inline friend bool increase_precision(MPFR &, mpfr_prec_t);
        
        // Decreases the precision to the input precision.
        // If there is only one input, the precision is halved
        inline friend bool decrease_precision(MPFR &);
        inline friend bool decrease_precision(MPFR &, mpfr_prec_t);
        
        // Get exponent of MPFR (mantissa is assumed to be between 1/2 and 1)
        inline friend long exponent(MPFR &);
        
        ///////////////////////
        // Operations
        ///////////////////////
        inline friend MPFR square_norm(MPFR);
        inline friend bool square_norm(MPFR &, const MPFR &);

      
        // Returns or assigns the absolute value of the MPFR
        inline friend MPFR abs(MPFR);
        inline friend bool abs(MPFR &,const MPFR &);
        
        // Returns or assigns exponential of an MPFR.  Rounds to nearest.
        // The power must be an unsigned integer.
        inline friend MPFR exp(MPFR,unsigned int);
        inline friend bool exp(MPFR &, const MPFR &, unsigned int);
        
        // Returns or assigns the square root of an MPFR.  Rounds to nearest.
        inline friend MPFR sqrt(MPFR);
        inline friend bool sqrt(MPFR &,const MPFR &);
        
        ///////////////////////
        // Friend Interval data
        ///////////////////////
        
        // Functions from MPFI which need to be friends (for efficiency).
        inline friend MPFR width(const MPFI &);
        inline friend bool width(MPFR &, const MPFI &);
        
        inline friend MPFR upper(const MPFI &);
        inline friend bool upper(MPFR &, const MPFI &);
        
        inline friend MPFR lower(const MPFI &);
        inline friend bool lower(MPFR &, const MPFI &);
        
        inline friend MPFR median(const MPFI &);
        inline friend bool median(MPFR &, const MPFI &);
        
        inline friend MPFR abs(const MPFI &);
        inline friend bool abs(MPFR &, const MPFI &);
        
        inline friend MPFI square_norm(const MPFI &);
        inline friend bool square_norm(MPFI &, const MPFI &);
        
        ///////////////////////
        // Friend conversion functions
        ///////////////////////
        
        inline friend bool convert(double &, const MPFR &);
        inline friend bool convert(Interval &, const MPFR &);

        ///////////////////////
        // Friend Operators
        ///////////////////////
        
        // Subtraction operator.  Using current precision
        inline friend MPFR operator-(int, MPFR);
        inline friend MPFR operator-(unsigned int, MPFR);
        inline friend MPFR operator-(float, MPFR);
        inline friend MPFR operator-(double, MPFR);
        
        // Subtraction operator.  Using current precision
        // Assumes the input can be cast to a double.
        template <typename Type>
        inline friend MPFR operator-(const Type &, MPFR);
        
        // Division operator.  Using current precision
        inline friend MPFR operator/(int, MPFR);
        inline friend MPFR operator/(unsigned int, MPFR);
        inline friend MPFR operator/(float, MPFR);
        inline friend MPFR operator/(double, MPFR);
        
        // Division operator.  Using current precision
        // Assumes the input can be cast to a double.
        template <typename Type>
        inline friend MPFR operator/(const Type &, MPFR);
        
        ///////////////////////
        // Friend comparison Operators
        ///////////////////////
        
        // Equality operator.
        inline friend bool operator==(const MPFR &,const MPFR &);
        
        // Equality operator.
        inline friend bool operator==(const MPFR &,int);
        inline friend bool operator==(int,const MPFR &);
        inline friend bool operator==(const MPFR &,unsigned int);
        inline friend bool operator==(unsigned int, const MPFR &);
        inline friend bool operator==(const MPFR &,float);
        inline friend bool operator==(float, const MPFR &);
        inline friend bool operator==(const MPFR &,double);
        inline friend bool operator==(double, const MPFR &);
        
        // Equality operator.
        // Assumes the input can be cast to a double.
        template <typename Type>
        inline friend bool operator==(const MPFR &,const Type &);
        
        template <typename Type>
        inline friend bool operator==(const Type &,const MPFR &);
        
        // Inequality operator.
        inline friend bool operator!=(const MPFR &,const MPFR &);
        
        // Inequality operator.
        inline friend bool operator!=(const MPFR &,int);
        inline friend bool operator!=(int,const MPFR &);
        inline friend bool operator!=(const MPFR &,unsigned int);
        inline friend bool operator!=(unsigned int, const MPFR &);
        inline friend bool operator!=(const MPFR &,float);
        inline friend bool operator!=(float, const MPFR &);
        inline friend bool operator!=(const MPFR &,double);
        inline friend bool operator!=(double, const MPFR &);
        
        // Inequality operator.
        // Assumes the input can be cast to a double.
        template <typename Type>
        inline friend bool operator!=(const MPFR &,const Type &);
        
        template <typename Type>
        inline friend bool operator!=(const Type &,const MPFR &);

        // Less than operator.
        inline friend bool operator<(const MPFR &,const MPFR &);
        
        // Less than operator.
        inline friend bool operator<(const MPFR &,int);
        inline friend bool operator<(int,const MPFR &);
        inline friend bool operator<(const MPFR &,unsigned int);
        inline friend bool operator<(unsigned int, const MPFR &);
        inline friend bool operator<(const MPFR &,float);
        inline friend bool operator<(float, const MPFR &);
        inline friend bool operator<(const MPFR &,double);
        inline friend bool operator<(double, const MPFR &);
        
        // Less than operator.
        // Assumes the input can be cast to a double.
        template <typename Type>
        inline friend bool operator<(const MPFR &,const Type &);
        
        template <typename Type>
        inline friend bool operator<(const Type &,const MPFR &);
        
        // Less than or equal operator.
        inline friend bool operator<=(const MPFR &,const MPFR &);
        
        // Less than or equal operator.
        inline friend bool operator<=(const MPFR &,int);
        inline friend bool operator<=(int,const MPFR &);
        inline friend bool operator<=(const MPFR &,unsigned int);
        inline friend bool operator<=(unsigned int, const MPFR &);
        inline friend bool operator<=(const MPFR &,float);
        inline friend bool operator<=(float, const MPFR &);
        inline friend bool operator<=(const MPFR &,double);
        inline friend bool operator<=(double, const MPFR &);
        
        // Less than or equal operator.
        // Assumes the input can be cast to a double.
        template <typename Type>
        inline friend bool operator<=(const MPFR &,const Type &);
        
        template <typename Type>
        inline friend bool operator<=(const Type &,const MPFR &);

        // Greater than operator.
        inline friend bool operator>(const MPFR &,const MPFR &);
        
        // Greater than operator.
        inline friend bool operator>(const MPFR &,int);
        inline friend bool operator>(int,const MPFR &);
        inline friend bool operator>(const MPFR &,unsigned int);
        inline friend bool operator>(unsigned int, const MPFR &);
        inline friend bool operator>(const MPFR &,float);
        inline friend bool operator>(float, const MPFR &);
        inline friend bool operator>(const MPFR &,double);
        inline friend bool operator>(double, const MPFR &);
        
        // Greater than operator.
        // Assumes the input can be cast to a double.
        template <typename Type>
        inline friend bool operator>(const MPFR &,const Type &);
        
        template <typename Type>
        inline friend bool operator>(const Type &,const MPFR &);
        
        // Greater than or equal operator.
        inline friend bool operator>=(const MPFR &,const MPFR &);
        
        // Greater than or equal operator.  
        inline friend bool operator>=(const MPFR &,int);
        inline friend bool operator>=(int,const MPFR &);
        inline friend bool operator>=(const MPFR &,unsigned int);
        inline friend bool operator>=(unsigned int, const MPFR &);
        inline friend bool operator>=(const MPFR &,float);
        inline friend bool operator>=(float, const MPFR &);
        inline friend bool operator>=(const MPFR &,double);
        inline friend bool operator>=(double, const MPFR &);
        
        // Greater than or equal operator.
        // Assumes the input can be cast to a double.
        template <typename Type>
        inline friend bool operator>=(const MPFR &,const Type &);
        
        template <typename Type>
        inline friend bool operator>=(const Type &,const MPFR &);
        
        ///////////////////////
        // Friend MPFI Operators
        ///////////////////////
        
        inline friend MPFI operator-(const MPFR &, MPFI);
        inline friend MPFI operator/(const MPFR &,MPFI);
        inline friend bool operator==(const MPFI &,const MPFR &);
        inline friend bool operator==(const MPFR &,const MPFI &);
        inline friend bool operator<(const MPFI &,const MPFR &);
        inline friend bool operator<(const MPFR &,const MPFI &);
        inline friend bool operator<=(const MPFI &,const MPFR &);
        inline friend bool operator<=(const MPFR &,const MPFI &);
        inline friend bool operator>(const MPFI &,const MPFR &);
        inline friend bool operator>(const MPFR &,const MPFI &);
        inline friend bool operator>=(const MPFI &,const MPFR &);
        inline friend bool operator>=(const MPFR &,const MPFI &);
        
        ///////////////////////
        // Friend MPFI Comparisons
        ///////////////////////
        
        inline friend bool contains(const MPFR&, const MPFI&);
        
        inline friend bool round_up(double &,const MPFR &);
        inline friend bool round_up(MPFR &, double);
        inline friend bool round_up(MPFR &, const MPFR &);
        
        inline friend bool round_down(double &,const MPFR &);
        inline friend bool round_down(MPFR &, double);
        inline friend bool round_down(MPFR &, const MPFR &);
        
        // A rewrite of mpfi.h's multiplication to use static mpfr's
        int friend my_mpfi_mul (MPFI&, const MPFI&, const MPFI&);

        inline friend MPFR floor(MPFR);
        
        inline friend MPFR log2(MPFR);
        
        // Returns true if the interval has finite endpoints and
        // false if at least one endpoint is infinite
        inline friend bool is_finite(const MPFR &);
        
    private:
        ///////////////////////
        // Data
        ///////////////////////
        
        // MPFR data
        mpfr_t point;
        
        // Has the data been initialized
        bool initialized;
    };
    
    ///////////////////////
    // Additional Operators
    ///////////////////////
    
    // Addition operator.  Uses maximum precision of inputs
    inline MPFR operator+(MPFR,const MPFR &);
    
    // Addition operator.  Using current precision
    inline MPFR operator+(MPFR, int);
    inline MPFR operator+(int, MPFR);
    inline MPFR operator+(MPFR,unsigned int);
    inline MPFR operator+(unsigned int, MPFR);
    inline MPFR operator+(MPFR, float);
    inline MPFR operator+(float, MPFR);
    inline MPFR operator+(MPFR, double);
    inline MPFR operator+(double, MPFR);
    
    // Addition operator.  Using current precision
    // Assumes the input can be cast to a double.
    template <typename Type>
    inline MPFR operator+(MPFR, const Type &);
    
    template <typename Type>
    inline MPFR operator+(const Type &, MPFR);
    
    // Subtraction operator.  Uses maximum precision of inputs
    inline MPFR operator-(MPFR,const MPFR &);
    
    // Subtraction operator.  Using current precision
    inline MPFR operator-(MPFR, int);
    inline MPFR operator-(MPFR,unsigned int);
    inline MPFR operator-(MPFR, float);
    inline MPFR operator-(MPFR, double);
    
    // Subtraction operator.  Using current precision
    // Assumes the input can be cast to a double.
    template <typename Type>
    inline MPFR operator-(MPFR, const Type &);
    
    // Multiplication operator.  Uses maximum precision of inputs
    inline MPFR operator*(MPFR,const MPFR &);
    
    // Multiplication operator.  Using current precision
    inline MPFR operator*(MPFR, int);
    inline MPFR operator*(int, MPFR);
    inline MPFR operator*(MPFR,unsigned int);
    inline MPFR operator*(unsigned int, MPFR);
    inline MPFR operator*(MPFR, float);
    inline MPFR operator*(float, MPFR);
    inline MPFR operator*(MPFR, double);
    inline MPFR operator*(double, MPFR);
    
    // Multiplication operator.  Using current precision
    // Assumes the input can be cast to a double.
    template <typename Type>
    inline MPFR operator*(MPFR, const Type &);
    
    template <typename Type>
    inline MPFR operator*(const Type &, MPFR);
    
    // Division operator.  Uses maximum precision of inputs
    inline MPFR operator/(MPFR,const MPFR &);
    
    // Division operator.  Using current precision
    inline MPFR operator/(MPFR, int);
    inline MPFR operator/(MPFR,unsigned int);
    inline MPFR operator/(MPFR, float);
    inline MPFR operator/(MPFR, double);
    
    // Division operator.  Using current precision
    // Assumes the input can be cast to a double.
    template <typename Type>
    inline MPFR operator/(MPFR, const Type &);
}

// The following is included because the functions are encouraged to be inline or templated.
#include "MPFR.tcc"

#endif /* MPFR_hpp */
