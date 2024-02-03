//
//  MPFI.hpp
//  Homotopy Continuation
//
//  Created by burr2 on 8/10/18.
//  Copyright © 2018 burr2. All rights reserved.
//

#ifndef MPFI_hpp
#define MPFI_hpp

#include <iostream>
#include <mpfi.h>
#include <stdio.h>
#include <mpfr.h>
#include "MPFR.hpp"

namespace homotopy
{
    class Interval;
    class MPFR;
    
    class MPFI
    {
    public:
        ///////////////////////
        // Constructors
        ///////////////////////
        
        // Default constructor.  Sets everything equal to 0 using default precision.
        inline MPFI();

        // Constructor with precision.  Sets everything equal to 0 using input precision.
        inline MPFI(mpfr_prec_t);

        // Copy constructor.  Sets upper and lower equal to values in input.
        // Sets precision level equal to that of the input.
        inline MPFI(const MPFI &);
        
        // Copy constructor with precision.  Copies value of one point to another.
        // Sets precision level equal to the user specified precision.
        // Rounds outward, if needed.
        inline MPFI(const MPFI &, mpfr_prec_t);
        
        // Move constructor.  Swaps the data between two inputs.
        inline MPFI(MPFI &&);

        // Parameter constructors.  Sets upper and lower values to input value.
        // MPFI is not templated, so we need individual functions.
        // MPFI provides more functionality, but this should be enough for our purposes.
        
        // Using default precision.
        inline MPFI(int);
        inline MPFI(unsigned int);
        inline MPFI(float);
        inline MPFI(double);
        inline MPFI(const MPFR &);

        // Using user defined precision.
        inline MPFI(int, mpfr_prec_t);
        inline MPFI(unsigned int, mpfr_prec_t);
        inline MPFI(float, mpfr_prec_t);
        inline MPFI(double, mpfr_prec_t);
        inline MPFI(const MPFR &, mpfr_prec_t);
        
        // Assumes that the type can be cast to double.  Using default precision.
        template <typename Type>
        inline MPFI(const Type &);
        
        // Assumes that the type can be cast to double.  Using user defined precision.
        template <typename Type>
        inline MPFI(const Type &, mpfr_prec_t);
        
        // Parameter constructors.  Sets upper and lower values to input values.
        // MPFI is not templated, so we need individual functions.
        // MPFI provides more functionality, but this should be enough for our purposes.
        // The smaller of the two inputs is always the left (order does not matter).
        
        // Using default precision.
        inline MPFI(int, int);
        inline MPFI(unsigned int, unsigned int);
        inline MPFI(float, float);
        inline MPFI(double, double);
        
        // Using maximum precision of the two endpoints.
        inline MPFI(const MPFR &,const MPFR &);
        
        // Using user defined precision.
        inline MPFI(int, int, mpfr_prec_t);
        inline MPFI(unsigned int, unsigned int, mpfr_prec_t);
        inline MPFI(float, float, mpfr_prec_t);
        inline MPFI(double, double, mpfr_prec_t);
        inline MPFI(const MPFR &, const MPFR &, mpfr_prec_t);
        
        // Assumes that the type can be cast to double.  Using default precision.
        template <typename Type1, typename Type2>
        inline MPFI(const Type1 &, const Type2 &);
        
        // Assumes that the type can be cast to double.  Using user defined precision.
        template <typename Type1, typename Type2>
        inline MPFI(const Type1 &, const Type2 &, mpfr_prec_t);

        ///////////////////////
        // Destructor
        ///////////////////////
        
        // Destructor.  Frees memory in interval.
        inline ~MPFI();

        ///////////////////////
        // Operators
        ///////////////////////

        // Assignment operators.  Changes precision to match input precision.
        inline MPFI & operator=(const MPFI &);
        inline MPFI & operator=(const MPFR &);
        
        // Assignment operator for Intervals.  Should be done via conversions.
        inline MPFI & operator=(const Interval &);
        
        // Move-assignment operator.  Swaps the data in the input.
        inline MPFI & operator=(MPFI &&);
        
        // Assignment operators.  Using user defined precision.
        inline MPFI & operator=(int);
        inline MPFI & operator=(unsigned int);
        inline MPFI & operator=(float);
        inline MPFI & operator=(double);
        
        // Assignment operator.
        // Assumes that the type can be cast to double.  Using default precision.
        template <typename Type>
        inline MPFI & operator=(const Type &);

        // In-place addition operator.  Increases precision to maximum of input precisions.
        inline MPFI & operator+=(const MPFI &);
        inline MPFI & operator+=(const MPFR &);
        
        // In-place addition operator for Intervals.  Should be done via conversions.
        inline MPFI & operator+=(const Interval &);
        
        // In-place addition operators.  Using user defined precision.
        inline MPFI & operator+=(int);
        inline MPFI & operator+=(unsigned int);
        inline MPFI & operator+=(float);
        inline MPFI & operator+=(double);
        
        // In-place addition operator.
        // Assumes that the type can be cast to double.  Using default precision.
        template <typename Type>
        inline MPFI & operator+=(const Type &);

        // In-place subtraction operator.  Increases precision to maximum of input precisions.
        inline MPFI & operator-=(const MPFI &);
        inline MPFI & operator-=(const MPFR &);
        
        // In-place subtraction operator for Intervals.  Should be done via conversions.
        inline MPFI & operator-=(const Interval &);
        
        // In-place subtraction operators.  Using user defined precision.
        inline MPFI & operator-=(int);
        inline MPFI & operator-=(unsigned int);
        inline MPFI & operator-=(float);
        inline MPFI & operator-=(double);
        
        // In-place subtraction operator.
        // Assumes that the type can be cast to double.  Using default precision.
        template <typename Type>
        inline MPFI & operator-=(const Type &);

        // In-place multiplication operator.  Increases precision to maximum of input precisions.
        inline MPFI & operator*=(const MPFI &);
        inline MPFI & operator*=(const MPFR &);
        
        // In-place multiplication operators.  Using user defined precision.
        inline MPFI & operator*=(int);
        inline MPFI & operator*=(unsigned int);
        inline MPFI & operator*=(float);
        inline MPFI & operator*=(double);
        
        // In-place multiplication operator.
        // Assumes that the type can be cast to double.  Using default precision.
        template <typename Type>
        inline MPFI & operator*=(const Type &);
 
        // In-place division operator.  Increases precision to maximum of input precisions.
        inline MPFI & operator/=(const MPFI &);
        inline MPFI & operator/=(const MPFR &);
        
        // In-place division operators.  Using user defined precision.
        inline MPFI & operator/=(int);
        inline MPFI & operator/=(unsigned int);
        inline MPFI & operator/=(float);
        inline MPFI & operator/=(double);
        
        // In-place division operator.
        // Assumes that the type can be cast to double.  Using default precision.
        template <typename Type>
        inline MPFI & operator/=(const Type &);
        
        ///////////////////////
        // Interval data
        ///////////////////////
        
        // Returns or sets the width of an MPFI
        // The second function requires fewer constructor calls
        inline friend MPFR width(const MPFI &);
        inline friend bool width(MPFR &, const MPFI &);
        
        // Returns or assigns the upper endpoint of an MPFI
        inline friend MPFR upper(const MPFI &);
        inline friend bool upper(MPFR &, const MPFI &);
        
        // Returns or assigns the lower endpoint of an MPFI
        inline friend MPFR lower(const MPFI &);
        inline friend bool lower(MPFR &, const MPFI &);
        
        // Redurns or sets the midpoint of an MPFI
        inline friend MPFR median(const MPFI &);
        inline friend bool median(MPFR &, const MPFI &);
        
        // Returns or sets largest absolute value of element of an MPFI
        inline friend MPFR abs(const MPFI &);
        inline friend bool abs(MPFR &, const MPFI &);
        
        // Get the square of all the elements of an MPFI
        // The first function save the resulting interval to itself.
        // the second function save the resulting interval to its 1st argument
        inline friend MPFI square_norm(MPFI);
        inline friend bool square_norm(MPFI&, const MPFI &);

        ///////////////////////
        // Output operator
        ///////////////////////
        
        // Display function for an MPFI.  Begins by converting MPFI to a double or pair of doubles.
        // If higher output precision is needed, then we will need a new function.
        inline friend std::ostream& operator<<(std::ostream &, const MPFI &);
        // Display function for an MPFI. Display the MPFI to a given length.
        inline friend bool print_more_digits(const MPFI &, unsigned int );
        ///////////////////////
        // Precision operators
        ///////////////////////
        
        // Returns the current precision.
        inline friend mpfr_prec_t get_precision(const MPFI &);
        
        // Sets the precision of the MPFR (without destroying data).
        // Does not set precision below default 64.
        // Also checks that the precision is within allowable amounts.
        inline friend bool set_precision(MPFI &, mpfr_prec_t);
        
        // Increases the precision to the input precision.
        // If there is only one input, the precision is doubled
        inline friend bool increase_precision(MPFI &);
        inline friend bool increase_precision(MPFI &, mpfr_prec_t);
        
        // Decreases the precision to the input precision.
        // If there is only one input, the precision is halved
        inline friend bool decrease_precision(MPFI &);
        inline friend bool decrease_precision(MPFI &, mpfr_prec_t);
        
        // Get exponent of MPFI
        // Versions for the largest and smallest
        inline friend long max_exponent(MPFI &);
        inline friend long min_exponent(MPFI &);

        ///////////////////////
        // Operations
        ///////////////////////
        
        // Returns or sets exponential of an MPFR.  Rounds to nearest.
        // The power must be an unsigned integer.
        // The first input is set in the second version.
        inline friend MPFI exp(MPFI,unsigned int);
        inline friend bool exp(MPFI &,const MPFI &,unsigned int);
        
        
        // Returns or sets the square root of an MPFR.  Rounds to nearest.
        // The first input is set in the second version.
        inline friend MPFI sqrt(MPFI);
        inline friend bool sqrt(MPFI &,const MPFI &);
        
        ///////////////////////
        // Friend conversion functions
        ///////////////////////
        
        inline friend bool convert(double &, const MPFI &);
        inline friend bool convert(Interval &, const MPFI &);
        inline friend bool convert(MPFI &, const Interval &);

        ///////////////////////
        // Friend Operators
        ///////////////////////
        
        // Subtraction operator.  Using current precision
        inline friend MPFI operator-(const MPFR &, MPFI);
        inline friend MPFI operator-(int, MPFI);
        inline friend MPFI operator-(unsigned int, MPFI);
        inline friend MPFI operator-(float, MPFI);
        inline friend MPFI operator-(double, MPFI);
        
        // Subtraction operator.  Using current precision
        // Assumes the input can be cast to a double.
        template <typename Type>
        inline friend MPFI operator-(const Type &, MPFI);

        // Division operator.  Uses maximum precision of inputs
        inline friend MPFI operator/(const MPFR &, MPFI);
        
        // Division operator.  Using current precision
        inline friend MPFI operator/(int, MPFI);
        inline friend MPFI operator/(unsigned int, MPFI);
        inline friend MPFI operator/(float, MPFI);
        inline friend MPFI operator/(double, MPFI);
        
        // Division operator.  Using current precision
        // Assumes the input can be cast to a double.
        template <typename Type>
        inline friend MPFI operator/(const Type &, MPFI);
        
        // A rewrite of mpfi.h's multiplication to use static mpfr's
        int friend my_mpfi_mul (MPFI&, const MPFI&, const MPFI&);

        ///////////////////////
        // Friend Comparison Operators
        ///////////////////////
        
        // Equality operator.
        inline friend bool operator==(const MPFI &,const MPFI &);
        inline friend bool operator==(const MPFI &,const MPFR &);
        inline friend bool operator==(const MPFR &,const MPFI &);
        
        // Equality operator.
        inline friend bool operator==(const MPFI &,int);
        inline friend bool operator==(int,const MPFI &);
        inline friend bool operator==(const MPFI &,unsigned int);
        inline friend bool operator==(unsigned int, const MPFI &);
        inline friend bool operator==(const MPFI &,float);
        inline friend bool operator==(float, const MPFI &);
        inline friend bool operator==(const MPFI &,double);
        inline friend bool operator==(double, const MPFI &);
        
        // Inequality operator.
        inline friend bool operator!=(const MPFI &,const MPFI &);
        inline friend bool operator!=(const MPFI &,const MPFR &);
        inline friend bool operator!=(const MPFR &,const MPFI &);
        
        // Inequality operator.  May have some problems if the MPFI is NaN.
        inline friend bool operator!=(const MPFI &,int);
        inline friend bool operator!=(int,const MPFI &);
        inline friend bool operator!=(const MPFI &,unsigned int);
        inline friend bool operator!=(unsigned int, const MPFI &);
        inline friend bool operator!=(const MPFI &,float);
        inline friend bool operator!=(float, const MPFI &);
        inline friend bool operator!=(const MPFI &,double);
        inline friend bool operator!=(double, const MPFI &);

        // Less than operator.
        inline friend bool operator<(const MPFI &,const MPFI &);
        inline friend bool operator<(const MPFI &,const MPFR &);
        inline friend bool operator<(const MPFR &,const MPFI &);
        
        // Less than operator.  May have some problems if the MPFR is NaN.
        inline friend bool operator<(const MPFI &,int);
        inline friend bool operator<(int,const MPFI &);
        inline friend bool operator<(const MPFI &,unsigned int);
        inline friend bool operator<(unsigned int, const MPFI &);
        inline friend bool operator<(const MPFI &,float);
        inline friend bool operator<(float, const MPFI &);
        inline friend bool operator<(const MPFI &,double);
        inline friend bool operator<(double, const MPFI &);
        
        // Less than operator.
        // Assumes the input can be cast to a double.
        template <typename Type>
        inline friend bool operator<(const MPFI &,const Type &);
        
        template <typename Type>
        inline friend bool operator<(const Type &,const MPFI &);

        // Less than or equal operator.
        inline friend bool operator<=(const MPFI &,const MPFI &);
        inline friend bool operator<=(const MPFI &,const MPFR &);
        inline friend bool operator<=(const MPFR &,const MPFI &);
 
        // Less than or equal operator.  May have some problems if the MPFI is NaN.
        inline friend bool operator<=(const MPFI &,int);
        inline friend bool operator<=(int,const MPFI &);
        inline friend bool operator<=(const MPFI &,unsigned int);
        inline friend bool operator<=(unsigned int, const MPFI &);
        inline friend bool operator<=(const MPFI &,float);
        inline friend bool operator<=(float, const MPFI &);
        inline friend bool operator<=(const MPFI &,double);
        inline friend bool operator<=(double, const MPFI &);
        
        // Less than or equal operator.
        // Assumes the input can be cast to a double.
        template <typename Type>
        inline friend bool operator<=(const MPFI &,const Type &);
        
        template <typename Type>
        inline friend bool operator<=(const Type &,const MPFI &);
        
        ///////////////////////
        // Friend Inclusion Operators
        ///////////////////////
        
        // Returns true if the first value is contained in the interval
        inline friend bool contains(int, const MPFI&);
        inline friend bool contains(unsigned int, const MPFI&);
        inline friend bool contains(float, const MPFI&);
        inline friend bool contains(double, const MPFI&);
        inline friend bool contains(const MPFR&, const MPFI&);
        
        // Templated version
        template <typename Type>
        inline friend bool contains(const Type&, const MPFI&);
        
        // Returns true if 0 is in the interval
        inline friend bool containsZero(const MPFI&);
        
        // Returns true if the interval has finite endpoints and
        // false if at least one endpoint is infinite
        inline friend bool is_finite(const MPFI &);
        
    private:
        ///////////////////////
        // Data
        ///////////////////////
        
        // MPFI data
        mpfi_t interval;
        
        // Has the data been initialized
        bool initialized;
        
        // static MPFR data for temporaries
        static MPFR temp1, temp2;
    };
    
    ///////////////////////
    // Additional operators
    ///////////////////////
    
    // Addition operator.  Uses maximum precision of inputs
    inline MPFI operator+(MPFI,const MPFI &);
    
    inline MPFI operator+(MPFI,const MPFR &);
    inline MPFI operator+(const MPFR &, MPFI);
    
    // Addition operator for use with Intervals
    inline MPFI operator+(MPFI, const Interval &);
    inline MPFI operator+(const Interval &, MPFI);
    
    // Addition operator.  Using current precision
    inline MPFI operator+(MPFI, int);
    inline MPFI operator+(int, MPFI);
    inline MPFI operator+(MPFI,unsigned int);
    inline MPFI operator+(unsigned int, MPFI);
    inline MPFI operator+(MPFI, float);
    inline MPFI operator+(float, MPFI);
    inline MPFI operator+(MPFI, double);
    inline MPFI operator+(double, MPFI);
    
    // Addition operator.  Using current precision
    // Assumes the input can be cast to a double.
    template <typename Type>
    inline MPFI operator+(MPFI, const Type &);
    
    template <typename Type>
    inline MPFI operator+(const Type &, MPFI);
    
    // Subtraction operator.  Uses maximum precision of inputs
    inline MPFI operator-(MPFI,const MPFI &);
    inline MPFI operator-(MPFI,const MPFR &);
    inline MPFI operator-(MPFI, int);
    inline MPFI operator-(MPFI,unsigned int);
    inline MPFI operator-(MPFI, float);
    inline MPFI operator-(MPFI, double);
    
    // Subtraction operator.  Using current precision
    // Assumes the input can be cast to a double.
    template <typename Type>
    inline MPFI operator-(MPFI, const Type &);
    
    // Multiplication operator.  Uses maximum precision of inputs
    inline MPFI operator*(MPFI,const MPFI &);
    inline MPFI operator*(MPFI,const MPFR &);
    inline MPFI operator*(const MPFR &, MPFI);
    
    // Multiplication operator.  Using current precision
    inline MPFI operator*(MPFI, int);
    inline MPFI operator*(int, MPFI);
    inline MPFI operator*(MPFI,unsigned int);
    inline MPFI operator*(unsigned int, MPFI);
    inline MPFI operator*(MPFI, float);
    inline MPFI operator*(float, MPFI);
    inline MPFI operator*(MPFI, double);
    inline MPFI operator*(double, MPFI);
    
    // Multiplication operator.  Using current precision
    // Assumes the input can be cast to a double.
    template <typename Type>
    inline MPFI operator*(MPFI, const Type &);
    
    template <typename Type>
    inline MPFI operator*(const Type &, MPFI);

    // Division operator.  Uses maximum precision of inputs
    inline MPFI operator/(MPFI,const MPFI &);
    inline MPFI operator/(MPFI,const MPFR &);
    inline MPFI operator/(const MPFR &, MPFI);
    
    // Division operator.  Using current precision
    inline MPFI operator/(MPFI, int);
    inline MPFI operator/(int, MPFI);
    inline MPFI operator/(MPFI,unsigned int);
    inline MPFI operator/(unsigned int, MPFI);
    inline MPFI operator/(MPFI, float);
    inline MPFI operator/(float, MPFI);
    inline MPFI operator/(MPFI, double);
    inline MPFI operator/(double, MPFI);
    
    // Division operator.  Using current precision
    // Assumes the input can be cast to a double.
    template <typename Type>
    inline MPFI operator/(MPFI, const Type &);
    
    ///////////////////////
    // Comparison Operators
    ///////////////////////
    
    // Equality operator.
    // Assumes the input can be cast to a double.
    template <typename Type>
    inline bool operator==(const MPFI &,const Type &);
    
    template <typename Type>
    inline bool operator==(const Type &,const MPFI &);
    
    // Inequality operator.
    // Assumes the input can be cast to a double.
    template <typename Type>
    inline bool operator!=(const MPFI &,const Type &);
    
    template <typename Type>
    inline bool operator!=(const Type &,const MPFI &);
    
    // Greater than operator.
    inline bool operator>(const MPFI &,const MPFI &);
    inline bool operator>(const MPFI &,const MPFR &);
    inline bool operator>(const MPFR &,const MPFI &);
    
    // Greater than operator.  May have some problems if the MPFI is NaN.
    inline bool operator>(const MPFI &,int);
    inline bool operator>(int,const MPFI &);
    inline bool operator>(const MPFI &,unsigned int);
    inline bool operator>(unsigned int, const MPFI &);
    inline bool operator>(const MPFI &,float);
    inline bool operator>(float, const MPFI &);
    inline bool operator>(const MPFI &,double);
    inline bool operator>(double, const MPFI &);
    
    // Greater than operator.
    // Assumes the input can be cast to a double.
    template <typename Type>
    inline bool operator>(const MPFI &,const Type &);
    
    template <typename Type>
    inline bool operator>(const Type &,const MPFI &);
    
    // Greater than or equal operator.
    inline bool operator>=(const MPFI &,const MPFI &);
    inline bool operator>=(const MPFI &,const MPFR &);
    inline bool operator>=(const MPFR &,const MPFI &);
    
    // Greater than or equal operator.  May have some problems if the MPFI is NaN.
    inline bool operator>=(const MPFI &,int);
    inline bool operator>=(int,const MPFI &);
    inline bool operator>=(const MPFI &,unsigned int);
    inline bool operator>=(unsigned int, const MPFI &);
    inline bool operator>=(const MPFI &,float);
    inline bool operator>=(float, const MPFI &);
    inline bool operator>=(const MPFI &,double);
    inline bool operator>=(double, const MPFI &);
    
    // Greater than or equal operator.
    // Assumes the input can be cast to a double.
    template <typename Type>
    inline bool operator>=(const MPFI &,const Type &);
    
    template <typename Type>
    inline bool operator>=(const Type &,const MPFI &);
}

// The following is included because the functions are encouraged to be inline or templated.
#include "MPFI.tcc"

#endif /* MPFI_hpp */
