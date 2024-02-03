//
//  Complex.hpp
//  Homotopy Continuation
//
//  Created by burr2 on 8/17/18.
//  Copyright © 2018 burr2. All rights reserved.
//

#ifndef Complex_hpp
#define Complex_hpp

#include <iostream>
#include <mpfr.h>
#include "NumberUtils.hpp"

namespace homotopy
{
    template <typename Type>
    class Complex
    {
        template <typename Type1> friend class Complex;
        
    public:
        ///////////////////////
        // Constructors
        ///////////////////////
        
        // Default constructor.  Sets complex number equal to 0.
        inline Complex();
        
        // Constructor with precision.  Sets complex number equal to 0 using input precision.
        inline Complex(mpfr_prec_t);
        
        // Copy constructor.  Copies values in input complex number
        // Uses copy constructor of Type
        template <typename Type1>
        inline Complex(const Complex<Type1> &);
        
        // Copy constructor with precision.  Copies values in input complex number
        // Uses copy constructor of Type
        template <typename Type1>
        inline Complex(const Complex<Type1> &, mpfr_prec_t);
        
        // Parameter constructors.
        // Sets complex number equal to input (as a real number)
        // Should be slightly faster than pass by reference
        inline Complex(int);
        inline Complex(unsigned int);
        inline Complex(float);
        inline Complex(double);
        
        // Templated constructor.
        template <typename Number>
        inline Complex(const Number &);
        
        // Sets complex number equal to input (as a real number)
        // Uses user defined precision
        inline Complex(int, mpfr_prec_t);
        inline Complex(unsigned int, mpfr_prec_t);
        inline Complex(float, mpfr_prec_t);
        inline Complex(double, mpfr_prec_t);
        
        // Templated constructor.
        template <typename Number>
        inline Complex(const Number &, mpfr_prec_t);
        
        // Parameter constructors.
        // Sets complex number equal to input as a pair
        // Should be slightly faster than pass by reference
        inline Complex(int, int);
        inline Complex(unsigned int, unsigned int);
        inline Complex(float, float);
        inline Complex(double, double);
        
        // Templated constructor.
        template <typename Number1, typename Number2>
        inline Complex(const Number1 &,const Number2 &);
        
        // Sets complex number equal to input (as a real number)
        // Uses user defined precision
        inline Complex(int, int, mpfr_prec_t);
        inline Complex(unsigned int, unsigned int, mpfr_prec_t);
        inline Complex(float, float, mpfr_prec_t);
        inline Complex(double, double, mpfr_prec_t);
        
        // Templated constructor.
        template <typename Number1, typename Number2>
        inline Complex(const Number1 &, const Number2 &, mpfr_prec_t);
        
        ///////////////////////
        // Operators
        ///////////////////////
        
        // Reassign the real part or the imag part
        inline Complex<Type> & set_real(int);
        inline Complex<Type> & set_real(unsigned int);
        inline Complex<Type> & set_real(float);
        inline Complex<Type> & set_real(double);

        template <typename Type1>
        inline Complex<Type> & set_real(Type1 &);

        inline Complex<Type> & set_imag(int);
        inline Complex<Type> & set_imag(unsigned int);
        inline Complex<Type> & set_imag(float);
        inline Complex<Type> & set_imag(double);

        template <typename Type1>
        inline Complex<Type> & set_imag(Type1 &);

        // Assignment operators.  Uses assignment of Type
        template <typename Type1>
        inline Complex<Type> & operator=(const Complex<Type1> &);
        
        // Assignment operators to (real) input.
        inline Complex<Type> & operator=(int);
        inline Complex<Type> & operator=(unsigned int);
        inline Complex<Type> & operator=(float);
        inline Complex<Type> & operator=(double);
        
        // Templated assignment operator
        // Assumes that the type of Number can be assigned to Type.
        template <typename Number>
        inline Complex<Type> & operator=(const Number &);
        
        // In-place addition operator.
        template <typename Number>
        inline Complex<Type> & operator+=(const Complex<Number> &);
        
        // Addition operator to (real) input.
        inline Complex<Type> & operator+=(int);
        inline Complex<Type> & operator+=(unsigned int);
        inline Complex<Type> & operator+=(float);
        inline Complex<Type> & operator+=(double);
        
        // Templated addition operator
        // Assumes that the type of Number can be added to Type.
        template <typename Number>
        inline Complex<Type> & operator+=(const Number &);
        
        // In-place subtraction operator.
        template <typename Number>
        inline Complex<Type> & operator-=(const Complex<Number> &);
        
        // Subtraction operator to (real) input.
        inline Complex<Type> & operator-=(int);
        inline Complex<Type> & operator-=(unsigned int);
        inline Complex<Type> & operator-=(float);
        inline Complex<Type> & operator-=(double);
        
        // Templated subtraction operator
        // Assumes that the type of Number can be added to Type.
        template <typename Number>
        inline Complex<Type> & operator-=(const Number &);
        
        // In-place multiplication operator.
        template <typename Number>
        inline Complex<Type> & operator*=(const Complex<Number> &);
        
        // Multiplication operator to (real) input.
        inline Complex<Type> & operator*=(int);
        inline Complex<Type> & operator*=(unsigned int);
        inline Complex<Type> & operator*=(float);
        inline Complex<Type> & operator*=(double);
        
        // Templated multiplication operator
        // Assumes that the type of Number can be added to Type.
        template <typename Number>
        inline Complex<Type> & operator*=(const Number &);
        
        // In-place division operator.
        template <typename Number>
        inline Complex<Type> & operator/=(const Complex<Number> &);
        
        // Division operator to (real) input.
        inline Complex<Type> & operator/=(int);
        inline Complex<Type> & operator/=(unsigned int);
        inline Complex<Type> & operator/=(float);
        inline Complex<Type> & operator/=(double);
        
        // Templated division operator
        // Assumes that the type of Number can be added to Type.
        template <typename Number>
        inline Complex<Type> & operator/=(const Number &);
        
        ///////////////////////
        // Output operator
        ///////////////////////
        
        // Display function for a Complex number.
        template <typename Number>
        inline friend std::ostream& operator<<(std::ostream &, const Complex<Number> &);
        
        ///////////////////////
        // Precision operators
        ///////////////////////
        
        // Returns the current precision.
        template <typename Number>
        inline friend mpfr_prec_t get_precision(const Complex<Number> &);
        
        // Sets the precision of the complex (without destroying data).
        template <typename Number>
        inline friend bool set_precision(Complex<Number> &, mpfr_prec_t);
        
        // Increases the precision to the input precision.
        // If there is only one input, the precision is doubled
        template <typename Number>
        inline friend bool increase_precision(Complex<Number> &);
        
        template <typename Number>
        inline friend bool increase_precision(Complex<Number> &, mpfr_prec_t);
        
        // Decreases the precision to the input precision.
        // If there is only one input, the precision is halved
        template <typename Number>
        inline friend bool decrease_precision(Complex<Number> &);
        
        template <typename Number>
        inline friend bool decrease_precision(Complex<Number> &, mpfr_prec_t);
        
        // Get exponent of largest element.  An approximation to the exponent of norm.
        template <typename Number>
        inline friend long exponent(Complex<Number> &);
        
        // Versions for largest and smallest.
        // max_exponent is the same as exponent.
        template <typename Number>
        inline friend long max_exponent(Complex<Number> &);
        
        template <typename Number>
        inline friend long min_exponent(Complex<Number> &);

        ///////////////////////
        // Operations
        ///////////////////////
 
        // Returns or sets the square of the norm of a complex number
        template <typename Number>
        inline friend Number norm(const Complex<Number> &);
        
        template <typename Number>
        inline friend bool norm(Number &, const Complex<Number> &);
        
        // Returns or sets the real part of a complex number
        template <typename Number>
        inline friend Number real(const Complex<Number> &);
        
        template <typename Number>
        inline friend bool real(Number &, const Complex<Number> &);
        
        // Returns or sets the imaginary part of a complex number
        template <typename Number>
        inline friend Number imag(const Complex<Number> &);
        
        template <typename Number>
        inline friend bool imag(Number &,const Complex<Number> &);

        // Returns or sets exponential of a complex number.
        // The power must be an unsigned integer.
        template <typename Number>
        inline friend Complex<Number> exp(Complex<Number>,unsigned int);
        
        template <typename Number>
        inline friend bool exp(Complex<Number> &,const Complex<Number> &, unsigned int);

        ///////////////////////
        // Friend Operators
        ///////////////////////

        // Division operator.
        template <typename Number>
        inline friend Complex<Number> operator/(int, Complex<Number>);
        
        template <typename Number>
        inline friend Complex<Number> operator/(unsigned int, Complex<Number>);
        
        template <typename Number>
        inline friend Complex<Number> operator/(float, Complex<Number>);
        
        template <typename Number>
        inline friend Complex<Number> operator/(double, Complex<Number>);
        
        // Templated division operator
        // Assumes the second type can be multiplied by and divided by the first.
        template <typename Number1, typename Number2>
        inline friend Complex<Number1> operator/(const Number2 &,Complex<Number1>);
        
        ///////////////////////
        // Comparison Operators
        ///////////////////////
        
        // Equality operator.  Assumes the first type can be compared with the second
        template <typename Number1,typename Number2>
        inline friend bool operator==(const Complex<Number1> &, const Complex<Number2> &);
        
        // Equality operator.  Compares to (real) variable
        template <typename Number>
        inline friend bool operator==(const Complex<Number> &,int);
        
        template <typename Number>
        inline friend bool operator==(int, const Complex<Number> &);
        
        template <typename Number>
        inline friend bool operator==(const Complex<Number> &,unsigned int);
        
        template <typename Number>
        inline friend bool operator==(unsigned int, const Complex<Number> &);
        
        template <typename Number>
        inline friend bool operator==(const Complex<Number> &,float);
        
        template <typename Number>
        inline friend bool operator==(float, const Complex<Number> &);
        
        template <typename Number>
        inline friend bool operator==(const Complex<Number> &, double);
        
        template <typename Number>
        inline friend bool operator==(double, const Complex<Number> &);
        
        // Equality operator.
        // Assumes the input can be compared.
        template <typename Number1, typename Number2>
        inline friend bool operator==(const Complex<Number1> &,const Number2 &);
        
        template <typename Number1, typename Number2>
        inline friend bool operator==(Number2, const Complex<Number1> &);
        
        ///////////////////////
        // Friend conversion operators
        ///////////////////////
        
        template <typename Type1, typename Type2>
        inline friend bool convert(Complex<Type1>&, const Type2&);
        
        template <typename Type1, typename Type2>
        inline friend bool convert(Complex<Type1> &, const Complex<Type2> &);
        
        ///////////////////////
        // Friend conversion operators
        ///////////////////////
        
        template <typename Type1>
        inline friend bool containsZero(const Complex<Type1> &);
        
        template <typename Type1>
        inline friend bool is_finite(const Complex<Type1> &);

    private:
        ///////////////////////
        // Data
        ///////////////////////
        
        // Real and complex parts
        Type real, imag;
        
        static Type temp1, temp2;
    };
    
    ///////////////////////
    // Operators
    ///////////////////////
    
    // Addition operator.  Assumes that the second type can be added to the first type.
    template <typename Number1, typename Number2>
    inline Complex<Number1> operator+(Complex<Number1>, const Complex<Number2> &);
    
    // Addition operator.
    template <typename Number>
    inline Complex<Number> operator+(Complex<Number>,int);
    
    template <typename Number>
    inline Complex<Number> operator+(int, Complex<Number>);
    
    template <typename Number>
    inline Complex<Number> operator+(Complex<Number>,unsigned int);
    
    template <typename Number>
    inline Complex<Number> operator+(unsigned int, Complex<Number>);
    
    template <typename Number>
    inline Complex<Number> operator+(Complex<Number>, float);
    
    template <typename Number>
    inline Complex<Number> operator+(float, Complex<Number>);
    
    template <typename Number>
    inline Complex<Number> operator+(Complex<Number>, double);
    
    template <typename Number>
    inline Complex<Number> operator+(double, Complex<Number>);
    
    // Templated addition operator
    // Assumes the second type can be added to the first.
    template <typename Number1, typename Number2>
    inline Complex<Number1> operator+(const Number2 &,Complex<Number1>);
    
    template <typename Number1, typename Number2>
    inline Complex<Number1> operator+(Complex<Number1>,const Number2 &);
    
    // Subtraction operator.  Assumes that the second type can be added to the first type.
    template <typename Number1, typename Number2>
    inline Complex<Number1> operator-(Complex<Number1>, const Complex<Number2> &);
    
    // Subtraction operator.
    template <typename Number>
    inline Complex<Number> operator-(Complex<Number>,int);
    
    template <typename Number>
    inline Complex<Number> operator-(int, Complex<Number>);
    
    template <typename Number>
    inline Complex<Number> operator-(Complex<Number>,unsigned int);
    
    template <typename Number>
    inline Complex<Number> operator-(unsigned int, Complex<Number>);
    
    template <typename Number>
    inline Complex<Number> operator-(Complex<Number>, float);
    
    template <typename Number>
    inline Complex<Number> operator-(float, Complex<Number>);
    
    template <typename Number>
    inline Complex<Number> operator-(Complex<Number>, double);
    
    template <typename Number>
    inline Complex<Number> operator-(double, Complex<Number>);
    
    // Templated subtraction operator
    // Assumes the second type can be subtracted from the first.
    template <typename Number1, typename Number2>
    inline Complex<Number1> operator-(const Number2 &,Complex<Number1>);
    
    template <typename Number1, typename Number2>
    inline Complex<Number1> operator-(Complex<Number1>,const Number2 &);
    
    // Multiplication operator.  Assumes that the second type can be added to the first type.
    template <typename Number1, typename Number2>
    inline Complex<Number1> operator*(Complex<Number1>, const Complex<Number2> &);
    
    // Multiplication operator.
    template <typename Number>
    inline Complex<Number> operator*(Complex<Number>,int);
    
    template <typename Number>
    inline Complex<Number> operator*(int, Complex<Number>);
    
    template <typename Number>
    inline Complex<Number> operator*(Complex<Number>,unsigned int);
    
    template <typename Number>
    inline Complex<Number> operator*(unsigned int, Complex<Number>);
    
    template <typename Number>
    inline Complex<Number> operator*(Complex<Number>, float);
    
    template <typename Number>
    inline Complex<Number> operator*(float, Complex<Number>);
    
    template <typename Number>
    inline Complex<Number> operator*(Complex<Number>, double);
    
    template <typename Number>
    inline Complex<Number> operator*(double, Complex<Number>);
    
    // Templated multiplication operator
    // Assumes the second type can be multiplied to the first.
    template <typename Number1, typename Number2>
    inline Complex<Number1> operator*(const Number2 &,Complex<Number1>);
    
    template <typename Number1, typename Number2>
    inline Complex<Number1> operator*(Complex<Number1>,const Number2 &);
    
    // Division operator.  Assumes that the second type can be added to the first type.
    template <typename Number1, typename Number2>
    inline Complex<Number1> operator/(Complex<Number1>, const Complex<Number2> &);
    
    // Division operator.
    template <typename Number>
    inline Complex<Number> operator/(Complex<Number>,int);
    
    template <typename Number>
    inline Complex<Number> operator/(Complex<Number>,unsigned int);
    
    template <typename Number>
    inline Complex<Number> operator/(Complex<Number>, float);
    
    template <typename Number>
    inline Complex<Number> operator/(Complex<Number>, double);
    
    // Templated division operator
    // Assumes the second type can be multiplied by and divided by the first.
    template <typename Number1, typename Number2>
    inline Complex<Number1> operator/(Complex<Number1>,const Number2 &);

    // Inequality operator.  Assumes the first type can be compared with the second
    template <typename Number1,typename Number2>
    inline bool operator!=(const Complex<Number1> &, const Complex<Number2> &);
    
    // Inequality operator.  Compares to (real) variable
    template <typename Number>
    inline bool operator!=(const Complex<Number> &,int);
    
    template <typename Number>
    inline bool operator!=(int, const Complex<Number> &);
    
    template <typename Number>
    inline bool operator!=(const Complex<Number> &,unsigned int);
    
    template <typename Number>
    inline bool operator!=(unsigned int, const Complex<Number> &);
    
    template <typename Number>
    inline bool operator!=(const Complex<Number> &,float);
    
    template <typename Number>
    inline bool operator!=(float, const Complex<Number> &);
    
    template <typename Number>
    inline bool operator!=(const Complex<Number> &, double);
    
    template <typename Number>
    inline bool operator!=(double, const Complex<Number> &);
    
    // Inequality operator.
    // Assumes the input can be compared.
    template <typename Number1, typename Number2>
    inline bool operator!=(const Complex<Number1> &,const Number2 &);
    
    template <typename Number1, typename Number2>
    inline bool operator!=(Number2, const Complex<Number1> &);
}

// The following is included because the functions are encouraged to be inline or templated.
#include "Complex.tcc"

#endif /* Complex_hpp */
