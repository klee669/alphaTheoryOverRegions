//
//  Complex.tcc
//  Homotopy Continuation
//
//  Created by burr2 on 8/17/18.
//  Copyright © 2018 burr2. All rights reserved.
//

#ifndef Complex_tcc
#define Complex_tcc

#include <algorithm>
#include "Complex.hpp"
#include "Interval.hpp"
#include "NumberUtils.hpp"

namespace homotopy
{
    ///////////////////////
    // Forward declaration
    ///////////////////////
    
    inline mpfr_prec_t get_precision(int);
    inline mpfr_prec_t get_precision(unsigned int);
    inline mpfr_prec_t get_precision(float);
    inline mpfr_prec_t get_precision(double);
    
    inline bool set_precision(int, mpfr_prec_t);
    inline bool set_precision(unsigned int, mpfr_prec_t);
    inline bool set_precision(float, mpfr_prec_t);
    inline bool set_precision(double, mpfr_prec_t);

    inline bool increase_precision(int);
    inline bool increase_precision(unsigned int);
    inline bool increase_precision(float);
    inline bool increase_precision(double);
    
    inline bool increase_precision(int, mpfr_prec_t);
    inline bool increase_precision(unsigned int, mpfr_prec_t);
    inline bool increase_precision(float, mpfr_prec_t);
    inline bool increase_precision(double, mpfr_prec_t); 

    inline bool decrease_precision(int);
    inline bool decrease_precision(unsigned int);
    inline bool decrease_precision(float);
    inline bool decrease_precision(double);
    
    inline bool decrease_precision(int, mpfr_prec_t);
    inline bool decrease_precision(unsigned int, mpfr_prec_t);
    inline bool decrease_precision(float, mpfr_prec_t);
    inline bool decrease_precision(double, mpfr_prec_t);

    bool is_interval(const Interval &);

    bool is_interval(const MPFR &);

    bool is_interval(const MPFI &); 

    bool is_interval(double ); 

    template <typename Number>
    bool is_interval(const Number &);

    MPFR square_norm(const MPFR &);
    bool square_norm(MPFR &, const MPFR &);

    double square_norm(double);
    bool square_norm(double &, double );
    
    
    ///////////////////////
    // Declaring static members
    ///////////////////////
    template <typename Type> Type Complex<Type>::temp1;
    template <typename Type> Type Complex<Type>::temp2;

    /*
    template <typename Number>
    Number square_norm(const Number &);
    template <typename Number>
    bool square_norm(Number &, const Number &);
    */
    ///////////////////////
    // Constructors
    ///////////////////////
    
    // Default constructor.  Sets everything equal to 0 using default precision.
    template <typename Type>
    Complex<Type>::Complex() : real(0), imag(0)
    {
    }
    
    // Constructor with precision.  Sets complex number equal to 0 using input precision.
    // Does not use precision in constructors due to non-uniform behavior
    template <typename Type>
    Complex<Type>::Complex(mpfr_prec_t prec) : real(0), imag(0)
    {
        set_precision(real,prec);
        set_precision(imag,prec);
    }
    
    // Copy constructor.  Copies values in input complex number
    // Uses copy constructor of Type
    template <typename Type>
    template <typename Type1>
    Complex<Type>::Complex(const Complex<Type1> &complex)
    {
        convert(real,complex.real);
        convert(imag,complex.imag);
    }
    
    // Copy constructor with precision.  Copies values in input complex number
    // Uses copy constructor of Type
    template <typename Type>
    template <typename Type1>
    Complex<Type>::Complex(const Complex<Type1> & complex, mpfr_prec_t prec)
    {
        convert(real,complex.real);
        convert(imag,complex.imag);
        set_precision(real,prec);
        set_precision(imag,prec);
    }
    
    // Sets complex number equal to input (as a real number)
    template <typename Type>
    Complex<Type>::Complex(int value) : real(value),imag(0.0)
    {
    }
    
    template <typename Type>
    Complex<Type>::Complex(unsigned int value) : real(value),imag(0.0)
    {
    }
    
    template <typename Type>
    Complex<Type>::Complex(float value) : real(value),imag(0.0)
    {
    }
    
    template <typename Type>
    Complex<Type>::Complex(double value) : real(value),imag(0.0)
    {
    }

    // Templated constructor.
    // Assumes the input can be typecast to a double
    template <typename Type>
    template <typename Number>
    Complex<Type>::Complex(const Number &value):imag(0.0)
    {
        convert(real,value);
    }
    
    // Sets complex number equal to input (as a real number)
    // Uses user defined precision
    template <typename Type>
    Complex<Type>::Complex(int value, mpfr_prec_t prec) : real(value),imag(0.0)
    {
        set_precision(real,prec);
        set_precision(imag,prec);
    }
    
    template <typename Type>
    Complex<Type>::Complex(unsigned int value, mpfr_prec_t prec) : real(value),imag(0.0)
    {
        set_precision(real,prec);
        set_precision(imag,prec);
    }
    
    template <typename Type>
    Complex<Type>::Complex(float value, mpfr_prec_t prec) : real(value),imag(0.0)
    {
        set_precision(real,prec);
        set_precision(imag,prec);
    }
    
    template <typename Type>
    Complex<Type>::Complex(double value, mpfr_prec_t prec) : real(value),imag(0.0)
    {
        set_precision(real,prec);
        set_precision(imag,prec);
    }
    
    // Templated constructor.
    // Assumes the input can be typecast to a double
    template <typename Type>
    template <typename Number>
    Complex<Type>::Complex(const Number &value, mpfr_prec_t prec) : real(value), imag(0.0)
    {
        set_precision(real,prec);
        set_precision(imag,prec);
    }
    
    // Parameter constructors.
    // Sets complex number equal to input as a pair
    // Should be slightly faster than pass by reference
    template <typename Type>
    Complex<Type>::Complex(int Real, int Imag) : real(Real), imag(Imag)
    {
    }
    
    template <typename Type>
    Complex<Type>::Complex(unsigned int Real, unsigned int Imag) : real(Real), imag(Imag)
    {
    }
    
    template <typename Type>
    Complex<Type>::Complex(float Real, float Imag) : real(Real), imag(Imag)
    {
    }
    
    template <typename Type>
    Complex<Type>::Complex(double Real, double Imag) : real(Real), imag(Imag)
    {
    }
    
    // Templated constructor.
    // Sets complex number equal to input as a pair
    // Assumes the input can be typecast to a double
    template <typename Type>
    template <typename Number1, typename Number2>
    Complex<Type>::Complex(const Number1 & Real,const Number2 & Imag) : real(Real),imag(Imag)
    {
    }
    
    // Sets complex number equal to input as a pair
    // Should be slightly faster than pass by reference
    // Uses user defined precision
    template <typename Type>
    Complex<Type>::Complex(int Real, int Imag, mpfr_prec_t prec) : real(Real), imag(Imag)
    {
        set_precision(real,prec);
        set_precision(imag,prec);
    }
    
    template <typename Type>
    Complex<Type>::Complex(unsigned int Real, unsigned int Imag, mpfr_prec_t prec) : real(Real), imag(Imag)
    {
        set_precision(real,prec);
        set_precision(imag,prec);
    }
    
    template <typename Type>
    Complex<Type>::Complex(float Real, float Imag, mpfr_prec_t prec) : real(Real), imag(Imag)
    {
        set_precision(real,prec);
        set_precision(imag,prec);
    }
    
    template <typename Type>
    Complex<Type>::Complex(double Real, double Imag, mpfr_prec_t prec) : real(Real), imag(Imag)
    {
        set_precision(real,prec);
        set_precision(imag,prec);
    }
    
    // Templated constructor.
    // Sets complex number equal to input as a pair
    // Assumes the input can be typecast to a double
    // Uses user defined precision
    template <typename Type>
    template <typename Number1, typename Number2>
    Complex<Type>::Complex(const Number1 & Real,const Number2 & Imag, mpfr_prec_t prec) : real(Real),imag(Imag)
    {
        set_precision(real,prec);
        set_precision(imag,prec);
    }
    
    ///////////////////////
    // Operators
    ///////////////////////
    
    // Reassign the real part or the imag part
    template <typename Type>
    Complex<Type> & Complex<Type>::set_real(int value)
    {
        real = value;
        return *this;
    }

    template <typename Type>
    Complex<Type> & Complex<Type>::set_real(unsigned int value)
    {
        real = value;
        return *this;
    }

    template <typename Type>
    Complex<Type> & Complex<Type>::set_real(float value)
    {
        real = value;
        return *this;
    }

    template <typename Type>
    Complex<Type> & Complex<Type>::set_real(double value)
    {
        real = value;
        return *this;
    }

    template <typename Type>
    template <typename Type1>
    Complex<Type> & Complex<Type>::set_real(Type1 & value)
    {
        real = value;
        return *this;
    }

    template <typename Type>
    Complex<Type> & Complex<Type>::set_imag(int value)
    {
        imag = value;
        return *this;
    }

    template <typename Type>
    Complex<Type> & Complex<Type>::set_imag(unsigned int value)
    {
        imag = value;
        return *this;
    }

    template <typename Type>
    Complex<Type> & Complex<Type>::set_imag(float value)
    {
        imag = value;
        return *this;
    }

    template <typename Type>
    Complex<Type> & Complex<Type>::set_imag(double value)
    {
        imag = value;
        return *this;
    }

    template <typename Type>
    template <typename Type1>
    Complex<Type> & Complex<Type>::set_imag(Type1 & value)
    {
        imag = value;
        return *this;
    }

    // Assignment operators.  Uses assignment of Type
    template <typename Type>
    template <typename Number>
    Complex<Type> & Complex<Type>::operator=(const Complex<Number> &complex)
    {
        convert(real,complex.real);
        convert(imag,complex.imag);
        return *this;
    }
    
    // Assignment operators to (real) input.
    template <typename Type>
    Complex<Type> & Complex<Type>::operator=(int value)
    {
        real = value;
        imag = 0;
        return *this;
    }
    
    template <typename Type>
    Complex<Type> & Complex<Type>::operator=(unsigned int value)
    {
        real = value;
        imag = 0;
        return *this;
    }
    
    template <typename Type>
    Complex<Type> & Complex<Type>::operator=(float value)
    {
        real = value;
        imag = 0;
        return *this;
    }
    
    template <typename Type>
    Complex<Type> & Complex<Type>::operator=(double value)
    {
        real = value;
        imag = 0;
        return *this;
    }
    
    // Templated assignment operator
    template <typename Type>
    template <typename Number>
    Complex<Type> & Complex<Type>::operator=(const Number &value)
    {
        real = value;
        imag = 0;
        return *this;
    }
    
    // In-place addition operator.
    template <typename Type>
    template <typename Number>
    Complex<Type> & Complex<Type>::operator+=(const Complex<Number> &complex)
    {
        real += complex.real;
        imag += complex.imag;
        return *this;
    }
    
    // Addition operator to (real) input.
    template <typename Type>
    Complex<Type> & Complex<Type>::operator+=(int value)
    {
        real += value;
        return *this;
    }
    
    template <typename Type>
    Complex<Type> & Complex<Type>::operator+=(unsigned int value)
    {
        real += value;
        return *this;
    }
    
    template <typename Type>
    Complex<Type> & Complex<Type>::operator+=(float value)
    {
        real += value;
        return *this;
    }
    
    template <typename Type>
    Complex<Type> & Complex<Type>::operator+=(double value)
    {
        real += value;
        return *this;
    }
    
    // Templated addition operator
    // Assumes that the type of Number can be added to Type.
    template <typename Type>
    template <typename Number>
    Complex<Type> & Complex<Type>::operator+=(const Number &value)
    {
        real += value;
        return *this;
    }
    
    // In-place subtraction operator.
    template <typename Type>
    template <typename Number>
    Complex<Type> & Complex<Type>::operator-=(const Complex<Number> &complex)
    {
        real -= complex.real;
        imag -= complex.imag;
        return *this;
    }
    
    // Subtraction operator to (real) input.
    template <typename Type>
    Complex<Type> & Complex<Type>::operator-=(int value)
    {
        real -= value;
        return *this;
    }
    
    template <typename Type>
    Complex<Type> & Complex<Type>::operator-=(unsigned int value)
    {
        real -= value;
        return *this;
    }
    
    template <typename Type>
    Complex<Type> & Complex<Type>::operator-=(float value)
    {
        real -= value;
        return *this;
    }
    
    template <typename Type>
    Complex<Type> & Complex<Type>::operator-=(double value)
    {
        real -= value;
        return *this;
    }
    
    // Templated subtraction operator
    // Assumes that the type of Number can be added to Type.
    template <typename Type>
    template <typename Number>
    Complex<Type> & Complex<Type>::operator-=(const Number &value)
    {
        real -= value;
        return *this;
    }
    
    // In-place multiplication operator.
    template <typename Type>
    template <typename Number>
    Complex<Type> & Complex<Type>::operator*=(const Complex<Number> &complex)
    {
        temp1 = real;
        temp1 *= complex.real;
        temp2 = imag;
        temp2 *= complex.imag;
        temp1 -= temp2;
        
        if(&complex == this)
        {
            imag *= real;
            imag *= 2;
        }
        else
        {
            imag *= complex.real;
            temp2 = real;
            temp2 *= complex.imag;
            imag += temp2;
        }
        
        real = temp1;
        return *this;
    }
    
    // Multiplication operator to (real) input.
    template <typename Type>
    Complex<Type> & Complex<Type>::operator*=(int value)
    {
        real *= value;
        imag *= value;
        return *this;
    }
    
    template <typename Type>
    Complex<Type> & Complex<Type>::operator*=(unsigned int value)
    {
        real *= value;
        imag *= value;
        return *this;
    }
    
    template <typename Type>
    Complex<Type> & Complex<Type>::operator*=(float value)
    {
        real *= value;
        imag *= value;
        return *this;
    }
    
    template <typename Type>
    Complex<Type> & Complex<Type>::operator*=(double value)
    {
        real *= value;
        imag *= value;
        return *this;
    }
    
    // Templated multiplication operator
    // Assumes that the type of Number can be added to Type.
    template <typename Type>
    template <typename Number>
    Complex<Type> & Complex<Type>::operator*=(const Number &value)
    {
        real *= value;
        imag *= value;
        return *this;
    }
    
    // In-place division operator.
    template <typename Type>
    template <typename Number>
    Complex<Type> & Complex<Type>::operator/=(const Complex<Number> &complex)
    {
        Type newNorm, newReal;
        newReal = real;
        newReal *= complex.real;
        newNorm = imag;
        newNorm *= complex.imag;
        newReal += newNorm;
        norm(newNorm,complex);
        
        if(&complex == this)
        {
            imag *= real;
            imag -= imag;
        }
        else
        {
            imag *= complex.real;
            real *= complex.imag;
            imag -= real;
        }
        
        real = newReal;
        real /= newNorm;
        imag /= newNorm;
        return *this;
    }
    
    // Division operator to (real) input.
    template <typename Type>
    Complex<Type> & Complex<Type>::operator/=(int value)
    {
        real /= value;
        imag /= value;
        return *this;
    }
    
    template <typename Type>
    Complex<Type> & Complex<Type>::operator/=(unsigned int value)
    {
        real /= value;
        imag /= value;
        return *this;
    }
    
    template <typename Type>
    Complex<Type> & Complex<Type>::operator/=(float value)
    {
        real /= value;
        imag /= value;
        return *this;
    }
    
    template <typename Type>
    Complex<Type> & Complex<Type>::operator/=(double value)
    {
        real /= value;
        imag /= value;
        return *this;
    }
    
    // Templated division operator
    // Assumes that the type of Number can be added to Type.
    template <typename Type>
    template <typename Number>
    Complex<Type> & Complex<Type>::operator/=(const Number &value)
    {
        real /= value;
        imag /= value;
        return *this;
    }
    
    ///////////////////////
    // Output operator
    ///////////////////////
    
    // Display function for an MPFI.  Begins by converting MPFI to a double or pair of doubles.
    // If higher output precision is needed, then we will need a new function.
    template <typename Number>
    inline std::ostream& operator<<(std::ostream & out, const Complex<Number> & complex)
    {
        if(complex.imag == 0)
            return out << complex.real;
        if(complex.real == 0)
        {
            if (complex.imag == 1)
                return out << "i";
            if (complex.imag == -1)
                return out << "-i";
            return out << complex.imag << "i";
        }
        if(complex.imag>0)
        {
            if (complex.imag == 1)
                return out << complex.real << "+" << "i";
            return out << complex.real << "+" << complex.imag << "i";
        }
        if (complex.imag == -1)
            return out << complex.real << "+" << complex.imag << "i";
        return out << complex.real << complex.imag << "i";
    }
    
    ///////////////////////
    // Precision operators
    ///////////////////////
    
    // Returns the current precision.
    template <typename Number>
    inline mpfr_prec_t get_precision(const Complex<Number> & complex)
    {
        return std::max<mpfr_prec_t>(get_precision(complex.real),get_precision(complex.imag));
    }
    
    // Sets the precision of the complex (without destroying data).
    template <typename Number>
    inline bool set_precision(Complex<Number> & complex, mpfr_prec_t prec)
    {
        return (set_precision(complex.real,prec)) && (set_precision(complex.imag,prec));
    }
    
    // Increases the precision to the input precision.
    // If there is only one input, the precision is doubled
    template <typename Number>
    inline bool increase_precision(Complex<Number> &complex)
    {
        return (increase_precision(complex.real)) && (increase_precision(complex.imag));
    }
    
    template <typename Number>
    inline bool increase_precision(Complex<Number> &complex, mpfr_prec_t prec)
    {
        return (increase_precision(complex.real,prec)) && (increase_precision(complex.imag,prec));
    }
    
    // Decreases the precision to the input precision.
    // If there is only one input, the precision is halved
    template <typename Number>
    inline bool decrease_precision(Complex<Number> &complex)
    {
        return (decrease_precision(complex.real)) && (decrease_precision(complex.imag));
    }
    
    template <typename Number>
    inline bool decrease_precision(Complex<Number> &complex, mpfr_prec_t prec)
    {
        return (decrease_precision(complex.real,prec)) && (decrease_precision(complex.imag,prec));
    }
    
    // Get exponent of largest element.  An approximation to the exponent of norm.
    template <typename Number>
    inline long exponent(Complex<Number> &complex)
    {
        return max_exponent(complex);
    }
    
    // Versions for largest and smallest.
    // max_exponent is the same as exponent.
    template <typename Number>
    inline long max_exponent(Complex<Number> &complex)
    {
        return std::max(exponent(complex.real),exponent(complex.imag));
    }
    
    template <typename Number>
    inline long min_exponent(Complex<Number> &complex)
    {
        return std::min(exponent(complex.real),exponent(complex.imag));
    }
    
    ///////////////////////
    // Operations
    ///////////////////////
    
    // Returns or sets the square of the norm of a complex number
    template <typename Number>
    inline Number norm(const Complex<Number> &complex)
    {
        return square_norm(complex.real) + square_norm(complex.imag);
    }
    
    template <typename Number>
    inline bool norm(Number &number, const Complex<Number> &complex)
    {
        Number value;
        square_norm(value,complex.real);
        square_norm(number,complex.imag);
        number += value;
        return true;
    }
    
    // Returns or sets the real part of a complex number
    template <typename Number>
    inline Number real(const Complex<Number> &complex)
    {
        return complex.real;
    }
    
    template <typename Number>
    inline bool real(Number &number, const Complex<Number> &complex)
    {
        number = complex.real;
        return true;
    }
    
    // Returns or sets the imaginary part of a complex number
    template <typename Number>
    inline Number imag(const Complex<Number> &complex)
    {
        return complex.imag;
    }
    
    template <typename Number>
    inline bool imag(Number &number,const Complex<Number> &complex)
    {
        number = complex.imag;
        return true;
    }
    
    // Returns or sets exponential of a complex number.
    // The power must be an unsigned integer.
    template <typename Number>
    inline Complex<Number> exp(Complex<Number> complex,unsigned int power)
    {
        Complex<Number> base(1,0);
        while (power != 0)
        {
            if ((power % 2)==1)
                base *= complex;
            power /= 2;
            complex *= complex;
        }
        return base;
    }
    
    template <typename Number>
    inline bool exp(Complex<Number> &complex,const Complex<Number> &value, unsigned int power)
    {
        Complex<Number> base = value;
        complex = 1;
        while (power != 0)
        {
            if ((power % 2)==1)
                complex *= base;
            power /= 2;
            base *= base;
        }
        return true;
    }
    
    ///////////////////////
    // Operators
    ///////////////////////
    
    // Addition operator.  Assumes that the second type can be added to the first type.
    template <typename Number1, typename Number2>
    inline Complex<Number1> operator+(Complex<Number1> complex1, const Complex<Number2> &complex2)
    {
        complex1 += complex2;
        return complex1;
    }
    
    // Addition operator.
    template <typename Number>
    inline Complex<Number> operator+(Complex<Number> complex,int value)
    {
        complex += value;
        return complex;
    }
    
    template <typename Number>
    inline Complex<Number> operator+(int value, Complex<Number> complex)
    {
        complex += value;
        return complex;
    }
    
    template <typename Number>
    inline Complex<Number> operator+(Complex<Number> complex,unsigned int value)
    {
        complex += value;
        return complex;
    }
    
    template <typename Number>
    inline Complex<Number> operator+(unsigned int value, Complex<Number> complex)
    {
        complex += value;
        return complex;
    }
    
    template <typename Number>
    inline Complex<Number> operator+(Complex<Number> complex, float value)
    {
        complex += value;
        return complex;
    }
    
    template <typename Number>
    inline Complex<Number> operator+(float value, Complex<Number> complex)
    {
        complex += value;
        return complex;
    }
    
    template <typename Number>
    inline Complex<Number> operator+(Complex<Number> complex, double value)
    {
        complex += value;
        return complex;
    }
    
    template <typename Number>
    inline Complex<Number> operator+(double value, Complex<Number> complex)
    {
        complex += value;
        return complex;
    }
    
    // Templated addition operator
    // Assumes the second type can be added to the first.
    template <typename Number1, typename Number2>
    inline Complex<Number1> operator+(const Number2 &value,Complex<Number1> complex)
    {
        complex += value;
        return complex;
    }
    
    template <typename Number1, typename Number2>
    inline Complex<Number1> operator+(Complex<Number1> complex,const Number2 &value)
    {
        complex += value;
        return complex;
    }
    
    // Subtraction operator.  Assumes that the second type can be added to the first type.
    template <typename Number1, typename Number2>
    inline Complex<Number1> operator-(Complex<Number1> complex1, const Complex<Number2> &complex2)
    {
        complex1 -= complex2;
        return complex1;
    }
    
    // Subtraction operator.
    template <typename Number>
    inline Complex<Number> operator-(Complex<Number> complex,int value)
    {
        complex -= value;
        return complex;
    }
    
    template <typename Number>
    inline Complex<Number> operator-(int value, Complex<Number> complex)
    {
        complex *= -1;
        complex += value;
        return complex;
    }
    
    template <typename Number>
    inline Complex<Number> operator-(Complex<Number> complex,unsigned int value)
    {
        complex -= value;
        return complex;
    }
    
    template <typename Number>
    inline Complex<Number> operator-(unsigned int value, Complex<Number> complex)
    {
        complex *= -1;
        complex -= value;
        return complex;
    }
    
    template <typename Number>
    inline Complex<Number> operator-(Complex<Number> complex, float value)
    {
        complex -= value;
        return complex;
    }
    
    template <typename Number>
    inline Complex<Number> operator-(float value, Complex<Number> complex)
    {
        complex *= -1;
        complex += value;
        return complex;
    }
    
    template <typename Number>
    inline Complex<Number> operator-(Complex<Number> complex, double value)
    {
        complex -= value;
        return complex;
    }
    
    template <typename Number>
    inline Complex<Number> operator-(double value, Complex<Number> complex)
    {
        complex *= -1;
        complex -= value;
        return complex;
    }
    
    // Templated subtraction operator
    // Assumes the second type can be subtracted from the first.
    template <typename Number1, typename Number2>
    inline Complex<Number1> operator-(const Number2 &value,Complex<Number1> complex)
    {
        complex -= value;
        return complex;
    }
    
    template <typename Number1, typename Number2>
    inline Complex<Number1> operator-(Complex<Number1> complex,const Number2 &value)
    {
        complex *= -1;
        complex += value;
        return complex;
    }

    // Multiplication operator.  Assumes that the second type can be multiplied by the first type.
    template <typename Number1, typename Number2>
    inline Complex<Number1> operator*(Complex<Number1> complex1, const Complex<Number2> &complex2)
    {
        complex1 *= complex2;
        return complex1;
    }

    // Multiplication operator.
    template <typename Number>
    inline Complex<Number> operator*(Complex<Number> complex,int value)
    {
        complex *= value;
        return complex;
    }
    
    template <typename Number>
    inline Complex<Number> operator*(int value, Complex<Number> complex)
    {
        complex *= value;
        return complex;
    }
    
    template <typename Number>
    inline Complex<Number> operator*(Complex<Number> complex,unsigned int value)
    {
        complex *= value;
        return complex;
    }
    
    template <typename Number>
    inline Complex<Number> operator*(unsigned int value, Complex<Number> complex)
    {
        complex *= value;
        return complex;
    }
    
    template <typename Number>
    inline Complex<Number> operator*(Complex<Number> complex, float value)
    {
        complex *= value;
        return complex;
    }
    
    template <typename Number>
    inline Complex<Number> operator*(float value, Complex<Number> complex)
    {
        complex *= value;
        return complex;
    }
    
    template <typename Number>
    inline Complex<Number> operator*(Complex<Number> complex, double value)
    {
        complex *= value;
        return complex;
    }
    
    template <typename Number>
    inline Complex<Number> operator*(double value, Complex<Number> complex)
    {
        complex *= value;
        return complex;
    }
    
    // Templated multiplication operator
    // Assumes the second type can be multiplied by the first.
    template <typename Number1, typename Number2>
    inline Complex<Number1> operator*(const Number2 &value,Complex<Number1> complex)
    {
        complex *= value;
        return complex;
    }
    
    template <typename Number1, typename Number2>
    inline Complex<Number1> operator*(Complex<Number1> complex,const Number2 &value)
    {
        complex *= value;
        return complex;
    }
    
    // Division operator.  Assumes that the second type can be multiplied and divided by the first type.
    template <typename Number1, typename Number2>
    inline Complex<Number1> operator/(Complex<Number1> complex1, const Complex<Number2> &complex2)
    {
        complex1 /= complex2;
        return complex1;
    }
    
    // Division operator.
    template <typename Number>
    inline Complex<Number> operator/(Complex<Number> complex,int value)
    {
        complex /= value;
        return complex;
    }
    
    template <typename Number>
    inline Complex<Number> operator/(int value, Complex<Number> complex)
    {
        Number newNorm;
        norm(newNorm,complex);
        complex.imag *= -1;
        complex *= value;
        complex /= newNorm;
        return complex;
    }
    
    template <typename Number>
    inline Complex<Number> operator/(Complex<Number> complex,unsigned int value)
    {
        complex /= value;
        return complex;
    }
    
    template <typename Number>
    inline Complex<Number> operator/(unsigned int value, Complex<Number> complex)
    {
        Number newNorm;
        norm(newNorm,complex);
        complex.imag *= -1;
        complex *= value;
        complex /= newNorm;
        return complex;
    }
    
    template <typename Number>
    inline Complex<Number> operator/(Complex<Number> complex, float value)
    {
        complex /= value;
        return complex;
    }
    
    template <typename Number>
    inline Complex<Number> operator/(float value, Complex<Number> complex)
    {
        Number newNorm;
        norm(newNorm,complex);
        complex.imag *= -1;
        complex *= value;
        complex /= newNorm;
        return complex;
    }
    
    template <typename Number>
    inline Complex<Number> operator/(Complex<Number> complex, double value)
    {
        complex /= value;
        return complex;
    }
    
    template <typename Number>
    inline Complex<Number> operator/(double value, Complex<Number> complex)
    {
        Number newNorm;
        norm(newNorm,complex);
        complex.imag *= -1;
        complex *= value;
        complex /= newNorm;
        return complex;
    }
    
    // Templated multiplication operator
    // Assumes the second type can be multiplied by the first.
    template <typename Number1, typename Number2>
    inline Complex<Number1> operator/(const Number2 &value,Complex<Number1> complex)
    {
        complex /= value;
        return complex;
    }
    
    template <typename Number1, typename Number2>
    inline Complex<Number1> operator/(Complex<Number1> complex,const Number2 &value)
    {
        Number1 newNorm;
        norm(newNorm,complex);
        complex.imag *= -1;
        complex *= value;
        complex /= newNorm;
        return complex;
    }

    ///////////////////////
    // Comparison Operators
    ///////////////////////
    
    // Equality operator.  Assumes the first type can be compared with the second
    template <typename Number1,typename Number2>
    inline bool operator==(const Complex<Number1> & complex1, const Complex<Number2> &complex2)
    {
        return (complex1.real == complex2.real) && (complex1.imag == complex2.imag);
    }
    
    // Equality operator.  Compares to (real) variable
    template <typename Number>
    inline bool operator==(const Complex<Number> & complex,int value)
    {
        return (complex.real==value) && (complex.imag == 0);
    }
    
    template <typename Number>
    inline bool operator==(int value, const Complex<Number> & complex)
    {
        return (complex.real==value) && (complex.imag == 0);
    }
    
    template <typename Number>
    inline bool operator==(const Complex<Number> &complex,unsigned int value)
    {
        return (complex.real==value) && (complex.imag == 0);
    }
    
    template <typename Number>
    inline bool operator==(unsigned int value, const Complex<Number> &complex)
    {
        return (complex.real==value) && (complex.imag == 0);
    }
    
    template <typename Number>
    inline bool operator==(const Complex<Number> & complex,float value)
    {
        return (complex.real==value) && (complex.imag == 0);
    }
    
    template <typename Number>
    inline bool operator==(float value, const Complex<Number> &complex)
    {
        return (complex.real==value) && (complex.imag == 0);
    }
    
    template <typename Number>
    inline bool operator==(const Complex<Number> &complex, double value)
    {
        return (complex.real==value) && (complex.imag == 0);
    }
    
    template <typename Number>
    inline bool operator==(double value, const Complex<Number> &complex)
    {
        return (complex.real==value) && (complex.imag == 0);
    }
    
    // Equality operator.
    // Assumes the input types can be compared.
    template <typename Number1, typename Number2>
    inline bool operator==(const Complex<Number1> &complex,const Number2 &value)
    {
        return (complex.real==value) && (complex.imag == 0);
    }
    
    template <typename Number1, typename Number2>
    inline bool operator==(Number2 value, const Complex<Number1> &complex)
    {
        return (complex.real==value) && (complex.imag == 0);
    }

    // Inequality operator.  Assumes the first type can be compared with the second
    template <typename Number1,typename Number2>
    inline bool operator!=(const Complex<Number1> & complex1, const Complex<Number2> &complex2)
    {
        return !(complex1==complex2);
    }
    
    // Inequality operator.  Compares to (real) variable
    template <typename Number>
    inline bool operator!=(const Complex<Number> & complex,int value)
    {
        return !(complex==value);
    }
    
    template <typename Number>
    inline bool operator!=(int value, const Complex<Number> & complex)
    {
        return !(complex==value);
    }
    
    template <typename Number>
    inline bool operator!=(const Complex<Number> &complex,unsigned int value)
    {
        return !(complex==value);
    }
    
    template <typename Number>
    inline bool operator!=(unsigned int value, const Complex<Number> &complex)
    {
        return !(complex==value);
    }
    
    template <typename Number>
    inline bool operator!=(const Complex<Number> & complex,float value)
    {
        return !(complex==value);
    }
    
    template <typename Number>
    inline bool operator!=(float value, const Complex<Number> &complex)
    {
        return !(complex==value);
    }
    
    template <typename Number>
    inline bool operator!=(const Complex<Number> &complex, double value)
    {
        return !(complex==value);
    }
    
    template <typename Number>
    inline bool operator!=(double value, const Complex<Number> &complex)
    {
        return !(complex==value);
    }
    
    // Equality operator.
    // Assumes the input types can be compared.
    template <typename Number1, typename Number2>
    inline bool operator!=(const Complex<Number1> &complex,const Number2 &value)
    {
        return !(complex==value);
    }
    
    template <typename Number1, typename Number2>
    inline bool operator!=(Number2 value, const Complex<Number1> &complex)
    {
        return !(complex==value);
    }
}

#endif /* Complex_tcc */
