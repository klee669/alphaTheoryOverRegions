//
//  NumberUtils.hpp
//  Homotopy Continuation
//
//  Created by burr2 on 9/28/18.
//  Copyright © 2018 burr2. All rights reserved.
//

#ifndef NumberUtil_hpp
#define NumberUtil_hpp

#include "Interval.hpp"
#include "MPFR.hpp"
#include "MPFI.hpp"
#include "Complex.hpp"

namespace homotopy
{
    //////////////////////////////////////////////////
    // Forward Definition of Complex class
    //////////////////////////////////////////////////
    template <typename Type>
    class Complex;
    
    //////////////////////////////////////////////////
    // Upper, lower, median, and width for non-intervals.
    // Default is to return the number itself or zero.
    // Other methods should be written for specific cases.
    //////////////////////////////////////////////////
    
    // Versions for built-in types for speed
    inline int upper(int);
    inline bool upper(int &,int);
    
    inline unsigned int upper(unsigned int);
    inline bool upper(unsigned int &,unsigned int);
    
    inline float upper(float);
    inline bool upper(float &, float);
    
    inline double upper(double);
    inline bool upper(double &, double);
    
    // Templated versions
    template <typename Type>
    inline Type upper(const Type &);
    
    template <typename Type>
    inline bool upper(Type &, const Type &);
    
    // Versions for built-in types for speed
    inline int lower(int);
    inline bool lower(int &,int);
    
    inline unsigned int lower(unsigned int);
    inline bool lower(unsigned int &,unsigned int);
    
    inline float lower(float);
    inline bool lower(float &, float);
    
    inline double lower(double);
    inline bool lower(double &, double);
    
    // Templated versions
    template <typename Type>
    inline Type lower(const Type &);
    
    template <typename Type>
    inline bool lower(Type &, const Type &);
    
    // Versions for built-in types for speed
    inline int median(int);
    inline bool median(int &,int);
    
    inline unsigned int median(unsigned int);
    inline bool median(unsigned int &,unsigned int);
    
    inline float median(float);
    inline bool median(float &, float);
    
    inline double median(double);
    inline bool median(double &, double);
    
    // Templated versions
    template <typename Type>
    inline Type median(const Type &);
    
    template <typename Type>
    inline bool median(Type &, const Type &);
    
    // Versions for built-in types for speed
    inline int width(int);
    inline bool width(int &,int);
    
    inline unsigned int width(unsigned int);
    inline bool width(unsigned int &,unsigned int);
    
    inline float width(float);
    inline bool width(float &, float);
    
    inline double width(double);
    inline bool width(double &, double);
    
    // Templated versions
    template <typename Type>
    inline Type width(const Type &);
    
    template <typename Type>
    inline bool width(Type &, const Type &);
    
    //////////////////////////////////////////////////
    // Abs, exp, sqrt for non-intervals.
    // Default is to call cmath.
    // Other methods should be written for specific cases.
    //////////////////////////////////////////////////
    
    // Versions for built-in types for speed
    inline int abs(int);
    inline bool abs(int &,int);
    
    inline unsigned int abs(unsigned int);
    inline bool abs(unsigned int &,unsigned int);
    
    inline float abs(float);
    inline bool abs(float &, float);
    
    inline double abs(double);
    inline bool abs(double &, double);
    
    // Templated versions
    template <typename Type>
    inline Type abs(const Type &);
    
    template <typename Type>
    inline bool abs(Type &, const Type &);
    
    // Versions for built-in types for speed
    inline int exp(int, unsigned int);
    inline bool exp(int &,int, unsigned int);
    
    inline unsigned int exp(unsigned int, unsigned int);
    inline bool exp(unsigned int &,unsigned int, unsigned int);
    
    inline float exp(float, unsigned int);
    inline bool exp(float &, float, unsigned int);
    
    inline double exp(double, unsigned int);
    inline bool exp(double &, double, unsigned int);
    
    // Templated versions
    template <typename Type>
    inline Type exp(const Type &, unsigned int);
    
    template <typename Type>
    inline bool exp(Type &, const Type &, unsigned int);

    // Versions for built-in types for speed
    inline int sqrt(int);
    inline bool sqrt(int &,int);
    
    inline unsigned int sqrt(unsigned int);
    inline bool sqrt(unsigned int &,unsigned int);
    
    inline float sqrt(float);
    inline bool sqrt(float &, float);
    
    inline double sqrt(double);
    inline bool sqrt(double &, double);
    
    // Templated versions
    template <typename Type>
    inline Type sqrt(const Type &);
    
    template <typename Type>
    inline bool sqrt(Type &, const Type &);
    
    ///////////////////////
    // Real, imag, and norm for non-complex numbers
    // Provides the methods for real numbers
    ///////////////////////
    
    // Versions for built-in types for speed
    inline int norm(int);
    inline bool norm(int &,int);
    
    inline unsigned int norm(unsigned int);
    inline bool norm(unsigned int &,unsigned int);
    
    inline float norm(float);
    inline bool norm(float &, float);
    
    inline double norm(double);
    inline bool norm(double &, double);
    
    // Templated versions
    template <typename Type>
    inline Type norm(const Type &);
    
    template <typename Type>
    inline bool norm(Type &, const Type &);
    
    // Versions for built-in types for speed
    inline int real(int);
    inline bool real(int &,int);
    
    inline unsigned int real(unsigned int);
    inline bool real(unsigned int &,unsigned int);
    
    inline float real(float);
    inline bool real(float &, float);
    
    inline double real(double);
    inline bool real(double &, double);
    
    // Templated versions
    template <typename Type>
    inline Type real(const Type &);
    
    template <typename Type>
    inline bool real(Type &, const Type &);
    
    // Versions for built-in types for speed
    inline int imag(int);
    inline bool imag(int &,int);
    
    inline unsigned int imag(unsigned int);
    inline bool imag(unsigned int &,unsigned int);
    
    inline float imag(float);
    inline bool imag(float &, float);
    
    inline double imag(double);
    inline bool imag(double &, double);
    
    // Templated versions
    template <typename Type>
    inline Type imag(const Type &);
    
    template <typename Type>
    inline bool imag(Type &, const Type &);
    
    //////////////////////////////////////////////////
    // Test if a number contains 0 or a given value
    // True if the values are equal
    // Tests if the first input is contained within the second
    // Should be implemented for intervals separately.
    //////////////////////////////////////////////////
    
    // Versions for built-in types for speed
    inline bool contains(int,int);
    inline bool contains(int,unsigned int);
    inline bool contains(int,float);
    inline bool contains(int,double);
    inline bool contains(int,const MPFR &);
    
    inline bool contains(unsigned int,int);
    inline bool contains(unsigned int,unsigned int);
    inline bool contains(unsigned int,float);
    inline bool contains(unsigned int,double);
    inline bool contains(unsigned int,const MPFR &);
    
    inline bool contains(float,int);
    inline bool contains(float,unsigned int);
    inline bool contains(float,float);
    inline bool contains(float,double);
    inline bool contains(float,const MPFR &);
    
    inline bool contains(double,int);
    inline bool contains(double,unsigned int);
    inline bool contains(double,float);
    inline bool contains(double,double);
    inline bool contains(double,const MPFR &);
    
    inline bool contains(const MPFR &,int);
    inline bool contains(const MPFR &,unsigned int);
    inline bool contains(const MPFR &,float);
    inline bool contains(const MPFR &,double);
    inline bool contains(const MPFR &,const MPFR &);
    
    // Templated versions
    template <typename Type1, typename Type2>
    inline bool contains(const Type1 &, const Type2 &);

    // Versions for built-in types for speed
    inline bool containsZero(int);
    inline bool containsZero(unsigned int);
    inline bool containsZero(float);
    inline bool containsZero(double);
    inline bool containsZero(const MPFR &);
    
    // Templated version
    template <typename Type>
    inline bool containsZero(const Type &);
    
    template <typename Type>
    inline bool containsZero(const Complex<Type> &);
    
    // Returns true if the data has finite endpoints and
    // false if at least one endpoint is infinite
    // Versions for built-in types for speed
    inline bool is_finite(int);
    inline bool is_finite(unsigned int);
    inline bool is_finite(float);
    inline bool is_finite(double);
    
    // Templated versions
    template <typename Type>
    inline bool is_finite(const Type &);
    
    template <typename Type>
    inline bool is_finite(const Complex<Type> &);
    
    //////////////////////////////////////////////////
    // Default precision operators.
    // Always returns false or zero.
    // Precision operators should be implemented for each type that can change precision
    //////////////////////////////////////////////////
    
    // Returns the current precision.
    // Always returns 0
    inline mpfr_prec_t get_precision(int);
    inline mpfr_prec_t get_precision(unsigned int);
    inline mpfr_prec_t get_precision(float);
    inline mpfr_prec_t get_precision(double);
    
    // Templated version
    template <typename Type>
    inline mpfr_prec_t get_precision(const Type &);
    
    // Sets the precision of the input (without destroying data).
    // Doesn't do anything for default fixed types
    inline bool set_precision(int, mpfr_prec_t);
    inline bool set_precision(unsigned int, mpfr_prec_t);
    inline bool set_precision(float, mpfr_prec_t);
    inline bool set_precision(double, mpfr_prec_t);
    
    // Templated version
    template <typename Type>
    inline bool set_precision(Type &, mpfr_prec_t);
    
    // Increases the precision to the input precision.
    // Doesn't do anything for default fixed types
    inline bool increase_precision(int);
    inline bool increase_precision(unsigned int);
    inline bool increase_precision(float);
    inline bool increase_precision(double);
    
    inline bool increase_precision(int, mpfr_prec_t);
    inline bool increase_precision(unsigned int, mpfr_prec_t);
    inline bool increase_precision(float, mpfr_prec_t);
    inline bool increase_precision(double, mpfr_prec_t);
    
    // Templated versions
    template <typename Type>
    inline bool increase_precision(Type &);
    
    template <typename Type>
    inline bool increase_precision(Type &, mpfr_prec_t);
    
    // Decreases the precision to the input precision.
    // Doesn't do anything for default fixed types
    inline bool decrease_precision(int);
    inline bool decrease_precision(unsigned int);
    inline bool decrease_precision(float);
    inline bool decrease_precision(double);
    
    inline bool decrease_precision(int, mpfr_prec_t);
    inline bool decrease_precision(unsigned int, mpfr_prec_t);
    inline bool decrease_precision(float, mpfr_prec_t);
    inline bool decrease_precision(double, mpfr_prec_t);
    
    // Templated versions
    template <typename Type>
    inline bool decrease_precision(Type &);
    
    template <typename Type>
    inline bool decrease_precision(Type &, mpfr_prec_t);
    
    // Get exponent of input.  Mantissa is assumed to be between 1/2 and 1.
    inline long exponent(float);
    inline long exponent(double);
    
    // Templated versions.
    // Exponent of interval type is the maximum exponent.
    template <typename Type>
    inline long exponent(Type &);
    
    // Maximum or minimum exponent of non-interval type is the exponent.
    template <typename Type>
    inline long min_exponent(Type &);
    
    template <typename Type>
    inline long max_exponent(Type &);
    
    //////////////////////////////////////////////////
    // Functions to change default precision.
    //////////////////////////////////////////////////
    
    // Functions to get the default precision of the specified type
    template <typename Type>
    inline mpfr_prec_t get_default_precision();
    
    template <>
    inline mpfr_prec_t get_default_precision<MPFI>();
    
    template <>
    inline mpfr_prec_t get_default_precision<MPFR>();
    
    // Functions to get the default precision of the specified type
    template <typename Type>
    inline bool set_default_precision(mpfr_prec_t);
    
    template <>
    inline bool set_default_precision<MPFI>(mpfr_prec_t);
    
    template <>
    inline bool set_default_precision<MPFR>(mpfr_prec_t);
    
    // Functions to increase default precision of the specified type
    template <typename Type>
    inline bool increase_default_precision();
    
    template <typename Type>
    inline bool increase_default_precision(mpfr_prec_t);
    
    template <>
    inline bool increase_default_precision<MPFR>();
    
    template <>
    inline bool increase_default_precision<MPFR>(mpfr_prec_t);
    
    template <>
    inline bool increase_default_precision<MPFI>();
    
    template <>
    inline bool increase_default_precision<MPFI>(mpfr_prec_t);
    
    // Functions to decrease default precision of the specified type
    template <typename Type>
    inline bool decrease_default_precision();
    
    template <typename Type>
    inline bool decrease_default_precision(mpfr_prec_t);
    
    template <>
    inline bool decrease_default_precision<MPFR>();
    
    template <>
    inline bool decrease_default_precision<MPFR>(mpfr_prec_t);
    
    template <>
    inline bool decrease_default_precision<MPFI>();
    
    template <>
    inline bool decrease_default_precision<MPFI>(mpfr_prec_t);

    //////////////////////////////////////////////////
    // Interval Tests
    //////////////////////////////////////////////////

    bool is_interval(const Interval &);

    bool is_interval(const MPFR &);

    bool is_interval(const MPFI &); 

    bool is_interval(double ); 

    template <typename Number>
    bool is_interval(const Number &);

    //////////////////////////////////////////////////
    // Get the square of all the elements of the interval 
    //////////////////////////////////////////////////
    double square_norm(double);
    bool square_norm(double &, double );

    //////////////////////////////////////////////////
    // Functions to convert data between specified types.
    //////////////////////////////////////////////////
    
    // Templated conversion functions
    template <typename Type1, typename Type2>
    inline Type1 convert(const Type2&);
    
    template <typename Type1, typename Type2>
    inline bool convert(Complex<Type1>&, const Type2&);
    
    template <typename Type1, typename Type2>
    inline Complex<Type1> convert(const Complex<Type2> &);
    
    template <typename Type1, typename Type2>
    inline bool convert(Complex<Type1> &, const Complex<Type2> &);
    
    // Conversions between supported types
    inline bool convert(double &, double);
    inline bool convert(double &, const Interval &);
    inline bool convert(double &, const MPFR &);
    inline bool convert(double &, const MPFI &);
    
    inline bool convert(Interval &, double);
    inline bool convert(Interval &, const Interval &);
    inline bool convert(Interval &, const MPFR &);
    inline bool convert(Interval &, const MPFI &);
    
    inline bool convert(MPFR &, double);
    inline bool convert(MPFR &, const Interval &);
    inline bool convert(MPFR &, const MPFR &);
    inline bool convert(MPFR &, const MPFI &);
    
    inline bool convert(MPFI &, double);
    inline bool convert(MPFI &, const Interval &);
    inline bool convert(MPFI &, const MPFR &);
    inline bool convert(MPFI &, const MPFI &);
    
    // Templated rounding functions
    template <typename Type1, typename Type2>
    inline Type1 round_up(const Type2&);
    
    template <typename Type1, typename Type2>
    inline Type1 round_down(const Type2&);
    
    // Rounding between supported types
    // Does not change the precision of the input MPFR.
    inline bool round_up(double &, double);
    inline bool round_up(double &,const MPFR &);
    inline bool round_up(MPFR &, double);
    inline bool round_up(MPFR &, const MPFR &);
    
    inline bool round_down(double &, double);
    inline bool round_down(double &,const MPFR &);
    inline bool round_down(MPFR &, double);
    inline bool round_down(MPFR &, const MPFR &);
    
    //////////////////////////////////////////////////
    // Functions to perform operations between data types.
    //////////////////////////////////////////////////
    
    inline double & operator+=(double &, const Interval &);
    inline double & operator+=(double &, const MPFR &);
    inline double & operator+=(double &, const MPFI &);
    
    inline double & operator-=(double &, const Interval &);
    inline double & operator-=(double &, const MPFR &);
    inline double & operator-=(double &, const MPFI &);
    
    inline double & operator*=(double &, const Interval &);
    inline double & operator*=(double &, const MPFR &);
    inline double & operator*=(double &, const MPFI &);
    
    inline double & operator/=(double &, const Interval &);
    inline double & operator/=(double &, const MPFR &);
    inline double & operator/=(double &, const MPFI &);
    
    inline Interval & operator+=(Interval &, const MPFR &);
    inline Interval & operator+=(Interval &, const MPFI &);
    
    inline Interval & operator-=(Interval &, const MPFR &);
    inline Interval & operator-=(Interval &, const MPFI &);
    
    inline Interval & operator*=(Interval &, const MPFR &);
    inline Interval & operator*=(Interval &, const MPFI &);
    
    inline Interval & operator/=(Interval &, const MPFR &);
    inline Interval & operator/=(Interval &, const MPFI &);
    //////////////////////////////////////////////////
    // Functions to compare between data types.
    //////////////////////////////////////////////////
    
    inline bool operator==(const MPFR &,const Interval &);
    inline bool operator==(const Interval &,const MPFR &);
    
    inline bool operator!=(const MPFR &,const Interval &);
    inline bool operator!=(const Interval &,const MPFR &);

    //////////////////////////////////////////////////
    // Get the position of the first non-zero bit
    // and the length of the nonzero bit
    //////////////////////////////////////////////////
    inline bool get_position_length(double, int , double);

    inline bool get_position_length(double, int , MPFR);

    //////////////////////////////////////////////////
    // Print more digits of the double (to make the code compile)
    //////////////////////////////////////////////////
    inline void print_more_digits(double, unsigned int);
    inline void print_more_digits(const Interval &, unsigned int);
}

// The following is included because the functions are encouraged to be inline or templated.
#include "NumberUtils.tcc"

#endif /* NumberUtil_hpp */
