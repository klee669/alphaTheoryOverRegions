//
//  NumberUtils.tcc
//  Homotopy Continuation
//
//  Created by burr2 on 9/28/18.
//  Copyright © 2018 burr2. All rights reserved.
//

#ifndef NumberUtils_tcc
#define NumberUtils_tcc

#include <cmath>
#include "NumberUtils.hpp"

namespace homotopy
{
    //////////////////////////////////////////////////
    // Upper, lower, median, and width for non-intervals.
    // Default is to return the number itself or zero.
    // Other methods should be written for specific cases.
    //////////////////////////////////////////////////
    
    // Versions for built-in types for speed
    inline int upper(int value)
    {
        return value;
    }
    
    inline bool upper(int &value, int input)
    {
        value = input;
        return true;
    }
    
    inline unsigned int upper(unsigned int value)
    {
        return value;
    }
    
    inline bool upper(unsigned int &value, unsigned int input)
    {
        value = input;
        return true;
    }
    
    inline float upper(float value)
    {
        return value;
    }
    
    inline bool upper(float &value, float input)
    {
        value = input;
        return true;
    }
    
    inline double upper(double value)
    {
        return value;
    }
    
    inline bool upper(double &value, double input)
    {
        value = input;
        return true;
    }

    // Templated versions
    template <typename Type>
    inline Type upper(const Type &value)
    {
        return value;
    }

    template <typename Type>
    inline bool upper(Type &value, const Type &input)
    {
        value = input;
        return true;
    }
    
    // Versions for built-in types for speed
    inline int lower(int value)
    {
        return value;
    }
    
    inline bool lower(int &value, int input)
    {
        value = input;
        return true;
    }
    
    inline unsigned int lower(unsigned int value)
    {
        return value;
    }
    
    inline bool lower(unsigned int &value, unsigned int input)
    {
        value = input;
        return true;
    }
    
    inline float lower(float value)
    {
        return value;
    }
    
    inline bool lower(float &value, float input)
    {
        value = input;
        return true;
    }
    
    inline double lower(double value)
    {
        return value;
    }
    
    inline bool lower(double &value, double input)
    {
        value = input;
        return true;
    }
    
    // Templated versions
    template <typename Type>
    inline Type lower(const Type &value)
    {
        return value;
    }
    
    template <typename Type>
    inline bool lower(Type &value, const Type &input)
    {
        value = input;
        return true;
    }
    
    // Versions for built-in types for speed
    inline int median(int value)
    {
        return value;
    }
    
    inline bool median(int &value, int input)
    {
        value = input;
        return true;
    }
    
    inline unsigned int median(unsigned int value)
    {
        return value;
    }
    
    inline bool median(unsigned int &value, unsigned int input)
    {
        value = input;
        return true;
    }
    
    inline float median(float value)
    {
        return value;
    }
    
    inline bool median(float &value, float input)
    {
        value = input;
        return true;
    }
    
    inline double median(double value)
    {
        return value;
    }
    
    inline bool median(double &value, double input)
    {
        value = input;
        return true;
    }
    
    // Templated versions
    template <typename Type>
    inline Type median(const Type &value)
    {
        return value;
    }
    
    template <typename Type>
    inline bool median(Type &value, const Type &input)
    {
        value = input;
        return true;
    }

    // Versions for built-in types for speed
    inline int width(int value)
    {
        return 0;
    }
    
    inline bool width(int &value, int input)
    {
        value = 0;
        return true;
    }
    
    inline unsigned int width(unsigned int value)
    {
        return 0;
    }
    
    inline bool width(unsigned int &value, unsigned int input)
    {
        value = 0;
        return true;
    }
    
    inline float width(float value)
    {
        return 0;
    }
    
    inline bool width(float &value, float input)
    {
        value = 0;
        return true;
    }
    
    inline double width(double value)
    {
        return 0;
    }
    
    inline bool width(double &value, double input)
    {
        value = 0;
        return true;
    }
    
    // Templated versions
    template <typename Type>
    inline Type width(const Type &value)
    {
        return 0;
    }
    
    template <typename Type>
    inline bool width(Type &value, const Type &input)
    {
        value = 0;
        return true;
    }
    
    //////////////////////////////////////////////////
    // Abs, exp, sqrt for non-intervals.
    // Default is to call cmath.
    // Other methods should be written for specific cases.
    //////////////////////////////////////////////////
    
    // Versions for built-in types for speed
    inline int abs(int value)
    {
        return std::abs(value);
    }
    
    inline bool abs(int &value, int input)
    {
        value = std::abs(input);
        return true;
    }
    
    inline unsigned int abs(unsigned int value)
    {
        return value;
    }
    
    inline bool abs(unsigned int &value, unsigned int input)
    {
        value = input;
        return true;
    }
    
    inline float abs(float value)
    {
        return std::abs(value);
    }
    
    inline bool abs(float &value, float input)
    {
        value = std::abs(input);
        return true;
    }
    
    inline double abs(double value)
    {
        return std::abs(value);
    }
    
    inline bool abs(double &value, double input)
    {
        value = std::abs(input);
        return true;
    }
    
    // Templated versions
    template <typename Type>
    inline Type abs(const Type &value)
    {
        return std::abs(value);
    }
    
    template <typename Type>
    inline bool abs(Type &value, const Type &input)
    {
        value = std::abs(input);
        return true;
    }
    
    // Versions for built-in types for speed
    inline int exp(int value, unsigned int power)
    {
        return std::pow(value,power);
    }
    
    inline bool exp(int &value, int input, unsigned int power)
    {
        value = std::pow(input,power);
        return true;
    }
    
    inline unsigned int exp(unsigned int value, unsigned int power)
    {
        return std::pow(value,power);
    }
    
    inline bool exp(unsigned int &value, unsigned int input, unsigned int power)
    {
        value = std::pow(input,power);
        return true;
    }
    
    inline float exp(float value, unsigned int power)
    {
        return std::pow(value,power);
    }
    
    inline bool exp(float &value, float input, unsigned int power)
    {
        value = std::pow(input,power);
        return true;
    }
    
    inline double exp(double value, unsigned int power)
    {
        return std::pow(value,power);
    }
    
    inline bool exp(double &value, double input, unsigned int power)
    {
        value = std::pow(input,power);
        return true;
    }
    
    // Templated versions
    template <typename Type>
    inline Type exp(const Type &value, unsigned int power)
    {
        return std::pow(value,power);
    }
    
    template <typename Type>
    inline bool exp(Type &value, const Type &input, unsigned int power)
    {
        value = std::pow(input,power);
        return true;
    }
    
    // Versions for built-in types for speed
    inline int sqrt(int value)
    {
        return std::sqrt(value);
    }
    
    inline bool sqrt(int &value, int input)
    {
        value = std::sqrt(input);
        return true;
    }
    
    inline unsigned int sqrt(unsigned int value)
    {
        return std::sqrt(value);
    }
    
    inline bool sqrt(unsigned int &value, unsigned int input)
    {
        value = std::sqrt(input);
        return true;
    }
    
    inline float sqrt(float value)
    {
        return std::sqrt(value);
    }
    
    inline bool sqrt(float &value, float input)
    {
        value = std::sqrt(input);
        return true;
    }
    
    inline double sqrt(double value)
    {
        return std::sqrt(value);
    }
    
    inline bool sqrt(double &value, double input)
    {
        value = std::sqrt(input);
        return true;
    }
    
    // Templated versions
    template <typename Type>
    inline Type sqrt(const Type &value)
    {
        return std::sqrt(value);
    }
    
    template <typename Type>
    inline bool sqrt(Type &value, const Type &input)
    {
        value = std::sqrt(input);
        return true;
    }

    ///////////////////////
    // Real, imag, and norm for non-complex numbers
    // Provides the methods for real numbers
    ///////////////////////
    
    // Versions for built-in types for speed
    inline int norm(int value)
    {
        return value*value;
    }
    
    inline bool norm(int &value, int input)
    {
        value = input * input;
        return true;
    }
    
    inline unsigned int norm(unsigned int value)
    {
        return value*value;
    }
    
    inline bool norm(unsigned int &value, unsigned int input)
    {
        value = input * input;
        return true;
    }
    
    inline float norm(float value)
    {
        return value*value;
    }
    
    inline bool norm(float &value, float input)
    {
        value = input * input;
        return true;
    }
    
    inline double norm(double value)
    {
        return value*value;
    }
    
    inline bool norm(double &value, double input)
    {
        value = input * input;
        return true;
    }
    
    // Templated versions
    template <typename Type>
    inline Type norm(const Type &value)
    {
        return value*value;
    }
    
    template <typename Type>
    inline bool norm(Type &value, const Type &input)
    {
        value = input;
        value *= input;
        return true;
    }
    
    // Versions for built-in types for speed
    inline int real(int value)
    {
        return value;
    }
    
    inline bool real(int &value, int input)
    {
        value = input;
        return true;
    }
    
    inline unsigned int real(unsigned int value)
    {
        return value;
    }
    
    inline bool real(unsigned int &value, unsigned int input)
    {
        value = input;
        return true;
    }
    
    inline float real(float value)
    {
        return value;
    }
    
    inline bool real(float &value, float input)
    {
        value = input;
        return true;
    }
    
    inline double real(double value)
    {
        return value;
    }
    
    inline bool real(double &value, double input)
    {
        value = input;
        return true;
    }
    
    // Templated versions
    template <typename Type>
    inline Type real(const Type &value)
    {
        return value;
    }
    
    template <typename Type>
    inline bool real(Type &value, const Type &input)
    {
        value = input;
        return true;
    }
    
    // Versions for built-in types for speed
    inline int imag(int value)
    {
        return 0;
    }
    
    inline bool imag(int &value, int input)
    {
        value = 0;
        return true;
    }
    
    inline unsigned int imag(unsigned int value)
    {
        return 0;
    }
    
    inline bool imag(unsigned int &value, unsigned int input)
    {
        value = 0;
        return true;
    }
    
    inline float imag(float value)
    {
        return 0;
    }
    
    inline bool imag(float &value, float input)
    {
        value = 0;
        return true;
    }
    
    inline double imag(double value)
    {
        return 0;
    }
    
    inline bool imag(double &value, double input)
    {
        value = 0;
        return true;
    }
    
    // Templated versions
    template <typename Type>
    inline Type imag(const Type &value)
    {
        return 0;
    }
    
    template <typename Type>
    inline bool imag(Type &value, const Type &input)
    {
        value = 0;
        return true;
    }

    //////////////////////////////////////////////////
    // Test if a number contains 0 or a given value
    // True if the values are equal
    // Tests if the first input is contained within the second
    // Should be implemented for intervals separately.
    //////////////////////////////////////////////////
    
    // Versions for built-in types for speed
    inline bool contains(int input,int value)
    {
        return input == value;
    }
    
    inline bool contains(int input,unsigned int value)
    {
        return input == value;
    }
    
    inline bool contains(int input,float value)
    {
        return input == value;
    }

    inline bool contains(int input,double value)
    {
        return input == value;
    }
    
    inline bool contains(int input,const MPFR &value)
    {
        return input == value;
    }
    
    inline bool contains(unsigned int input,int value)
    {
        return input == value;
    }

    inline bool contains(unsigned int input,unsigned int value)
    {
        return input == value;
    }

    inline bool contains(unsigned int input,float value)
    {
        return input == value;
    }

    inline bool contains(unsigned int input,double value)
    {
        return input == value;
    }
    
    inline bool contains(unsigned int input,const MPFR &value)
    {
        return input == value;
    }

    inline bool contains(float input,int value)
    {
        return input == value;
    }

    inline bool contains(float input,unsigned int value)
    {
        return input == value;
    }

    inline bool contains(float input,float value)
    {
        return input == value;
    }

    inline bool contains(float input,double value)
    {
        return input == value;
    }
    
    inline bool contains(float input,const MPFR &value)
    {
        return input == value;
    }
    
    inline bool contains(double input,int value)
    {
        return input == value;
    }

    inline bool contains(double input,unsigned int value)
    {
        return input == value;
    }

    inline bool contains(double input,float value)
    {
        return input == value;
    }

    inline bool contains(double input,double value)
    {
        return input == value;
    }
    
    inline bool contains(double input,const MPFR &value)
    {
        return input == value;
    }
    
    inline bool contains(const MPFR &input,int value)
    {
        return input == value;
    }

    inline bool contains(const MPFR &input,unsigned int value)
    {
        return input == value;
    }

    inline bool contains(const MPFR &input,float value)
    {
        return input == value;
    }

    inline bool contains(const MPFR &input,double value)
    {
        return input == value;
    }

    inline bool contains(const MPFR &input,const MPFR & value)
    {
        return input == value;
    }

    // Templated version
    template <typename Type1, typename Type2>
    inline bool contains(const Type1 & input, const Type2 & value)
    {
        return ((input <= upper(value))&&(input >= lower(value)));
    }
    
    // Versions for built-in types for speed
    inline bool containsZer0(int value)
    {
        return (value == 0);
    }
    
    inline bool containsZero(unsigned int value)
    {
        return (value == 0);
    }
    
    inline bool containsZero(float value)
    {
        return (value == 0);
    }
    
    inline bool containsZero(double value)
    {
        return (value == 0);
    }
    
    inline bool containsZero(const MPFR &value)
    {
        return (value == 0);
    }
    
    // Templated version
    template <typename Type>
    inline bool containsZero(const Type &value)
    {
        return ((0 <= upper(value))&&(0 >= lower(value)));
    }
    
    template <typename Type>
    inline bool containsZero(const Complex<Type> &value)
    {
        return (containsZero(value.real) && containsZero(value.imag));
    }
    
    // Returns true if the data has finite endpoints and
    // false if at least one endpoint is infinite
    // Versions for built-in types for speed
    inline bool is_finite(int value)
    {
        return std::isfinite(value);
    }
    
    inline bool is_finite(unsigned int value)
    {
        return std::isfinite(value);
    }
    
    inline bool is_finite(float value)
    {
        return std::isfinite(value);
    }
    
    inline bool is_finite(double value)
    {
        return std::isfinite(value);
    }
    
    // Templated versions
    template <typename Type>
    inline bool is_finite(const Type &value)
    {
        return std::isfinite(value);
    }
    
    template <typename Type>
    inline bool is_finite(const Complex<Type> &value)
    {
        return (is_finite(value.real) && is_finite(value.imag));
    }
    
    //////////////////////////////////////////////////
    // Default precision operators.
    // Always returns false or zero.
    // Precision operators should be implemented for each type that can change precision
    //////////////////////////////////////////////////
    
    // Returns the current precision.
    // Always returns 0
    inline mpfr_prec_t get_precision(int value)
    {
        return 0;
    }
    
    inline mpfr_prec_t get_precision(unsigned int value)
    {
        return 0;
    }

    inline mpfr_prec_t get_precision(float value)
    {
        return 0;
    }

    inline mpfr_prec_t get_precision(double value)
    {
        return 0;
    }
    
    // Templated version
    template <typename Type>
    inline mpfr_prec_t get_precision(const Type &value)
    {
        return 0;
    }

    // Sets the precision of the input (without destroying data).
    // Doesn't do anything for default fixed types
    inline bool set_precision(int value, mpfr_prec_t prec)
    {
        return false;
    }
    
    inline bool set_precision(unsigned int value, mpfr_prec_t prec)
    {
        return false;
    }

    inline bool set_precision(float value, mpfr_prec_t prec)
    {
        return false;
    }

    inline bool set_precision(double value, mpfr_prec_t prec)
    {
        return false;
    }
    
    // Templated version
    template <typename Type>
    inline bool set_precision(Type &value, mpfr_prec_t prec)
    {
        return false;
    }

    // Increases the precision to the input precision.
    // Doesn't do anything for default fixed types
    // Increases the precision to the input precision.
    // Doesn't do anything for default fixed types
    inline bool increase_precision(int value)
    {
        return false;
    }
    
    inline bool increase_precision(unsigned int value)
    {
        return false;
    }

    inline bool increase_precision(float value)
    {
        return false;
    }

    inline bool increase_precision(double value)
    {
        return false;
    }
    
    inline bool increase_precision(int value, mpfr_prec_t prec)
    {
        return false;
    }
    
    inline bool increase_precision(unsigned int value, mpfr_prec_t prec)
    {
        return false;
    }

    inline bool increase_precision(float value, mpfr_prec_t prec)
    {
        return false;
    }

    inline bool increase_precision(double value, mpfr_prec_t prec)
    {
        return false;
    }
    
    // Templated versions
    template <typename Type>
    inline bool increase_precision(Type &value)
    {
        return false;
    }
    
    template <typename Type>
    inline bool increase_precision(Type &value, mpfr_prec_t prec)
    {
        return false;
    }
    
    // Decreases the precision to the input precision.
    // Doesn't do anything for default fixed types
    inline bool decrease_precision(int value)
    {
        return false;
    }
    
    inline bool decrease_precision(unsigned int value)
    {
        return false;
    }

    inline bool decrease_precision(float value)
    {
        return false;
    }

    inline bool decrease_precision(double value)
    {
        return false;
    }
    
    inline bool decrease_precision(int value, mpfr_prec_t prec)
    {
        return false;
    }
    
    inline bool decrease_precision(unsigned int value, mpfr_prec_t prec)
    {
        return false;
    }

    inline bool decrease_precision(float value, mpfr_prec_t prec)
    {
        return false;
    }

    inline bool decrease_precision(double value, mpfr_prec_t prec)
    {
        return false;
    }
    
    // Templated versions
    template <typename Type>
    inline bool decrease_precision(Type & value)
    {
        return false;
    }
    
    template <typename Type>
    inline bool decrease_precision(Type & value, mpfr_prec_t prec)
    {
        return false;
    }
    
    // Get exponent of input.  Mantissa is assumed to be between 1/2 and 1.
    inline long exponent(float value)
    {
        return (long)ilogb(value);
    }
    
    inline long exponent(double value)
    {
        return (long)ilogb(value);
    }
    
    // Templated versions.
    // Exponent of interval type is the maximum exponent.
    template <typename Type>
    inline long exponent(Type &value)
    {
        return max_exponent(value);
    }
    
    // Maximum or minimum exponent of non-interval type is the exponent.
    template <typename Type>
    inline long min_exponent(Type &value)
    {
        return exponent(value);
    }
    
    template <typename Type>
    inline long max_exponent(Type &value)
    {
        return exponent(value);
    }

    //////////////////////////////////////////////////
    // Functions to change default precision.
    //////////////////////////////////////////////////
    
    // Functions to get the default precision of the specified type
    template <typename Type>
    inline mpfr_prec_t get_default_precision()
    {
        return 0;
    }
    
    template <>
    inline mpfr_prec_t get_default_precision<MPFI>()
    {
        return mpfr_get_default_prec();
    }
    
    template <>
    inline mpfr_prec_t get_default_precision<MPFR>()
    {
        return mpfr_get_default_prec();
    }
    
    // Functions to set the default precision of the specified type
    template <typename Type>
    inline bool set_default_precision(mpfr_prec_t prec)
    {
        return false;
    }
    
    template <>
    inline bool set_default_precision<MPFI>(mpfr_prec_t prec)
    {
        if ((prec<MPFR_PREC_MIN)||(prec>MPFR_PREC_MAX))
            return false;
        if (prec<64)
            prec = 64;
        mpfr_set_default_prec(prec);
        return true;
    }
    
    template <>
    inline bool set_default_precision<MPFR>(mpfr_prec_t prec)
    {
        if ((prec<MPFR_PREC_MIN)||(prec>MPFR_PREC_MAX))
            return false;
        if (prec<64)
            prec = 64;
        mpfr_set_default_prec(prec);
        return true;
    }
    
    // Functions to increase default precision of the specified type
    template <typename Type>
    inline bool increase_default_precision()
    {
        return false;
    }
    
    template <typename Type>
    inline bool increase_default_precision(mpfr_prec_t prec)
    {
        return false;
    }
    
    template <>
    inline bool increase_default_precision<MPFR>()
    {
        return set_default_precision<MPFR>(get_default_precision<MPFR>()*2);
    }
    
    template <>
    inline bool increase_default_precision<MPFR>(mpfr_prec_t prec)
    {
        if(prec>get_default_precision<MPFR>())
            return set_default_precision<MPFR>(prec);
        return false;
    }
    
    template <>
    inline bool increase_default_precision<MPFI>()
    {
        return set_default_precision<MPFI>(get_default_precision<MPFI>()*2);
    }
    
    template <>
    inline bool increase_default_precision<MPFI>(mpfr_prec_t prec)
    {
        if(prec>get_default_precision<MPFI>())
            return set_default_precision<MPFI>(prec);
        return false;
    }
    
    // Functions to decrease default precision of the specified type
    template <typename Type>
    inline bool decrease_default_precision()
    {
        return false;
    }
    
    template <typename Type>
    inline bool decrease_default_precision(mpfr_prec_t prec)
    {
        return false;
    }
    
    template <>
    inline bool decrease_default_precision<MPFR>()
    {
        return set_default_precision<MPFR>(get_default_precision<MPFR>()/2);
    }
    
    template <>
    inline bool decrease_default_precision<MPFR>(mpfr_prec_t prec)
    {
        if(prec<get_default_precision<MPFR>())
            return set_default_precision<MPFR>(prec);
        return false;
    }
    
    template <>
    inline bool decrease_default_precision<MPFI>()
    {
        return set_default_precision<MPFI>(get_default_precision<MPFI>()/2);
    }
    
    template <>
    inline bool decrease_default_precision<MPFI>(mpfr_prec_t prec)
    {
        if(prec<get_default_precision<MPFI>())
            return set_default_precision<MPFI>(prec);
        return false;
    }

    //////////////////////////////////////////////////
    // Interval Tests
    //////////////////////////////////////////////////

    inline bool is_interval(const Interval &value)
    {
        return true;
    }

    inline bool is_interval(const MPFR &value)
    {
        return false;
    }

    inline bool is_interval(const MPFI &value)
    {
        return true;
    }

    inline bool is_interval(double )
    {
        return false;
    }

    template <typename Number>
    inline bool is_interval(const Number &value)
    {
        return false;
    }    

    //////////////////////////////////////////////////
    // Get the square of all the elements of the interval 
    //////////////////////////////////////////////////

    inline double square_norm(double input)
    {
        return input*input;
    }
    inline bool square_norm(double &value, double input)
    {
        value = input;
        value *= input;
        return true;
    }

    //////////////////////////////////////////////////
    // Functions to convert data between specified types.
    //////////////////////////////////////////////////
    
    // TODO: Make conversion functions for all data (including complex)

    // Templated conversion functions
    template <typename Type1, typename Type2>
    inline Type1 convert(const Type2 &value)
    {
        Type1 output;
        convert(output,value);
        return output;
    }
    
    template <typename Type1, typename Type2>
    inline bool convert(Complex<Type1> &value, const Type2 &input)
    {
        convert(value.real,input);
        value.imag = 0;
        return true;
    }
    
    template <typename Type1, typename Type2>
    inline Complex<Type1> convert(const Complex<Type2> &value)
    {
        Complex<Type1> output;
        convert(output, value);
        return output;
    }
    
    template <typename Type1, typename Type2>
    inline bool convert(Complex<Type1> &value, const Complex<Type2> &input)
    {
        convert(value.real,input.real);
        convert(value.imag,input.imag);
        return true;
    }
    
    // Conversions between supported types
    inline bool convert(double &value, double input)
    {
        value = input;
        return true;
    }
    
    inline bool convert(double &value, const Interval &input)
    {
        value = median(input);
        return true;
    }
    
    inline bool convert(double &value, const MPFR &input)
    {
        value = mpfr_get_d(input.point,MPFR_RNDN);
        return true;
    }
    
    inline bool convert(double &value, const MPFI &input)
    {
        value = mpfr_get_d(&(*input.interval).left,MPFR_RNDN);
        value += mpfr_get_d(&(*input.interval).right,MPFR_RNDN);
        value /= 2;
        return true;
    }
    
    inline bool convert(Interval &value, double input)
    {
        value = input;
        return true;
    }
    
    inline bool convert(Interval &value, const Interval &input)
    {
        value = input;
        return true;
    }
    
    inline bool convert(Interval &value, const MPFR &input)
    {
        value.lower = mpfr_get_d(input.point,MPFR_RNDU);
        value.lower *= -1;
        value.upper = mpfr_get_d(input.point,MPFR_RNDU);
        return true;
    }
    
    inline bool convert(Interval &value, const MPFI &input)
    {
        value.lower = mpfr_get_d(&(*input.interval).left,MPFR_RNDU);
        value.lower *= -1;
        value.upper = mpfr_get_d(&(*input.interval).right,MPFR_RNDU);
        return true;
    }
    
    inline bool convert(MPFR &value, double input)
    {
        value = input;
        return true;
    }
    
    inline bool convert(MPFR &value, const Interval &input)
    {
        value = median(input);
        return true;
    }
    
    inline bool convert(MPFR &value, const MPFR &input)
    {
        value = input;
        return true;
    }
    
    inline bool convert(MPFR &value, const MPFI &input)
    {
        median(value,input);
        return true;
    }
    
    inline bool convert(MPFI &value, double input)
    {
        value = input;
        return true;
    }
    
    inline bool convert(MPFI &value, const Interval &input)
    {
        mpfr_set_d(&(*value.interval).left,-input.lower,MPFR_RNDD);
        mpfr_set_d(&(*value.interval).right,input.upper,MPFR_RNDU);
        return true;
    }
    
    inline bool convert(MPFI &value, const MPFR &input)
    {
        value = input;
        return true;
    }
    
    inline bool convert(MPFI &value, const MPFI &input)
    {
        value = input;
        return true;
    }
    
    // Templated rounding functions
    template <typename Type1, typename Type2>
    inline Type1 round_up(const Type2 &value)
    {
        Type1 output;
        round_up(output,value);
        return value;
    }
    
    template <typename Type1, typename Type2>
    inline Type1 round_down(const Type2 &value)
    {
        Type1 output;
        round_down(output,value);
        return value;
    }
    
    // Rounding between supported types
    inline bool round_up(double &value, double input)
    {
        value = input;
        return true;
    }
    
    inline bool round_up(double &value,const MPFR &input)
    {
        value = mpfr_get_d(input.point,MPFR_RNDU);
        return true;
    }
    
    inline bool round_up(MPFR &value, double input)
    {
        mpfr_set_d(value.point,input,MPFR_RNDU);
        return true;
    }
    
    inline bool round_up(MPFR &value, const MPFR &input)
    {
        mpfr_set(value.point,input.point,MPFR_RNDU);
        return true;
    }
    
    inline bool round_down(double &value, double input)
    {
        value = input;
        return true;
    }
    
    inline bool round_down(double &value,const MPFR &input)
    {
        value = mpfr_get_d(input.point,MPFR_RNDD);
        return true;
    }
    
    inline bool round_down(MPFR &value, double input)
    {
        value = mpfr_get_d(value.point,MPFR_RNDD);
        return true;
    }
    
    inline bool round_down(MPFR &value, const MPFR &input)
    {
        mpfr_set(value.point,input.point,MPFR_RNDD);
        return true;
    }
    
    //////////////////////////////////////////////////
    // Functions to perform operations between data types.
    //////////////////////////////////////////////////
    
    // Addition operators
    inline double & operator+=(double &value, const Interval & input)
    {
        value += median(input);
        return value;
    }
    
    inline double & operator+=(double &value, const MPFR &input)
    {
        value += convert<double>(input);
        return value;
    }
    
    inline double & operator+=(double &value, const MPFI &input)
    {
        value += convert<double>(input);
        return value;
    }
    
    // Subtraction operators
    inline double & operator-=(double &value, const Interval & input)
    {
        value -= median(input);
        return value;
    }
    
    inline double & operator-=(double &value, const MPFR &input)
    {
        value -= convert<double>(input);
        return value;
    }
    
    inline double & operator-=(double &value, const MPFI &input)
    {
        value -= convert<double>(input);
        return value;
    }
    
    // Multiplication operators
    inline double & operator*=(double &value, const Interval & input)
    {
        value *= median(input);
        return value;
    }
    
    inline double & operator*=(double &value, const MPFR &input)
    {
        value *= convert<double>(input);
        return value;
    }
    
    inline double & operator*=(double &value, const MPFI &input)
    {
        value *= convert<double>(input);
        return value;
    }
    
    // Division operators
    inline double & operator/=(double &value, const Interval & input)
    {
        value /= median(input);
        return value;
    }
    
    inline double & operator/=(double &value, const MPFR &input)
    {
        value /= convert<double>(input);
        return value;
    }
    
    inline double & operator/=(double &value, const MPFI &input)
    {
        value /= convert<double>(input);
        return value;
    }
    
    // Addition operators
    inline Interval & operator+=(Interval &value, const MPFR &input)
    {
        value += convert<Interval>(input);
        return value;
    }
    
    inline Interval & operator+=(Interval &value, const MPFI &input)
    {
        value += convert<Interval>(input);
        return value;
    }
    
    // Subtraction operators
    inline Interval & operator-=(Interval &value, const MPFR &input)
    {
        value -= convert<Interval>(input);
        return value;
    }
    
    inline Interval & operator-=(Interval &value, const MPFI &input)
    {
        value -= convert<Interval>(input);
        return value;
    }
    
    // Multiplication operators
    inline Interval & operator*=(Interval &value, const MPFR &input)
    {
        value *= convert<Interval>(input);
        return value;
    }
    
    inline Interval & operator*=(Interval &value, const MPFI &input)
    {
        value *= convert<Interval>(input);
        return value;
    }
    
    // Division operators
    inline Interval & operator/=(Interval &value, const MPFR &input)
    {
        value /= convert<Interval>(input);
        return value;
    }
    
    inline Interval & operator/=(Interval &value, const MPFI &input)
    {
        value /= convert<Interval>(input);
        return value;
    }
    //////////////////////////////////////////////////
    // Functions to compare between data types.
    //////////////////////////////////////////////////
    
    inline bool operator==(const MPFR &value,const Interval &input)
    {
        if(width(input)!=0)
            return false;
        return value == upper(input);
    }
    
    inline bool operator==(const Interval &input,const MPFR &value)
    {
        if(width(input)!=0)
            return false;
        return value == upper(input);
    }
    
    inline bool operator!=(const MPFR &value, const Interval &input)
    {
        return !(value == input);
    }
    
    inline bool operator!=(const Interval &input, const MPFR &value)
    {
        return !(value == input);
    }

    //////////////////////////////////////////////////
    // Get the position of the first non-zero bit
    // and the length of the nonzero bit
    //////////////////////////////////////////////////
    inline bool get_position_length(double pos, int len, double value)
    {
        double pos_number = std::floor(std::log2(value));
        if(pos >= 0){
            for(int i=0; i < pos; i++)
                value /= 2;
        }
        else{ 
            for(int i=0; i> pos; i--)
                value *= 2;
        }
        len = 0;
        while(value != 0)
        {
            value -= floor(value);
            value *= 2;
            len++;
        }
        return true;
    }

    inline bool get_position_length(double pos, int len, MPFR value)
    {
        len = 0;
        pos = 0;
        if(value == 0)
            return true;
        MPFR pos_number = floor(log2(value));
        convert(pos, pos_number);
        value -= floor(value);
        if(pos >= 0){
            for(int i=0; i < pos; i++)
                value /= 2;
        }
        else{ 
            for(int i=0; i> pos; i--)
                value *= 2;
        }
        while(value != 0)
        {
            value -= floor(value);
            value *= 2;
            len++;
        }
        return true;
    }

    //////////////////////////////////////////////////
    // Print more digits of the double (to make the code compile)
    //////////////////////////////////////////////////
    inline void print_more_digits(double value, unsigned int length)
    {
        std::cout << value << std::endl;
    }

    inline void print_more_digits(const Interval & value, unsigned int length)
    {
        std::cout << value << std::endl;
    }
}

#endif /* NumberUtils_tcc */
