#include <iostream>
#include <list>
#include <complex>
#include <vector>
#include <cmath>
#include <sstream>
#include <fstream>
#include <chrono>
#include "Complex.hpp"
#include "Interval.hpp"
#include "MPFI.hpp"
#include "MPFR.hpp"
#include "LUdecomp.hpp"
//#include "lutest.hpp"

homotopy::Complex<homotopy::Interval> polyTerm::evaluate(std::vector<homotopy::Complex<homotopy::Interval> > point) {
  homotopy::Complex<homotopy::Interval> value(real(coefficient),imag(coefficient));
  if (coefficient == 0) {
    value = zeroI;
  } else {
  std::list<int>::iterator j = exponent.begin(); // initialize the exponent
  int itern = 0;
  homotopy::Complex<homotopy::Interval> dummy = oneI;

  //  for (std::list<std::complex<double> >::iterator my_iter = point.begin(); my_iter != point.end(); my_iter++, j++) {
  //    value *= pow(*my_iter, *j); // Extract values from iterators and use them as arguments
  //  }
  for(std::vector<homotopy::Complex<homotopy::Interval> >::iterator my_iter = point.begin();my_iter != point.end(); my_iter++) {
      while (itern < *j){
    dummy = dummy*(*my_iter); // Extract values from iterators and use them as arguments
    itern++;
  }
    j++;
    itern = 0;
  }
  value = value*dummy;
  }
  return value;
}



homotopy::Complex<homotopy::Interval> evaluate(std::list<polyTerm> &polynomial, std::vector<homotopy::Complex<homotopy::Interval> > &point) {
  homotopy::Complex<homotopy::Interval> value(0,0);
  std::list<polyTerm>::iterator my_iter;

  for (my_iter = polynomial.begin(); my_iter != polynomial.end(); my_iter++)
    value += my_iter->evaluate(point);

  return value;
}


std::vector<std::vector<homotopy::Complex<homotopy::Interval> > > evaluate(std::list<std::list<polyTerm> > &system, std::vector<homotopy::Complex<homotopy::Interval> > &point) {
  std::vector<std::vector<homotopy::Complex<homotopy::Interval> > > evalSys;
  std::list<std::list<polyTerm> >::iterator polynomial = system.begin();
    std::vector<homotopy::Complex<homotopy::Interval> > evalPoly;

  for (polynomial = system.begin(); polynomial != system.end(); ++polynomial) {
    evalPoly.push_back(evaluate(*polynomial, point));
    evalSys.push_back(evalPoly);
    evalPoly.clear();
  }

  return evalSys;
}





homotopy::Complex<double> polyTerm::getCoefficient() const {
  return coefficient;
}

const std::list<int>& polyTerm::getExponent() const {
  return exponent;
}

std::vector<std::vector<std::list<polyTerm> > > computeJacobian(std::list<std::list<polyTerm> > &system) {
    
  std::list<std::list<polyTerm> >::iterator polynomial = system.begin();
  std::list<polyTerm>::iterator term = polynomial->begin();
  std::vector<std::vector<std::list<polyTerm> > > jacobian;
  std::vector<std::list<polyTerm> > jacRow;
  std::list<std::list<polyTerm> >::iterator my_iter = system.begin();
  homotopy::Complex<double> coeff;

  std::list<polyTerm> jacElement;
  std::list<int> jacElementExp;
  polyTerm jacElementTerm(0,jacElementExp);

  for (polynomial = system.begin(); polynomial != system.end(); polynomial++) {
    int iter_counter = 0;
    for (my_iter = system.begin(); my_iter != system.end(); my_iter++) {
      for (term = polynomial->begin(); term != polynomial->end(); term++) {
        coeff = term ->getCoefficient();
	const std::list<int>& exp_list = term -> getExponent();
	std::list<int>::const_iterator exp = exp_list.begin();
	int var_counter = 0;

      for (exp = exp_list.begin(); exp != exp_list.end(); exp++) {
	if (iter_counter == var_counter) {
	  if (*exp>0) {
	    coeff = homotopy::Complex<double>(*exp,0.0) * coeff;
	    jacElementExp.push_back((*exp)-1);
	    var_counter++;
	  } else {
	    coeff = 0;
	       break;
	  }
	} else {
	  jacElementExp.push_back(*exp);
	  var_counter++;
	}
      }

      jacElementTerm = polyTerm(coeff, jacElementExp);
      
      jacElement.push_back(jacElementTerm);
      jacElementExp.clear();
      }
      iter_counter++;
      jacRow.push_back(jacElement);
      jacElement.clear();

    }
    jacobian.push_back(jacRow);
    jacRow.clear();
  }
  

  
  return jacobian;
}


std::vector<std::vector<homotopy::Complex<homotopy::Interval> > > evalJacobian(std::vector<std::vector<std::list<polyTerm> > > &jac, std::vector<homotopy::Complex<homotopy::Interval> > &point) {

  int n = point.size();
  std::vector<std::vector<homotopy::Complex<homotopy::Interval> > > result;

  // Initialize L and U
  for (int i = 0; i < n; ++i) {
    std::vector<homotopy::Complex<homotopy::Interval> > lRow(n, 0);
    result.push_back(lRow);
  }

  for (int i =0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      result[i][j] = evaluate(jac[i][j],point);
    }
  }




  
  return result;
}





homotopy::Complex<double> polyTerm::evaluate(std::vector<homotopy::Complex<double> > point) {
  homotopy::Complex<double> value(real(coefficient),imag(coefficient));
  std::list<int>::iterator j = exponent.begin(); // initialize the exponent
  int itern = 0;
  homotopy::Complex<double> dummy = 1;

  //  for (std::list<std::complex<double> >::iterator my_iter = point.begin(); my_iter != point.end(); my_iter++, j++) {
  //    value *= pow(*my_iter, *j); // Extract values from iterators and use them as arguments
  //  }
  for(std::vector<homotopy::Complex<double> >::iterator my_iter = point.begin();my_iter != point.end(); my_iter++) {
      while (itern < *j){
    dummy = dummy*(*my_iter); // Extract values from iterators and use them as arguments
    itern++;
  }
    j++;
    itern = 0;
  }
  value = value*dummy;

  return value;
}


homotopy::Complex<double> evaluate(std::list<polyTerm> polynomial, std::vector<homotopy::Complex<double> > point) {
  homotopy::Complex<double> value(0,0);
  std::list<polyTerm>::iterator my_iter;

  for (my_iter = polynomial.begin(); my_iter != polynomial.end(); my_iter++)
    value += my_iter->evaluate(point);

  return value;
}


std::vector<std::vector<homotopy::Complex<double> > > evaluate(std::list<std::list<polyTerm> > system, std::vector<homotopy::Complex<double> > point) {
  std::vector<std::vector<homotopy::Complex<double> > > evalSys;
  std::list<std::list<polyTerm> >::iterator polynomial = system.begin();

  for (polynomial = system.begin(); polynomial != system.end(); ++polynomial) {
    std::vector<homotopy::Complex<double> > evalPoly;
    evalPoly.push_back(evaluate(*polynomial, point));
    evalSys.push_back(evalPoly);
  }

  return evalSys;
}

/*

std::vector<std::vector<homotopy::Complex<double> > > computeJacobian(std::list<std::list<polyTerm> > system, std::vector<homotopy::Complex<double> > point) {
  std::list<std::list<polyTerm> >::iterator polynomial;
  std::vector<std::vector<homotopy::Complex<double> > > jacobian;
  std::vector<homotopy::Complex<double> >::iterator my_iter = point.begin();

  for (polynomial = system.begin(); polynomial != system.end(); polynomial++) {
    int iter_counter = 0;
    std::vector<homotopy::Complex<double> > jacRow;
    for (my_iter = point.begin(); my_iter != point.end(); my_iter++) {
      std::list<polyTerm> jacElement;
      std::list<polyTerm>::iterator term;
      for (term = polynomial->begin(); term != polynomial->end(); term++) {
	homotopy::Complex<double> coeff = term ->getCoefficient();
	const std::list<int>& exp_list = term ->getExponent();
	std::list<int> jacElementExp;
	std::list<int>::const_iterator exp = exp_list.begin();
	int var_counter = 0;
	homotopy::Complex<double> derivative = coeff; 
	
      for (exp = exp_list.begin(); exp != exp_list.end(); exp++) {
	if (iter_counter == var_counter) {
	  if (*exp>0) {
	    coeff = homotopy::Complex<double>(*exp,0.0) * coeff;
	    jacElementExp.push_back((*exp)-1);
	    var_counter++;
	  } else {
	    coeff = 0;
	  }
	} else {
	  jacElementExp.push_back(*exp);
	  var_counter++;
	}
      }
      polyTerm jacElementTerm(coeff, jacElementExp);
	
	jacElement.push_back(jacElementTerm);
      }
      iter_counter++;
      jacRow.push_back(evaluate(jacElement, point));
    }
    jacobian.push_back(jacRow);
  }
  

  
  return jacobian;
}
*/

