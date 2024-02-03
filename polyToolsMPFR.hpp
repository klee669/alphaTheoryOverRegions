#include <iostream>
#include <list>
#include <complex>
#include <vector>
#include <cmath>
#include <sstream>
#include <fstream>
#include "Complex.hpp"
#include "Interval.hpp"
#include "MPFI.hpp"
#include "MPFR.hpp"


homotopy::MPFI zeroInt(0,0);
homotopy::MPFI oneInt(1,1);
homotopy::Complex<homotopy::MPFI> zeroI(zeroInt,zeroInt);
homotopy::Complex<homotopy::MPFI> oneI(oneInt,zeroInt);




homotopy::MPFR abs(homotopy::Complex<homotopy::MPFI> intv) {
  homotopy::MPFI realPart(real(intv));
  homotopy::MPFI imagPart(imag(intv));

  homotopy::MPFR realAbs(abs(realPart));
  homotopy::MPFR imagAbs(abs(imagPart));
  
  homotopy::MPFR value;
  value = exp(realAbs,2) + exp(imagAbs,2);

  return value;
}

homotopy::MPFR euc_int_norm(std::vector<std::vector<homotopy::Complex<homotopy::MPFI> > > intv_vec) {

  homotopy::MPFR value = 0;
  int n = intv_vec[0].size();
  
  for (int i = 0; i < n; ++i) {
    value += abs(intv_vec[0][i]);
  }

  return value;
}

homotopy::MPFR min_euc_int_norm(std::vector<std::vector<homotopy::Complex<homotopy::MPFI> > > intv_vec) {

  homotopy::MPFR value = 0;
  int n = intv_vec[0].size();
  
  for (int i = 0; i < n; ++i) {
    homotopy::Complex<homotopy::MPFI> entry = intv_vec[0][i];
    homotopy::MPFI realPart(real(entry));
    homotopy::MPFI imagPart(imag(entry));
    //    homotopy::Complex<homotopy::MPFR> e(std::min(upper(realPart),lower(realPart)),std::min(upper(imagPart),lower(imagPart)));
    homotopy::MPFR e(exp(std::min(upper(realPart),lower(realPart)),2)+exp(std::min(upper(imagPart),lower(imagPart)),2));
    
    value += sqrt(e);
  }

  return value;
}


homotopy::MPFR pointNorm(std::vector<homotopy::Complex<homotopy::MPFI> > point) {
  std::vector<std::vector<homotopy::Complex<homotopy::MPFI> > > intv_vec;
  intv_vec.push_back(point);
  homotopy::MPFR value;
  value = 1 + euc_int_norm(intv_vec);
  return value;
}

homotopy::MPFR minpointNorm(std::vector<homotopy::Complex<homotopy::MPFI> > point) {
  std::vector<std::vector<homotopy::Complex<homotopy::MPFI> > > intv_vec;
  intv_vec.push_back(point);
  homotopy::MPFR value;
  value = 1 + min_euc_int_norm(intv_vec);
  return value;
}


homotopy::MPFR operatorNorm(std::vector<std::vector<homotopy::Complex<homotopy::MPFI> > > M) {
  int n = M.size();
  homotopy::MPFR value = 0;
  
  for (int i = 0; i < n; ++i) {
    std::vector<std::vector<homotopy::Complex<homotopy::MPFI> > > row;
    row.push_back(M[i]);
    value += euc_int_norm(row);
  }

  return value;
}
    




class polyTerm {
public:
  polyTerm(homotopy::Complex<homotopy::MPFR>, std::list<int>);
  homotopy::Complex<homotopy::MPFI> evaluate(std::vector<homotopy::Complex<homotopy::MPFI> >);
  homotopy::Complex<homotopy::MPFR> evaluate(std::vector<homotopy::Complex<homotopy::MPFR> >);
  homotopy::Complex<homotopy::MPFR> getCoefficient() const;
  const std::list<int>& getExponent() const;
private:
  std::list<int> exponent;
  homotopy::Complex<homotopy::MPFR> coefficient;
};

polyTerm::polyTerm(homotopy::Complex<homotopy::MPFR> coeff, std::list<int> expon) {
  coefficient = coeff;
  exponent = expon;
}

int maxDegree(std::list<polyTerm> polynomial) {
  std::list<polyTerm>::iterator term;
  int max_degree = 0;

  for (term = polynomial.begin(); term != polynomial.end(); term++) {
    int deg = 0;
    const std::list<int>& exp_list = term -> getExponent();
    std::list<int>::const_iterator exp_iter = exp_list.begin();
    for (exp_iter = exp_list.begin(); exp_iter != exp_list.end(); exp_iter++) {
      deg += *exp_iter;
    }
    if (deg > max_degree) {
      max_degree = deg;
    }
  }
  return max_degree;
}
  

homotopy::MPFR polyNorm(std::list<polyTerm> polynomial) {
  int n = polynomial.size();
  homotopy::MPFR absVal = 0;
  std::list<polyTerm>::iterator term;
  int total_degree;
  int d = maxDegree(polynomial);
  homotopy::MPFR dFac = 1;
  for (int i = 0; i < d; ++i) {
    dFac *= d - i;
  }
  homotopy::MPFR invdFac = 1/dFac;
  homotopy::Complex<homotopy::MPFR> coeff;

  for (term = polynomial.begin(); term != polynomial.end(); term++) {
    coeff = term -> getCoefficient();
    const std::list<int>& exp_list = term -> getExponent();
    std::list<int>::const_iterator exp_iter = exp_list.begin();
    total_degree = 0;
    homotopy::MPFR total_degFac = 1; 
    for (exp_iter = exp_list.begin(); exp_iter != exp_list.end(); exp_iter++) {
      homotopy::MPFR degFac = 1;
      for (int i = 0; i < *exp_iter; ++i) {
	degFac *= (*exp_iter) - i;
      }
      total_degree += *exp_iter;
      total_degFac *= degFac;
    }
    for (int i = 0; i < (d - total_degree); ++i) {
      total_degFac *= (d - total_degree) - i;
    }
    absVal += (exp(real(coeff),2) + exp(imag(coeff),2))*total_degFac;
  }


  return (invdFac*absVal);
}

homotopy::MPFR sysNorm(std::list<std::list<polyTerm> > system) {
  std::list<std::list<polyTerm> >::iterator sys_iter;
  homotopy::MPFR value = 0; 

  for (sys_iter = system.begin(); sys_iter != system.end(); sys_iter++) {
    value += polyNorm(*sys_iter);
  }

  return value;
}


std::vector<std::vector<homotopy::Complex<homotopy::MPFI> > > delta_matrix(std::list<std::list<polyTerm> > &system, std::vector<homotopy::Complex<homotopy::MPFI> > &point, std::vector<int> &degList) {
  int n = system.size();
  homotopy::MPFR pn = sqrt(pointNorm(point));
  std::vector<std::vector<homotopy::Complex<homotopy::MPFI> > > D;

  
  for (int i = 0; i < n; ++i) {
    std::vector<homotopy::Complex<homotopy::MPFI> > dRow(n, 0);
    D.push_back(dRow);
  }

  for (int i = 0; i < n; ++i) {
    int deg = degList[i];
    homotopy::MPFI real_ientry(sqrt(deg),sqrt(deg));
    homotopy::Complex<homotopy::MPFI> ientry(real_ientry,zeroInt);
    homotopy::Complex<homotopy::MPFI> monom(exp(pn,deg-1));
    D[i][i] = ientry*monom;
  }

  return D;
}
