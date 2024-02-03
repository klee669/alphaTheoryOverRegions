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


homotopy::Interval zeroInt(0,0);
homotopy::Interval oneInt(1,1);
homotopy::Complex<homotopy::Interval> zeroI(zeroInt,zeroInt);
homotopy::Complex<homotopy::Interval> oneI(oneInt,zeroInt);




double abs(homotopy::Complex<homotopy::Interval> intv) {
  homotopy::Interval realPart(real(intv));
  homotopy::Interval imagPart(imag(intv));

  double realAbs(abs(realPart));
  double imagAbs(abs(imagPart));
  
  double value;
  value = pow(realAbs,2) + pow(imagAbs,2);

  return value;
}

double euc_int_norm(std::vector<std::vector<homotopy::Complex<homotopy::Interval> > > intv_vec) {

  double value = 0;
  int n = intv_vec[0].size();
  
  for (int i = 0; i < n; ++i) {
    value += abs(intv_vec[0][i]);
  }

  return value;
}

double min_euc_int_norm(std::vector<std::vector<homotopy::Complex<homotopy::Interval> > > intv_vec) {

  double value = 0;
  int n = intv_vec[0].size();
  
  for (int i = 0; i < n; ++i) {
    homotopy::Complex<homotopy::Interval> entry = intv_vec[0][i];
    homotopy::Interval realPart(real(entry));
    homotopy::Interval imagPart(imag(entry));
    homotopy::Complex<double> e(std::min(upper(realPart),lower(realPart)),std::min(upper(imagPart),lower(imagPart)));
    
    value += sqrt(norm(e));
  }

  return value;
}


double pointNorm(std::vector<homotopy::Complex<homotopy::Interval> > point) {
  std::vector<std::vector<homotopy::Complex<homotopy::Interval> > > intv_vec;
  intv_vec.push_back(point);
  double value;
  value = 1 + euc_int_norm(intv_vec);
  return value;
}

double minpointNorm(std::vector<homotopy::Complex<homotopy::Interval> > point) {
  std::vector<std::vector<homotopy::Complex<homotopy::Interval> > > intv_vec;
  intv_vec.push_back(point);
  double value;
  value = 1 + min_euc_int_norm(intv_vec);
  return value;
}


double operatorNorm(std::vector<std::vector<homotopy::Complex<homotopy::Interval> > > M) {
  int n = M.size();
  double value = 0;
  
  for (int i = 0; i < n; ++i) {
    std::vector<std::vector<homotopy::Complex<homotopy::Interval> > > row;
    row.push_back(M[i]);
    value += euc_int_norm(row);
  }

  return value;
}
    




class polyTerm {
public:
  polyTerm(homotopy::Complex<double>, std::list<int>);
  homotopy::Complex<homotopy::Interval> evaluate(std::vector<homotopy::Complex<homotopy::Interval> >);
  homotopy::Complex<double> evaluate(std::vector<homotopy::Complex<double> >);
  homotopy::Complex<double> getCoefficient() const;
  const std::list<int>& getExponent() const;
private:
  std::list<int> exponent;
  homotopy::Complex<double> coefficient;
};

polyTerm::polyTerm(homotopy::Complex<double> coeff, std::list<int> expon) {
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
  

double polyNorm(std::list<polyTerm> polynomial) {
  int n = polynomial.size();
  double absVal = 0;
  std::list<polyTerm>::iterator term;
  int total_degree;
  int d = maxDegree(polynomial);
  double dFac = 1;
  for (int i = 0; i < d; ++i) {
    dFac *= d - i;
  }
  double invdFac = 1/dFac;
  homotopy::Complex<double> coeff;

  for (term = polynomial.begin(); term != polynomial.end(); term++) {
    coeff = term -> getCoefficient();
    const std::list<int>& exp_list = term -> getExponent();
    std::list<int>::const_iterator exp_iter = exp_list.begin();
    total_degree = 0;
    double total_degFac = 1; 
    for (exp_iter = exp_list.begin(); exp_iter != exp_list.end(); exp_iter++) {
      double degFac = 1;
      for (int i = 0; i < *exp_iter; ++i) {
	degFac *= (*exp_iter) - i;
      }
      total_degree += *exp_iter;
      total_degFac *= degFac;
    }
    for (int i = 0; i < (d - total_degree); ++i) {
      total_degFac *= (d - total_degree) - i;
    }
    absVal += (pow(real(coeff),2) + pow(imag(coeff),2))*total_degFac;
  }


  return (invdFac*absVal);
}

double sysNorm(std::list<std::list<polyTerm> > system) {
  std::list<std::list<polyTerm> >::iterator sys_iter;
  double value = 0; 

  for (sys_iter = system.begin(); sys_iter != system.end(); sys_iter++) {
    value += polyNorm(*sys_iter);
  }

  return value;
}


std::vector<std::vector<homotopy::Complex<homotopy::Interval> > > delta_matrix(std::list<std::list<polyTerm> > &system, std::vector<homotopy::Complex<homotopy::Interval> > &point, std::vector<int> &degList) {
  int n = system.size();
  double pn = sqrt(pointNorm(point));
  std::vector<std::vector<homotopy::Complex<homotopy::Interval> > > D;

  
  for (int i = 0; i < n; ++i) {
    std::vector<homotopy::Complex<homotopy::Interval> > dRow(n, 0);
    D.push_back(dRow);
  }

  for (int i = 0; i < n; ++i) {
    int deg = degList[i];
    homotopy::Interval real_ientry(sqrt(deg),sqrt(deg));
    homotopy::Complex<homotopy::Interval> ientry(real_ientry,zeroInt);
    D[i][i] = ientry*pow(pn,deg-1);
  }

  return D;
}
