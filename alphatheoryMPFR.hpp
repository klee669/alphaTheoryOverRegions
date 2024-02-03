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
#include "evaluateMPFR.hpp"
#include <cfenv>


homotopy::MPFR mu(std::list<std::list<polyTerm> > &system, std::vector<homotopy::Complex<homotopy::MPFI> > &point, std::vector<std::vector<homotopy::Complex<homotopy::MPFI> > > &jacobian_inverse, std::vector<int> &degList) {
  homotopy::MPFR sys_norm = sysNorm(system);

  std::vector<std::vector<homotopy::Complex<homotopy::MPFI> > > D = delta_matrix(system, point, degList);
  std::vector<std::vector<homotopy::Complex<homotopy::MPFI> > > mu_matrix;
  mu_matrix = jacobian_inverse;
  matrix_mult(mu_matrix, jacobian_inverse, D);
  homotopy::MPFR oper_norm = operatorNorm(mu_matrix);
  homotopy::MPFR value = 0;
  value = sys_norm * oper_norm;
  if (value < 1) {
    value = 1;
  }

  return value;
}

homotopy::MPFR beta(std::list<std::list<polyTerm> > &system, std::vector<homotopy::Complex<homotopy::MPFI> > &point, std::vector<std::vector<homotopy::Complex<homotopy::MPFI> > > &jacobian_inverse) {
  std::vector<std::vector<homotopy::Complex<homotopy::MPFI> > > evalSys = evaluate(system, point);
  std::vector<std::vector<homotopy::Complex<homotopy::MPFI> > > beta_matrix;
  beta_matrix = evalSys;
  matrix_mult(beta_matrix, jacobian_inverse, evalSys);

  homotopy::MPFR value = operatorNorm(beta_matrix);

  return value;
}



homotopy::MPFR gamma(std::list<std::list<polyTerm> > &system, std::vector<homotopy::Complex<homotopy::MPFI> > &point, std::vector<std::vector<homotopy::Complex<homotopy::MPFI> > > &jacobian_inverse) {
  std::vector<int> degList;
  std::list<std::list<polyTerm> >::iterator poly_iter;
  for (poly_iter = system.begin(); poly_iter != system.end(); ++poly_iter) {
    int deg = maxDegree(*poly_iter);
    degList.push_back(deg);
  }
  homotopy::MPFR mu_value = mu(system, point, jacobian_inverse,degList);
  int max_deg = 0;
  homotopy::MPFR pn = minpointNorm(point);
  int n = system.size();

  for (int i = 0; i < n; ++i) {
    int deg = degList[i];
    if (deg > max_deg) {
      max_deg = deg;
    }
  }

  homotopy::MPFR value = (mu_value * pow(max_deg,3))/(4*pn);

  return value;
}



