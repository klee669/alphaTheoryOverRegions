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
#include "evaluate.hpp"
#include <cfenv>


double mu(std::list<std::list<polyTerm> > &system, std::vector<homotopy::Complex<homotopy::Interval> > &point, std::vector<std::vector<homotopy::Complex<homotopy::Interval> > > &jacobian_inverse, std::vector<int> &degList) {
  double sys_norm = sysNorm(system);

  std::vector<std::vector<homotopy::Complex<homotopy::Interval> > > D = delta_matrix(system, point, degList);
  std::vector<std::vector<homotopy::Complex<homotopy::Interval> > > mu_matrix;
  mu_matrix = jacobian_inverse;
  matrix_mult(mu_matrix, jacobian_inverse, D);
  double oper_norm = operatorNorm(mu_matrix);
  double value = 0;
  value = sys_norm * oper_norm;
  if (value < 1) {
    value = 1;
  }

  return value;
}

double beta(std::list<std::list<polyTerm> > &system, std::vector<homotopy::Complex<homotopy::Interval> > &point, std::vector<std::vector<homotopy::Complex<homotopy::Interval> > > &jacobian_inverse) {
  std::vector<std::vector<homotopy::Complex<homotopy::Interval> > > evalSys = evaluate(system, point);
  std::vector<std::vector<homotopy::Complex<homotopy::Interval> > > beta_matrix;
  beta_matrix = evalSys;
  matrix_mult(beta_matrix, jacobian_inverse, evalSys);

  double value = operatorNorm(beta_matrix);

  return value;
}



double gamma(std::list<std::list<polyTerm> > &system, std::vector<homotopy::Complex<homotopy::Interval> > &point, std::vector<std::vector<homotopy::Complex<homotopy::Interval> > > &jacobian_inverse) {
  std::vector<int> degList;
  std::list<std::list<polyTerm> >::iterator poly_iter;
  for (poly_iter = system.begin(); poly_iter != system.end(); ++poly_iter) {
    int deg = maxDegree(*poly_iter);
    degList.push_back(deg);
  }
  double mu_value = mu(system, point, jacobian_inverse,degList);
  int max_deg = 0;
  double pn = minpointNorm(point);
  int n = system.size();

  for (int i = 0; i < n; ++i) {
    int deg = degList[i];
    if (deg > max_deg) {
      max_deg = deg;
    }
  }

  double value = (mu_value * pow(max_deg,3))/(4*pn);

  return value;
}



double divideFraction(const std::string& fraction) {
    // Split the fraction by '/'
    std::istringstream iss(fraction);
    std::string numeratorStr, denominatorStr;

    if (std::getline(iss, numeratorStr, '/') && std::getline(iss, denominatorStr)) {
        try {
	  double numerator = std::stod(numeratorStr);
	  double denominator = std::stod(denominatorStr);

            if (denominator != 0) {
	      double value = 0;
                return value = numerator / denominator;
            } else {
                std::cerr << "Error: Division by zero." << std::endl;
                return 0.0;  // Return 0.0 for division by zero
            }
        } catch (const std::invalid_argument& e) {
            std::cerr << "Error: Invalid input format." << std::endl;
            return 0.0;  // Return 0.0 for invalid input
        }
    } else {
      double value = std::stod(fraction);
      return value;
    }
}
