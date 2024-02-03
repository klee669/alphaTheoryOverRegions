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
#include "alphatheoryMPFR.hpp"
#include <cfenv>


homotopy::MPFR divideFraction(const std::string& fraction) {
  // Split the fraction by '/'
    std::istringstream iss(fraction);
    std::string numeratorStr, denominatorStr;
    
    if (std::getline(iss, numeratorStr, '/') && std::getline(iss, denominatorStr)) {
      try {
	homotopy::MPFR numerator(numeratorStr);
	homotopy::MPFR denominator(denominatorStr);
	
	if (denominator != 0) {
	  homotopy::MPFR value = 0;
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
      homotopy::MPFR value(fraction);
      return value;
    }
}


int main(int argc, char *arg[]) {
  mpfr_set_default_prec(256);
  
  // Check if a filename is provided as a command-line argument
  if (argc != 2) {
    std::cerr << "Usage: " << arg[0] << " <filename>" << std::endl;
    return 1;  // Return an error code
  }
  
  
  using clock = std::chrono::system_clock;
  using sec = std::chrono::duration<double>;
  char *fileName = NULL;
  fileName = arg[1];
  std::fesetround(FE_UPWARD);
  std::ifstream inputFile(fileName);
  int failCount = 0;
  
  if (!inputFile) {
    std::cerr << "Failed to open input file." << std::endl;
    return 1;
  }
  
  std::list<std::list<polyTerm> > system;
  std::list<polyTerm> polynomial;
  std::list<std::vector<homotopy::Complex<homotopy::MPFI> > > points;
  std::vector<homotopy::Complex<homotopy::MPFI> > point;
  std::string line;
  bool readingPolynomial = true;
  bool encounteredLineBreak = false; // Track if we've encountered a line break
  
  
  
  homotopy::MPFR realPart, imagPart;
  std::list<int> exponents;
  
  
  while (std::getline(inputFile, line)) {
    std::istringstream iss(line);
    
    if (line.empty()) {
      encounteredLineBreak = true;
    }
    
    if (line.find("points") != std::string::npos) {
      readingPolynomial = false;
      continue;
    }
    
    if (readingPolynomial && !encounteredLineBreak) {
      
      std::string str1, str2;
      
      iss >> str1 >> str2;
      
      
      realPart = divideFraction(str1);
      imagPart = divideFraction(str2);
      homotopy::Complex<homotopy::MPFR> coeff(realPart, imagPart);

      int exp;
      while (iss >> exp) {
        exponents.push_back(exp);
      }
      
      
      polyTerm term(coeff, exponents);
      
      polynomial.push_back(term);
      exponents.clear();
    }
    if (readingPolynomial && encounteredLineBreak) {
      system.push_back(polynomial);
      polynomial.clear();
      encounteredLineBreak = false;
      continue;
    }
    if (!readingPolynomial) {
      std::string str1, str2;
      while (iss >> str1 >> str2) {
	realPart = divideFraction(str1);
	imagPart = divideFraction(str2);
	homotopy::MPFI realP(realPart-pow(10,-20),realPart+pow(10,-20));
	homotopy::MPFI imagP(imagPart-pow(10,-20),imagPart+pow(10,-20));
        homotopy::Complex<homotopy::MPFI> value(realP, imagP);
        point.push_back(value);
      }
      points.push_back(point);
      point.clear();
    }
    
  }
  
  int n = points.size();
  inputFile.close();
  
  std::vector<std::vector<std::list<polyTerm> > > jac = computeJacobian(system);
  std::vector<std::vector<homotopy::Complex<homotopy::MPFI> > > jacinv;
  
  
  
  // Open a file for writing
  std::ofstream outputFile("cyclic6/output_cyclic_MPFR_256.txt");
  
  // Check if the file is successfully opened
  if (!outputFile.is_open()) {
    std::cerr << "Failed to open the output file." << std::endl;
    return 1;  // Return an error code
  }
  
  int a =1;
  
  for (std::list<std::vector<homotopy::Complex<homotopy::MPFI> > >::iterator pt_iter = points.begin(); pt_iter != points.end(); ++pt_iter) {
    const auto beforeall = clock::now();
    jacinv = evalJacobian(jac, *pt_iter);
    jacinv = matrix_inv(jacinv);
    homotopy::MPFR beta_value;
    homotopy::MPFR gamma_value;
  
    std::vector<std::vector<homotopy::Complex<homotopy::MPFI> > > evalSys = evaluate(system, *pt_iter);


    beta_value = beta(system, *pt_iter, jacinv);
    gamma_value = gamma(system, *pt_iter, jacinv);

    homotopy::MPFR alpha;
    alpha = sqrt(beta_value * gamma_value);


    if (alpha >.157671) {
      ++failCount;
      outputFile << "failcount is " << failCount << "\n";
    }
    
    outputFile << "values for " << a << " are : \n";
    outputFile << "beta^2 : " << beta_value << "\n";
    outputFile << "gamma^2 : " << gamma_value << "\n";
    outputFile << "alpha^2 : " << beta_value * gamma_value << "\n";
    const sec durationall = clock::now() - beforeall;
    std::cout<< "all took" << durationall.count() << std::endl;  
    
    std::cout<< "solution " <<  a++ << std::endl;
    
  }  
  
  outputFile << failCount ;
  
  // Close the file
  outputFile.close();
  
  
  return 0;
}
