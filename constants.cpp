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
#include "alphatheory.hpp"
#include <cfenv>


int main(int argc, char *arg[]) {
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
  std::list<std::vector<homotopy::Complex<homotopy::Interval> > > points;
  std::vector<homotopy::Complex<homotopy::Interval> > point;
  std::string line;
  bool readingPolynomial = true;
  bool encounteredLineBreak = false; // Track if we've encountered a line break
  
  
  
  double realPart, imagPart;
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
      homotopy::Complex<double> coeff(realPart, imagPart);
      
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
	homotopy::Interval realP(realPart-pow(10,-20),realPart+pow(10,-20));
	homotopy::Interval imagP(imagPart-pow(10,-20),imagPart+pow(10,-20));
        homotopy::Complex<homotopy::Interval> value(realP, imagP);
        point.push_back(value);
      }
      points.push_back(point);
      point.clear();
    }
    
  }
  
  int n = points.size();
  inputFile.close();

  std::vector<std::vector<std::list<polyTerm> > > jac = computeJacobian(system);
  std::vector<std::vector<homotopy::Complex<homotopy::Interval> > > jacinv;



  // Open a file for writing
  std::ofstream outputFile("output.txt");
  
  // Check if the file is successfully opened
  if (!outputFile.is_open()) {
    std::cerr << "Failed to open the output file." << std::endl;
    return 1;  // Return an error code
  }

  int a =1;
  
  for (std::list<std::vector<homotopy::Complex<homotopy::Interval> > >::iterator pt_iter = points.begin(); pt_iter != points.end(); ++pt_iter) {
    const auto beforeall = clock::now();
    jacinv = evalJacobian(jac, *pt_iter);
    jacinv = matrix_inv(jacinv);
    double beta_value;
    double gamma_value;
    
    std::vector<std::vector<homotopy::Complex<homotopy::Interval> > > evalSys = evaluate(system, *pt_iter);


    beta_value = beta(system, *pt_iter, jacinv);
    gamma_value = gamma(system, *pt_iter, jacinv);

    double alpha;
    alpha = beta_value * gamma_value;
    
    
    if (alpha >(.157671)*(.157671)) {
      ++failCount;
      outputFile << "failcount is " << failCount << "\n";
    }
    
    outputFile << "values for " << a << " are : \n";
    outputFile << "beta^2 : " << beta_value << "\n";
    outputFile << "gamma^2 : " << gamma_value << "\n";
    outputFile << "alpha^2 : " << alpha << "\n";
    const sec durationall = clock::now() - beforeall;
    std::cout<< "all took" << durationall.count() << std::endl;      
    std::cout<< "solution " <<  a++ << std::endl;
    
  }  
  
  outputFile << failCount ;
  
  // Close the file
  outputFile.close();
  
  
  return 0;
}
