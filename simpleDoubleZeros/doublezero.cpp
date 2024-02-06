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
  mpfr_set_default_prec(128);
  
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

  homotopy::MPFI realV1(.894427,.894427);
  homotopy::MPFI imagV1(0,0);
  homotopy::Complex<homotopy::MPFI> v1(realV1, imagV1);
  homotopy::MPFI realV2(-.447214,-.447214);
  homotopy::MPFI imagV2(0,0);
  homotopy::Complex<homotopy::MPFI> v2(realV2, imagV2);

  
  std::list<std::list<polyTerm> > system;
  std::list<polyTerm> polynomial;
  std::list<std::vector<homotopy::Complex<homotopy::MPFI> > > points;
  std::vector<homotopy::Complex<homotopy::MPFI> > point;
  std::vector<homotopy::Complex<homotopy::MPFI> > rowv;
  std::vector<std::vector<homotopy::Complex<homotopy::MPFI> > >rv;
  std::vector<std::vector<homotopy::Complex<homotopy::MPFI> > >v;
  std::vector<homotopy::Complex<homotopy::MPFI> > vrow1;
  vrow1.push_back(v1);
  std::vector<homotopy::Complex<homotopy::MPFI> > vrow2;
  vrow2.push_back(v2);
  std::string line;
  bool readingPolynomial = true;
  bool encounteredLineBreak = false; // Track if we've encountered a line break
  v.push_back(vrow1);
  v.push_back(vrow2);
  rowv.push_back(v1);
  rowv.push_back(v2);
  rv.push_back(rowv);
  

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
      
      
      //      iss >> realPart >> imagPart; // Read real and imaginary parts
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
	homotopy::MPFI realP(realPart-pow(10,-7),realPart+pow(10,-7));
	homotopy::MPFI imagP(imagPart-pow(10,-7),imagPart+pow(10,-7));
        homotopy::Complex<homotopy::MPFI> value(realP, imagP);
        point.push_back(value);
      }
      points.push_back(point);
      point.clear();
    }
    
  }
  
  int n = points.size();
  inputFile.close();

  std::cout << 1 << " ";


  std::vector<std::vector<std::list<polyTerm> > > jaco = computeJacobian(system);
  std::vector<std::vector<homotopy::Complex<homotopy::MPFI> > > jacinv;
  std::vector<std::vector<homotopy::Complex<homotopy::MPFI> > > jacV;
  homotopy::MPFI h = .5;
  homotopy::Complex<homotopy::MPFI> hC(h,0);
  homotopy::MPFI gammavalue;
  homotopy::MPFR min_norm_jac;
  homotopy::MPFR min_norm_sys;
  homotopy::MPFR Amat_norm;
    homotopy::MPFR c=0.1983;
    homotopy::MPFR t=2;

    std::vector<std::vector<homotopy::Complex<homotopy::MPFI> > > Amat;
  // Initialize Aoper
  for (int i = 0; i < 2; ++i) {
    std::vector<homotopy::Complex<homotopy::MPFI> > lRow(2, 0);
    Amat.push_back(lRow);
  }
  jacV = Amat;
    std::vector<std::vector<homotopy::Complex<homotopy::MPFI> > > Vproj;
  // Initialize Vproj
  for (int i = 0; i < 2; ++i) {
    std::vector<homotopy::Complex<homotopy::MPFI> > lRow(2, 0);
    Vproj.push_back(lRow);
  }
        matrix_mult(Vproj,v,rv);

     // Open a file for writing
    std::ofstream outputFile("output_Double_rad_7.txt");

    // Check if the file is successfully opened
    if (!outputFile.is_open()) {
        std::cerr << "Failed to open the output file." << std::endl;
        return 1;  // Return an error code
    }

    int a =1;
    int i = 0;

  for (std::list<std::vector<homotopy::Complex<homotopy::MPFI> > >::iterator pt_iter = points.begin(); pt_iter != points.end(); ++pt_iter) {
    const auto beforeall = clock::now();


    for (std::list<std::list<polyTerm> >::iterator poly = system.begin(); poly != system.end(); ++poly) {
    std::vector<std::vector<homotopy::Complex<homotopy::MPFI> > > Aoper;
    std::vector<std::vector<homotopy::Complex<homotopy::MPFI> > > dummyAoper;
      Aoper=Amat;
  dummyAoper = Amat;

      std::list<std::list<polyTerm> > polyList;
      polyList.push_back(*poly);
      // std::vector<std::vector<homotopy::Complex<homotopy::Interval> > > jac;
      std::vector<std::vector<std::list<polyTerm> > > jac = computeJacobian2(polyList,2);
      std::list<std::list<polyTerm> > jacList = jacToList(jac,2);
      jac = computeJacobian(jacList);

    //    const auto beforejac = clock::now();
    //    jacinv = computeJacobian(system, *pt_iter);
    jacinv = evalJacobian(jac, *pt_iter);
    //const sec durationjac = clock::now() - beforejac;
    //  std::cout<< "jacobian took" << durationjac.count() << std::endl;
        matrix_mult(Aoper,rv,jacinv);
        matrix_mult(dummyAoper,Aoper,Vproj);
  	for (int j = 0; j < 2; ++j) {
	  Amat[i][j]= hC * dummyAoper[0][j];
	}
	i++;
  }
           jacinv = evalJacobian(jaco, *pt_iter);
	   matrix_mult(jacV,jacinv,v);
	   min_norm_jac = min_euc_int_norm(jacV);
  
  for (int i = 0; i < 2; ++i) {
    for (int j = 0; j < 2; ++j) {
      Amat[i][j]= jacinv[i][j]+Amat[i][j];
    }
  }

  jacinv = matrix_inv(Amat);
  Amat_norm = min_euc_int_norm(jacinv);

  min_norm_sys = min_euc_int_norm(evaluate(system, *pt_iter));
  

  gammavalue = gamma(system, *pt_iter, jacinv);

  outputFile << "values for " << a << " are : \n";
    outputFile << "gamma : " << gammavalue << "\n";
   outputFile << "lhs : " << min_norm_sys+ min_norm_jac* c/(t*gammavalue) << "\n";
   outputFile << "rhs : " << exp(c,4)/(t*t*Amat_norm*gammavalue*gammavalue) << "\n";
   outputFile << "radius : " << c/(t*gammavalue) << "\n";

  }

      // Print inverse A matrix
  std::cout << "Linv (A^(-1)):\n";
  for (const auto& row : jacinv) {
    for (const auto& element : row) {
      std::cout << element << " ";
    }
    std::cout << "\n";
  }

  outputFile << failCount ;

      // Close the file
    outputFile.close();


  return 0;
}
