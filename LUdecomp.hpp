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
#include "polyTools.hpp"



homotopy::Complex<homotopy::Interval> complex_interval_div(homotopy::Complex<homotopy::Interval> &a, homotopy::Complex<homotopy::Interval> &b) {
  homotopy::Interval areal(real(a));
  homotopy::Interval aimag(imag(a));
  homotopy::Interval breal(real(b));
  homotopy::Interval bimag(imag(b));

  homotopy::Interval realPart = (areal*breal+aimag*bimag)/(breal*breal+bimag*bimag);
  homotopy::Interval imagPart = (aimag*breal-areal*bimag)/(breal*breal+bimag*bimag);

  homotopy::Complex<homotopy::Interval> value(realPart,imagPart);

  return value;
}

homotopy::Complex<homotopy::Interval> complex_interval_div_LU(homotopy::Complex<homotopy::Interval> &a, homotopy::Complex<homotopy::Interval> b) {
  homotopy::Complex<homotopy::Interval> value;

  if (a == b) {
    value = oneI;
  } else {    
    homotopy::Interval areal(real(a));
    homotopy::Interval aimag(imag(a));
    homotopy::Interval breal(real(b));
    homotopy::Interval bimag(imag(b));
    
    homotopy::Interval realPart = (areal*breal+aimag*bimag)/(breal*breal+bimag*bimag);
    homotopy::Interval imagPart = (aimag*breal-areal*bimag)/(breal*breal+bimag*bimag);
    homotopy::Complex<homotopy::Interval> intv(realPart,imagPart); 

    value = intv;
  }

  return value;
}

homotopy::Complex<homotopy::Interval> complex_interval_sub_LU(homotopy::Complex<homotopy::Interval> &a, homotopy::Complex<homotopy::Interval> b) {
  homotopy::Complex<homotopy::Interval> value;

  if (a == b) {
    value = zeroI;
  } else {    
    value = a - b;
  }
  return value;
}





int detectMax(std::vector<std::vector<homotopy::Complex<homotopy::Interval> > >& A, int i) {
  double maxValue = abs(A[i][i]);  // Initialize 
  int n = A.size();

  //  int row = i;
  int position = i;
  homotopy::Complex<homotopy::Interval> complexNumber;
  double normValue;
  
  for (int j = i; j < n; ++j) {
    complexNumber = A[j][i];
    normValue = abs(complexNumber);
    if (normValue > maxValue) {
      maxValue = normValue;
      position = j;
    }
    //   ++row;
  }


  return position;
}


void rowSwap(std::vector<std::vector<homotopy::Complex<homotopy::Interval> > >& A, int row1, int row2) {

  // Get iterators to the desired positions
  auto it1 = std::next(A.begin(), row1);
  auto it2 = std::next(A.begin(), row2);
  
  // Swap the elements at the specified positions
  std::swap(*it1, *it2);

  //  return A;
}


void matrix_mult(std::vector<std::vector<homotopy::Complex<homotopy::Interval> > >& result, std::vector<std::vector<homotopy::Complex<homotopy::Interval> > >& A,std::vector<std::vector<homotopy::Complex<homotopy::Interval> > >& B) {
  //  std::vector<std::vector<homotopy::Complex<homotopy::Interval> > > result;
  //std::vector<homotopy::Complex<homotopy::Interval> > resultRow;
  //  homotopy::Complex<homotopy::Interval> sum;

  int numrowsA = A.size();
  int numcolsA = A[0].size();

  int numrowsB = B.size();
  int numcolsB = B[0].size();

  
  for (int i = 0; i < numrowsA; ++i) {
      for (int j = 0; j < numcolsB; ++j) {
	//	sum = zeroI;
	result[i][j] = zeroI;
	for (int k = 0; k < numcolsA; ++k) {
	  result[i][j] += A[i][k]*B[k][j];
	}
	//	std::cout << result[i][j] << std::endl;
	//	result[i][j] = sum;
	//resultRow.push_back(sum);
      }
      //      result.push_back(resultRow);
      //    resultRow.clear();
  }
  //  return result;
}


bool pivotCheck(std::vector<std::vector<homotopy::Complex<homotopy::Interval> > >& A, int i) {
  double abs_pivot;
  abs_pivot = abs(A[i][i]);
  return abs_pivot < .0001;
}


void luDecomposition(std::vector<std::vector<homotopy::Complex<homotopy::Interval> > >& A, std::vector<std::vector<homotopy::Complex<homotopy::Interval> > >& L, std::vector<std::vector<homotopy::Complex<homotopy::Interval> > >& U, std::vector<std::vector<homotopy::Complex<homotopy::Interval> > >& P) {
  int n = A.size();

  homotopy::Complex<homotopy::Interval> sum = zeroI;
  std::vector<std::vector<homotopy::Complex<homotopy::Interval> > > dummyP, iD, Linv, Ldummy; 
  std::vector<std::vector<std::vector<homotopy::Complex<homotopy::Interval> > > > dummyPlist, LdummyList; 
  
  
  // Initialize L and U
  for (int i = 0; i < n; ++i) {
    std::vector<homotopy::Complex<homotopy::Interval> > lRow(n, 0);
    L.push_back(lRow);
  }
  dummyP = L;
  
  for (int i = 0; i < n; ++i) {
    dummyP[i][i] = oneI;
  }
  iD = dummyP;
  L = dummyP;
  P = dummyP;
  Linv = dummyP;
  U = A;
  Ldummy = dummyP;

  


  
  for (int j = 0; j < n-1; j++) {
    int maxEntry;
    maxEntry = detectMax(U,j);
    // dummyP =
      rowSwap(dummyP,j,maxEntry);
    dummyPlist.push_back(dummyP);

    matrix_mult(Ldummy,dummyP,U);

    U = Ldummy;
    Ldummy = iD;


    for (int i = j; i < n-1; i++) {
      Ldummy[i+1][j] = complex_interval_sub_LU(zeroI,complex_interval_div_LU(U[i+1][j],U[j][j]));
      Linv[i+1][j] = complex_interval_div_LU(U[i+1][j],U[j][j]);
    }
    LdummyList.push_back(Linv);

        

    
    matrix_mult(L,dummyP, Linv);
    Linv = L;
    L = iD;

    matrix_mult(dummyP, Ldummy, U);
    U = dummyP;
    dummyP = iD;
    Ldummy = iD;
    Linv = iD;
  }
  for (int j = 0; j < n-1; j++) {
    Ldummy = LdummyList[j];
    for (int i = j+1; i < n-1; i++) {
      matrix_mult(Linv,dummyPlist[i],Ldummy);
      matrix_mult(Ldummy, Linv, dummyPlist[i]);
    }
    matrix_mult(Linv,L, Ldummy);
    L = Linv;
    Linv = iD;
    matrix_mult(dummyP,dummyPlist[j],P);
    P = dummyP;
    dummyP = iD;
  }

}


std::vector<std::vector<homotopy::Complex<homotopy::Interval> > > matrix_inv(std::vector<std::vector<homotopy::Complex<homotopy::Interval> > >& A) {
  int n = A.size();
  std::vector<std::vector<homotopy::Complex<homotopy::Interval> > > result;
  std::vector<std::vector<homotopy::Complex<homotopy::Interval> > > Adummy;

  std::vector<std::vector<homotopy::Complex<homotopy::Interval> > > L;
  std::vector<std::vector<homotopy::Complex<homotopy::Interval> > > U;
  std::vector<std::vector<homotopy::Complex<homotopy::Interval> > > P;
  //    std::vector<std::vector<homotopy::Complex<homotopy::Interval> > > LU;
  //  std::vector<std::vector<homotopy::Complex<homotopy::Interval> > > LUP;

  int maxEntry;

  Adummy = A;
  
  maxEntry = 0;


  // Perform LU decomposition
  luDecomposition(A, L, U, P);
  result = A;


    int numPerms = P.size();
  std::vector<std::vector<homotopy::Complex<homotopy::Interval> > > Linv;
  std::vector<std::vector<homotopy::Complex<homotopy::Interval> > > Uinv;

  // Initialize Linv and Uinv
  for (int i = 0; i < n; ++i) {
    std::vector<homotopy::Complex<homotopy::Interval> > lRow(n, 0);
    Linv.push_back(lRow);
    Uinv.push_back(lRow);
  }


  
  for (int i = 0; i < n; ++i) {
    Uinv[i][i] = oneI/U[i][i];
  }

  homotopy::Complex<homotopy::Interval> sum;
  
  for (int i = 1; i < n; ++i) {
    for (int k = i; k < n; ++k) {
      sum = zeroI;
      for (int j = 0; j < k; ++j) {
	sum -= U[j][k]*Uinv[k-i][j];
      }
      Uinv[k-i][k] = complex_interval_div_LU(sum,U[k][k]);
      //          std::cout << Uinv[k-i][k] << std::endl;
    }
  }

  for (int i = 0; i < n; ++i) {
    Linv[i][i] = 1;
  }

  
  for (int i = 1; i < n; ++i) {
    for (int k = i; k < n; ++k) {
      sum = zeroI;
      for (int j = 0; j < k; ++j) {
	sum -= L[j+1][k-i]*Linv[k][j+1];
      }
      Linv[k][k-i] = complex_interval_div_LU(sum,L[k-i][k-i]);
    }
  }



  


  
  

  
  matrix_mult(result,Uinv, Linv);
matrix_mult(A,result, P);
 result = A;
  A = Adummy;
  


  

  return result;
}
