# alphaTheoryOverRegions

### How to Run Codes ###

* go to the directory of the repository

* open 'constants.cpp' or 'constantsMPFR.cpp' file depending on the precision required to use.
* For 'constantsMPFR.cpp' file, go to the line of 'mpfr_set_default_prec(256);' and change 256 to a desired precision.
* go to the line of 'realP(realPart-pow(10,-20),realPart+pow(10,-20));' and 'imagP(imagPart-pow(10,-20),imagPart+pow(10,-20));' and change pow(10,-20) to a desired radius.
* go up the line of 'std::ofstream outputFile("output.txt");' and change what is inside the quote mark to the desired name of the output.

* run code in the terminal

```
clang++ -std=c++11 -o output constants.cpp -I /path/to/mpfr/mpfi/gmp/include/directory -L /path/to/mpfr/mpfi/gmp/library/directory -lmpfi -lmpfr -lgmp -Os
```
This makes an executable file named 'output'.

* run code with the input file 'input.txt'

```
./output input.txt
```


### Input file format ###

The code is run by the file listing polynomials and points.

** polynomials

* The first block consists of polynomials described in their coefficients and the degrees of each monomial.
* The file contains a block for each polynomial which contains the real and imaginary parts of its coefficient followed by the degrees of each variable in the monomial. For example, if the polynomial system depends upon three variables, say $x, y$ and $z$, the term $(1/2+3i)*x*y^2*z$ is represented by the row $1/2 3 1 2 1$
* Likewise, the polynomial $x^2-3xy+z^4$ is represented by the row $\begin{array}{ccccc} 1 & 0 & 2 & 0& 0 \\ -3 & 0 & 1 & 0 & 0 \\ 1 & 0 & 0 & 0 & 4\end{array}$
