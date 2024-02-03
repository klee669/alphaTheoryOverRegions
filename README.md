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
