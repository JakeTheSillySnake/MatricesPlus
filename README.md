# MatricesPlus

Implementation of basic operations with matrices: create_matrix (creation), remove_matrix (cleaning and destruction), eq_matrix (comparison), sum_matrix (addition), sub_matrix (subtraction), mult_matrix (multiplication), mult_number (multiplication by number), transpose (transpose), determinant (calculation of determinant), calc_complements (calculation of matrix of algebraic complements), inverse_matrix (finding inverse of the matrix).

## The project adheres to the following guidelines:

- The library must be developed in C language of C11 standard using gcc compiler;
- The library code must be located in the src folder on the develop branch;
- Do not use outdated and legacy language constructions and library functions. Pay attention to the legacy and obsolete marks in the official documentation on the language and the libraries used. Use the POSIX.1-2017 standard;
When writing code it is necessary to follow the Google style;
- Make it as a static library named s21_matrix.a (with the s21_matrix.h header file);
- The library must be developed according to the principles of structured programming;
- Use prefix s21_ before each function;
- Prepare full coverage of library functions code with unit-tests using the Check library;
- Unit tests must cover at least 80% of each function (checked using gcov);
- Provide a Makefile for building the library and tests (with targets all, clean, test, s21_matrix.a, gcov_report);
- The gcov_report target should generate a gcov report in the form of an html page. Unit tests must be run with gcov flags to do this;
- The matrix must be implemented as the structure described above;
- Verifiable accuracy of the fractional part is up to 6 decimal places.
