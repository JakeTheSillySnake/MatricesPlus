#ifndef MATRIX_H
#define MATRIX_H

#include <cmath>
#include <cstring>
#include <iostream>

// error messages
#define INVALID "Incorrect input: Matrix is invalid"
#define NOT_SQUARE "Incorrect input: Matrix should be square"
#define NOT_EQ "Incorrect input: Matrices should be the same size"
#define NO_MULT "Incorrect input: A's columns not equal to B's rows"
#define NULL_DET "Matrix determinant is zero"
#define NO_INDEX "Incorrect input: Index out of range"

class S21Matrix {
 private:
  // attributes
  int rows_, cols_;
  double** matrix_;

 public:
  S21Matrix();                        // default constructor
  explicit S21Matrix(int num);        // square matrix constructor
  S21Matrix(int rows, int cols);      // constructor with params
  S21Matrix(const S21Matrix& other);  // copy constructor
  S21Matrix(S21Matrix&& other);       // move constructor
  ~S21Matrix();                       // destructor

  // public methods
  bool EqMatrix(const S21Matrix& other);
  void SumMatrix(const S21Matrix& other);
  void SubMatrix(const S21Matrix& other);
  void MulNumber(const double num);
  void MulMatrix(const S21Matrix& other);
  S21Matrix Transpose();
  S21Matrix CalcComplements();
  double Determinant();
  S21Matrix InverseMatrix();
  S21Matrix minorMatrix(int row, int col);
  void printMatrix();
  void initMatrix();
  void initMatrix(double num);

  // accessors & mutators
  int getRows();
  int getCols();
  S21Matrix setRows(int rows);
  S21Matrix setCols(int cols);

  // operator overloads
  S21Matrix operator+(const S21Matrix& other);
  S21Matrix operator-(const S21Matrix& other);
  S21Matrix operator*(const double num);
  S21Matrix operator*(const S21Matrix& other);
  S21Matrix& operator+=(const S21Matrix& other);
  S21Matrix& operator-=(const S21Matrix& other);
  S21Matrix& operator*=(const double num);
  S21Matrix& operator*=(const S21Matrix& other);
  S21Matrix& operator=(const S21Matrix& other);
  bool operator==(const S21Matrix& other);
  double& operator()(int row, int col);
};

#endif