#include "s21_matrix_oop.h"

S21Matrix::S21Matrix() {
  rows_ = 3;
  cols_ = 3;
  matrix_ = new double*[rows_];
  for (int i = 0; i < rows_; i++) {
    matrix_[i] = new double[cols_];
  }
}

S21Matrix::S21Matrix(int num) : rows_(num), cols_(num) {
  if (rows_ <= 0 || cols_ <= 0) throw std::out_of_range(INVALID);
  matrix_ = new double*[rows_];
  for (int i = 0; i < rows_; i++) {
    matrix_[i] = new double[cols_];
  }
}

S21Matrix::S21Matrix(int rows, int cols) : rows_(rows), cols_(cols) {
  if (rows_ <= 0 || cols_ <= 0) throw std::out_of_range(INVALID);
  matrix_ = new double*[rows_];
  for (int i = 0; i < rows_; i++) {
    matrix_[i] = new double[cols_];
  }
}

S21Matrix::S21Matrix(const S21Matrix& other)
    : rows_(other.rows_), cols_(other.cols_) {
  if (rows_ <= 0 || cols_ <= 0) throw std::out_of_range(INVALID);
  matrix_ = new double*[other.rows_];
  for (int i = 0; i < other.rows_; i++) {
    matrix_[i] = new double[other.cols_];
    for (int j = 0; j < cols_; j++) {
      matrix_[i][j] = other.matrix_[i][j];
    }
  }
}

S21Matrix::S21Matrix(S21Matrix&& other)
    : rows_(other.rows_), cols_(other.cols_) {
  if (rows_ <= 0 || cols_ <= 0) throw std::out_of_range(INVALID);
  std::swap(matrix_, other.matrix_);
  other.matrix_ = NULL;
}

S21Matrix::~S21Matrix() {
  if (matrix_) {
    for (int i = 0; i < rows_; i++) {
      if (matrix_[i]) {
        delete[] matrix_[i];
      }
    }
    delete[] matrix_;
  }
}

bool S21Matrix::EqMatrix(const S21Matrix& other) {
  bool result = true;
  if (rows_ != other.rows_ || cols_ != other.cols_) result = false;
  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < cols_; j++) {
      if (fabs(matrix_[i][j] - other.matrix_[i][j]) >= 1e-7) result = false;
    }
  }
  return result;
}

void S21Matrix::SumMatrix(const S21Matrix& other) {
  if (rows_ != other.rows_ || cols_ != other.cols_)
    throw std::out_of_range(NOT_EQ);
  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < cols_; j++) matrix_[i][j] += other.matrix_[i][j];
  }
}

void S21Matrix::SubMatrix(const S21Matrix& other) {
  if (rows_ != other.rows_ || cols_ != other.cols_)
    throw std::out_of_range(NOT_EQ);
  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < cols_; j++) matrix_[i][j] -= other.matrix_[i][j];
  }
}

void S21Matrix::MulNumber(const double num) {
  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < cols_; j++) matrix_[i][j] *= num;
  }
}

void S21Matrix::MulMatrix(const S21Matrix& other) {
  if (cols_ != other.rows_) throw std::out_of_range(NO_MULT);
  S21Matrix result(rows_, other.cols_);
  for (int i = 0; i < result.rows_; i++) {
    for (int j = 0; j < result.cols_; j++) {
      result.matrix_[i][j] = 0;
      for (int k = 0; k < cols_; k++)
        result.matrix_[i][j] += matrix_[i][k] * other.matrix_[k][j];
    }
  }
  std::swap(matrix_, result.matrix_);
  int temp = cols_;
  cols_ = result.cols_;
  result.cols_ = temp;
}

S21Matrix S21Matrix::Transpose() {
  S21Matrix result(rows_, cols_);
  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < cols_; j++) result.matrix_[j][i] = matrix_[i][j];
  }
  return result;
}

S21Matrix S21Matrix::CalcComplements() {
  if (rows_ != cols_) throw std::out_of_range(NOT_SQUARE);
  S21Matrix result(rows_, cols_);
  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < cols_; j++) {
      S21Matrix minor_matrix(minorMatrix(i, j));
      double det = minor_matrix.Determinant();
      result.matrix_[i][j] = det * ((i + j) % 2 == 0 ? 1 : -1);
    }
  }
  return result;
}  // LCOV_EXCL_LINE

double S21Matrix::Determinant() {
  if (rows_ != cols_) throw std::out_of_range(NOT_SQUARE);
  double result = 0;
  if (rows_ == 1)
    result = matrix_[0][0];
  else if (rows_ == 2)
    result = matrix_[0][0] * matrix_[1][1] - matrix_[0][1] * matrix_[1][0];
  else {
    for (int i = 0; i < rows_; i++) {
      S21Matrix minor_matrix(minorMatrix(0, i));
      result +=
          matrix_[0][i] * minor_matrix.Determinant() * (i % 2 == 0 ? 1 : -1);
    }
  }
  return result;
}

S21Matrix S21Matrix::InverseMatrix() {
  if (rows_ != cols_) throw std::out_of_range(NOT_SQUARE);
  double det = this->Determinant();
  if (fabs(det) < 1e-7) throw std::out_of_range(NULL_DET);
  S21Matrix comps(this->CalcComplements());
  comps.MulNumber(1 / det);
  S21Matrix result(comps.Transpose());
  return result;
}

S21Matrix S21Matrix::minorMatrix(int row, int col) {
  S21Matrix result(rows_ - 1, cols_ - 1);
  for (int i = 0, minor_i = 0; i < rows_; i++) {
    if (i == row) continue;
    for (int j = 0, minor_j = 0; j < cols_; j++) {
      if (j == col) continue;
      result.matrix_[minor_i][minor_j] = matrix_[i][j];
      minor_j++;
    }
    minor_i++;
  }
  return result;
}

void S21Matrix::printMatrix() {
  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < cols_; j++) {
      std::cout << matrix_[i][j] << " ";
    }
    std::cout << "\n";
  }
}

void S21Matrix::initMatrix() {
  int count = 1;
  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < cols_; j++) {
      matrix_[i][j] = count;
      count++;
    }
  }
}

void S21Matrix::initMatrix(double num) {
  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < cols_; j++) {
      matrix_[i][j] = num;
    }
  }
}

int S21Matrix::getRows() { return rows_; }

int S21Matrix::getCols() { return cols_; }

S21Matrix S21Matrix::setRows(int rows) {
  S21Matrix temp(rows, cols_);
  for (int i = 0; i < temp.rows_; i++) {
    for (int j = 0; j < temp.cols_; j++) {
      if (i < rows_)
        temp.matrix_[i][j] = matrix_[i][j];
      else
        temp.matrix_[i][j] = 0;
    }
  }
  std::swap(matrix_, temp.matrix_);
  int tmp = rows_;
  rows_ = rows;
  temp.rows_ = tmp;
  return *this;
}

S21Matrix S21Matrix::setCols(int cols) {
  S21Matrix temp(rows_, cols);
  for (int i = 0; i < temp.rows_; i++) {
    for (int j = 0; j < temp.cols_; j++) {
      if (j < cols_)
        temp.matrix_[i][j] = matrix_[i][j];
      else
        temp.matrix_[i][j] = 0;
    }
  }
  std::swap(matrix_, temp.matrix_);
  int tmp = cols_;
  cols_ = cols;
  temp.cols_ = tmp;
  return *this;
}

S21Matrix S21Matrix::operator+(const S21Matrix& other) {
  this->SumMatrix(other);
  return *this;
}

S21Matrix S21Matrix::operator-(const S21Matrix& other) {
  this->SubMatrix(other);
  return *this;
}

S21Matrix S21Matrix::operator*(const double num) {
  this->MulNumber(num);
  return *this;
}

S21Matrix S21Matrix::operator*(const S21Matrix& other) {
  this->MulMatrix(other);
  return *this;
}

S21Matrix& S21Matrix::operator+=(const S21Matrix& other) {
  this->SumMatrix(other);
  return *this;
}

S21Matrix& S21Matrix::operator-=(const S21Matrix& other) {
  this->SubMatrix(other);
  return *this;
}

S21Matrix& S21Matrix::operator*=(const double num) {
  this->MulNumber(num);
  return *this;
}

S21Matrix& S21Matrix::operator*=(const S21Matrix& other) {
  this->MulMatrix(other);
  return *this;
}

S21Matrix& S21Matrix::operator=(const S21Matrix& other) {
  S21Matrix temp(other);
  std::swap(matrix_, temp.matrix_);
  rows_ = temp.rows_;
  cols_ = temp.cols_;
  return *this;
}

bool S21Matrix::operator==(const S21Matrix& other) {
  return this->EqMatrix(other);
}

double& S21Matrix::operator()(int row, int col) {
  if (row >= rows_ || col >= cols_) throw std::out_of_range(NO_INDEX);
  return matrix_[row][col];
}