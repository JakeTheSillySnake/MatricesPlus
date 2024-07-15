#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "s21_matrix_oop.h"

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  ::testing::InitGoogleMock(&argc, argv);

  return RUN_ALL_TESTS();
}

TEST(test_Constr, constr_1) {
  S21Matrix A(2);
  A.initMatrix();
  S21Matrix B = std::move(A);
  ASSERT_TRUE(B.getRows() == 2);
}

TEST(test_Constr, constr_2) {
  S21Matrix *A = new S21Matrix();
  int result = A->getCols();
  delete A;
  ASSERT_TRUE(result == 3);
}

TEST(test_EqMatrix, eq_matrix_1) {
  S21Matrix A(3);
  S21Matrix B(3);
  A.initMatrix();
  B.initMatrix();
  ASSERT_TRUE(A.EqMatrix(B));
}

TEST(test_EqMatrix, eq_matrix_2) {
  S21Matrix A(3);
  S21Matrix B(4);
  A.initMatrix();
  B.initMatrix();
  ASSERT_TRUE(!A.EqMatrix(B));
}

TEST(test_SumMatrix, sum_matrix_1) {
  S21Matrix A(3);
  S21Matrix B(3);
  S21Matrix C(3);
  A.initMatrix(1);
  B.initMatrix(2);
  C.initMatrix(3);
  A.SumMatrix(B);
  ASSERT_TRUE(A.EqMatrix(C));
}

TEST(test_SumMatrix, sum_matrix_2) {
  S21Matrix A(3);
  S21Matrix B(4);
  A.initMatrix();
  B.initMatrix();
  EXPECT_ANY_THROW(A.SumMatrix(B));
}

TEST(test_SubMatrix, sub_matrix_1) {
  S21Matrix A(3);
  S21Matrix B(3);
  S21Matrix C(3);
  A.initMatrix(1);
  B.initMatrix(2);
  C.initMatrix(-1);
  A.SubMatrix(B);
  ASSERT_TRUE(A.EqMatrix(C));
}

TEST(test_SubMatrix, sub_matrix_2) {
  S21Matrix A(3);
  S21Matrix B(4);
  A.initMatrix();
  B.initMatrix();
  EXPECT_ANY_THROW(A.SubMatrix(B));
}

TEST(test_MulNumber, mul_num_1) {
  S21Matrix A(3);
  S21Matrix B(3);
  A.initMatrix(2);
  B.initMatrix(6);
  A.MulNumber(3);
  ASSERT_TRUE(A.EqMatrix(B));
}

TEST(test_MulMatrix, mul_matrix_1) {
  S21Matrix A(1, 2);
  S21Matrix B(2, 1);
  A.initMatrix(2);
  B.initMatrix(3);
  A.MulMatrix(B);
  ASSERT_TRUE(A(0, 0) == 12.0);
}

TEST(test_MulMatrix, mul_matrix_2) {
  S21Matrix A(1);
  S21Matrix B(2);
  A.initMatrix(2);
  B.initMatrix(3);
  EXPECT_ANY_THROW(A.MulMatrix(B));
}

TEST(test_Transpose, transpose_1) {
  S21Matrix A(2);
  A(0, 0) = 1;
  A(0, 1) = 1;
  A(1, 0) = 2;
  A(1, 1) = 2;
  S21Matrix B = A.Transpose();
  ASSERT_TRUE(B(0, 1) == 2.0);
}

TEST(test_CalcComplements, calc_complements_1) {
  S21Matrix A(2, 1);
  A.initMatrix();
  EXPECT_ANY_THROW(A.CalcComplements());
}

TEST(test_Determinant, determinant_1) {
  S21Matrix A(2, 1);
  A.initMatrix();
  EXPECT_ANY_THROW(A.Determinant());
}

TEST(test_Inverse, inverse_1) {
  S21Matrix A(2);
  A.initMatrix();
  S21Matrix B = A.InverseMatrix();
  A(0, 0) = -2;
  A(0, 1) = 1;
  A(1, 0) = 1.5;
  A(1, 1) = -0.5;
  ASSERT_TRUE(B.EqMatrix(A));
}

TEST(test_Inverse, inverse_2) {
  S21Matrix A(2, 1);
  A.initMatrix();
  EXPECT_ANY_THROW(A.InverseMatrix());
}

TEST(test_Inverse, inverse_3) {
  S21Matrix A(3);
  A.initMatrix();
  EXPECT_ANY_THROW(A.InverseMatrix());
}

TEST(test_GetRows, get_rows_1) {
  S21Matrix A(2, 1);
  ASSERT_TRUE(A.getRows() == 2);
}

TEST(test_GetCols, get_cols_1) {
  S21Matrix A(2, 1);
  ASSERT_TRUE(A.getCols() == 1);
}

TEST(test_SetRows, set_rows_1) {
  S21Matrix A(2, 1);
  A.initMatrix();
  A.setRows(4);
  ASSERT_TRUE(A(3, 0) == 0);
}

TEST(test_SetRows, set_rows_2) {
  S21Matrix A(2, 1);
  A.initMatrix();
  EXPECT_ANY_THROW(A.setRows(0));
}

TEST(test_SetCols, set_cols_1) {
  S21Matrix A(2, 1);
  A.initMatrix();
  A.setCols(4);
  ASSERT_TRUE(A(0, 3) == 0);
}

TEST(test_SetCols, set_cols_2) {
  S21Matrix A(2, 1);
  A.initMatrix();
  EXPECT_ANY_THROW(A.setCols(0));
}

TEST(test_Operator, operator_1) {
  S21Matrix A(2);
  S21Matrix B(2);
  A.initMatrix(1);
  B.initMatrix(2);
  A + B;
  A.printMatrix();
  ASSERT_TRUE(A(0, 0) == 3.0);
}

TEST(test_Operator, operator_2) {
  S21Matrix A(2);
  S21Matrix B(2);
  A.initMatrix(1);
  B.initMatrix(2);
  A - B;
  ASSERT_TRUE(A(0, 0) == -1.0);
}

TEST(test_Operator, operator_3) {
  S21Matrix A(1, 2);
  S21Matrix B(2, 1);
  A.initMatrix(2);
  B.initMatrix(3);
  A *B;
  ASSERT_TRUE(A(0, 0) == 12.0);
}

TEST(test_Operator, operator_4) {
  S21Matrix A(2);
  A.initMatrix(2);
  A * 2.0;
  ASSERT_TRUE(A(0, 0) == 4.0);
}

TEST(test_Operator, operator_5) {
  S21Matrix A(2);
  S21Matrix B(2);
  A.initMatrix(1);
  B.initMatrix(2);
  A += B;
  ASSERT_TRUE(A(0, 0) == 3.0);
}

TEST(test_Operator, operator_6) {
  S21Matrix A(2);
  S21Matrix B(2);
  A.initMatrix(1);
  B.initMatrix(2);
  A -= B;
  A.printMatrix();
  ASSERT_TRUE(A(0, 0) == -1.0);
}

TEST(test_Operator, operator_7) {
  S21Matrix A(1, 2);
  S21Matrix B(2, 1);
  A.initMatrix(2);
  B.initMatrix(3);
  A *= B;
  ASSERT_TRUE(A(0, 0) == 12.0);
}

TEST(test_Operator, operator_8) {
  S21Matrix A(2);
  A.initMatrix(2);
  A *= 2.0;
  ASSERT_TRUE(A(0, 0) == 4.0);
}

TEST(test_Operator, operator_9) {
  S21Matrix A(2);
  S21Matrix B(2);
  A.initMatrix(2);
  B.initMatrix(2);
  ASSERT_TRUE(A == B);
}

TEST(test_Operator, operator_10) {
  S21Matrix A(2);
  S21Matrix B(2);
  B.initMatrix(2);
  A = B;
  ASSERT_TRUE(A(0, 0) == 2.0);
}