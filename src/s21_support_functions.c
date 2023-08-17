#include "s21_matrix.h"

int check_incorrect_matrix(matrix_t *target) {
  return (target && target->columns > 0 && target->rows > 0);
}

int check_valid_eq_matrix(matrix_t *A, matrix_t *B) {
  int status = OK;
  if (!(check_incorrect_matrix(A) && check_incorrect_matrix(B))) {
    status = INCORRECT_MATRIX;
  } else if ((A->rows - B->rows) || (A->columns - B->columns)) {
    status = ERROR_CALCULATING;
  }
  return status;
}

int check_valid_arithmetic_matrix(matrix_t *A, matrix_t *B, matrix_t *result) {
  int status = OK;
  if (!(check_incorrect_matrix(A) && check_incorrect_matrix(B) &&
        result != NULL))
    status = INCORRECT_MATRIX;
  else if (!(A->rows == B->rows && A->columns == B->columns))
    status = ERROR_CALCULATING;
  return status;
}

int check_valid_multiplication_matrix(matrix_t *A, matrix_t *B,
                                      matrix_t *result) {
  int status = OK;
  if (!(check_incorrect_matrix(A) && check_incorrect_matrix(B) &&
        result != NULL))
    status = INCORRECT_MATRIX;
  else if (A->columns != B->rows)
    status = ERROR_CALCULATING;
  return status;
}

double find_determinant(matrix_t temp) {
  int count = 0;
  double first = 0.0;
  double determinant = 1.0;
  for (; count < temp.rows; count++) {
    int i = count;
    if (fabs(temp.matrix[count][count]) < 1e-7) {
      for (; i < temp.rows && fabs(temp.matrix[i][count]) < 1e-7; i++) {
      }
      if (i == temp.rows) {
        return 0.0;
      } else {
        determinant *= -1;
        swap_rows(&temp, i, count);
      }
    }
    determinant *= temp.matrix[count][count];
    for (int row = count; row < temp.rows; row++) {
      first = temp.matrix[row][count];
      for (int col = count; col < temp.columns; col++) {
        if (row == count) {
          temp.matrix[row][col] /= first;
        } else {
          temp.matrix[row][col] -= temp.matrix[count][col] * first;
        }
      }
    }
  }
  return determinant;
}

void swap_rows(matrix_t *temp, int i, int k) {
  double d = 0.0;
  for (int j = k; j < temp->rows; j++) {
    d = temp->matrix[i][j];
    temp->matrix[i][j] = temp->matrix[k][j];
    temp->matrix[k][j] = d;
  }
}

void create_minor(int erase_row, int erase_col, matrix_t *temp, matrix_t *A) {
  int x = 0, y = 0;
  for (int i = 0; i < A->rows; i++) {
    for (int j = 0; j < A->rows; j++) {
      if (i != erase_row && j != erase_col) {
        temp->matrix[x][y] = A->matrix[i][j];
        y++;
      }
    }
    y = 0;
    if (i != erase_row) x++;
  }
}
