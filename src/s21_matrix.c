#include "s21_matrix.h"

int s21_create_matrix(int rows, int columns, matrix_t *result) {
  int status = OK;
  if (!(columns > 0 && rows > 0 && result != NULL)) {
    status = INCORRECT_MATRIX;
  }
  if (!status) {
    result->matrix = NULL;
    result->rows = rows;
    result->columns = columns;
  }
  if (!status) {
    if (!(result->matrix = calloc(rows, sizeof(double *)))) {
      status = INCORRECT_MATRIX;
      result->matrix = NULL;
    }
  }
  for (int i = 0; i < rows && !status; i++)
    if (!(result->matrix[i] = calloc(columns, sizeof(double)))) {
      result->rows = i;
      status = 1;
      s21_remove_matrix(result);
    }
  return status;
}

void s21_remove_matrix(matrix_t *A) {
  if (A) {
    for (int i = 0; i < A->rows; i++) free(A->matrix[i]);
    free(A->matrix);
    A->matrix = NULL;
    A->rows = 0;
    A->columns = 0;
  }
}

int s21_eq_matrix(matrix_t *A, matrix_t *B) {
  int status = SUCCESS;
  if ((check_valid_eq_matrix(A, B)) != 0) {
    status = FAILURE;
  } else {
    for (int i = 0; i < A->rows; i++)
      for (int j = 0; j < A->columns; j++)
        if (fabs(A->matrix[i][j] - B->matrix[i][j]) > EPS) status = FAILURE;
  }
  return status;
}

int s21_sum_matrix(matrix_t *A, matrix_t *B, matrix_t *result) {
  int status = OK;
  if (!(status = check_valid_arithmetic_matrix(A, B, result))) {
    s21_create_matrix(A->rows, A->columns, result);
    for (int i = 0; i < A->rows; i++)
      for (int j = 0; j < A->columns; j++)
        result->matrix[i][j] = A->matrix[i][j] + B->matrix[i][j];
  }
  return status;
}

int s21_sub_matrix(matrix_t *A, matrix_t *B, matrix_t *result) {
  int status = OK;
  if (!(status = check_valid_arithmetic_matrix(A, B, result))) {
    s21_create_matrix(A->rows, A->columns, result);
    for (int i = 0; i < A->rows; i++)
      for (int j = 0; j < A->columns; j++)
        result->matrix[i][j] = A->matrix[i][j] - B->matrix[i][j];
  }
  return status;
}

int s21_mult_number(matrix_t *A, double number, matrix_t *result) {
  int status = OK;
  if (!check_incorrect_matrix(A) || result == NULL) {
    status = INCORRECT_MATRIX;
  } else {
    s21_create_matrix(A->rows, A->columns, result);
    for (int i = 0; i < result->rows; i++)
      for (int j = 0; j < result->columns; j++)
        result->matrix[i][j] = A->matrix[i][j] * number;
  }
  return status;
}

int s21_mult_matrix(matrix_t *A, matrix_t *B, matrix_t *result) {
  int status = OK;
  if (!(status = check_valid_multiplication_matrix(A, B, result)) != 0) {
    s21_create_matrix(A->rows, B->columns, result);
    for (int i = 0; i < A->rows; i++)
      for (int j = 0; j < B->columns; j++)
        for (int k = 0; k < B->rows; k++)
          result->matrix[i][j] += A->matrix[i][k] * B->matrix[k][j];
  }
  return status;
}

int s21_transpose(matrix_t *A, matrix_t *result) {
  int status = OK;
  if (!check_incorrect_matrix(A) || result == NULL) {
    status = INCORRECT_MATRIX;
  } else {
    s21_create_matrix(A->columns, A->rows, result);
    for (int i = 0; i < A->rows; i++)
      for (int j = 0; j < A->columns; j++)
        result->matrix[j][i] = A->matrix[i][j];
  }
  return status;
}

int s21_calc_complements(matrix_t *A, matrix_t *result) {
  int status = OK;
  if (!check_incorrect_matrix(A) || result == NULL) {
    status = INCORRECT_MATRIX;
  } else if (A->rows != A->columns) {
    status = ERROR_CALCULATING;
  } else {
    matrix_t temp = {0};
    s21_create_matrix(A->rows - 1, A->columns - 1, &temp);
    s21_create_matrix(A->rows, A->columns, result);
    for (int i = 0; i < A->rows; i++)
      for (int j = 0; j < A->rows; j++) {
        create_minor(i, j, &temp, A);
        result->matrix[i][j] = find_determinant(temp) * pow(-1, i + j + 2);
      }
    s21_remove_matrix(&temp);
  }
  return status;
}

int s21_determinant(matrix_t *A, double *result) {
  int status = OK;
  if (!check_incorrect_matrix(A) || result == NULL) {
    status = INCORRECT_MATRIX;
  } else if (A->rows != A->columns) {
    status = ERROR_CALCULATING;
  } else {
    *result = 0.0;
    matrix_t temp = {0};
    s21_create_matrix(A->rows, A->columns, &temp);
    for (int i = 0; i < A->rows; i++)
      for (int j = 0; j < A->rows; j++) temp.matrix[i][j] = A->matrix[i][j];
    *result = find_determinant(temp);
    s21_remove_matrix(&temp);
  }
  return status;
}

int s21_inverse_matrix(matrix_t *A, matrix_t *result) {
  int status = OK;
  if (!check_incorrect_matrix(A) || result == NULL) {
    status = INCORRECT_MATRIX;
  } else if (A->rows != A->columns) {
    status = ERROR_CALCULATING;
  } else {
    matrix_t temp1 = {0}, temp2 = {0};
    double determinant = 0.0;
    s21_determinant(A, &determinant);
    if (determinant == 0.0) return 2;
    s21_calc_complements(A, &temp1);
    s21_transpose(&temp1, &temp2);
    s21_mult_number(&temp2, 1 / determinant, result);
    s21_remove_matrix(&temp1);
    s21_remove_matrix(&temp2);
  }
  return status;
}
