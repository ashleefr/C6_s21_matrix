#ifndef SRC_S21_MATRIX_H_
#define SRC_S21_MATRIX_H_

//// LIBRARIES ////

#include <math.h>
#include <stdlib.h>

//// CONSTANTS ////

#define EPS 1.e-7

#define OK 0
#define INCORRECT_MATRIX 1
#define ERROR_CALCULATING 2

#define SUCCESS 1
#define FAILURE 0

//// STRUCTURE ////

typedef struct matrix_struct {
  double **matrix;
  int rows;
  int columns;
} matrix_t;

//// MAIN FUNCTIONS ////

int s21_create_matrix(int rows, int columns, matrix_t *result);
void s21_remove_matrix(matrix_t *A);
int s21_eq_matrix(matrix_t *A, matrix_t *B);
int s21_sum_matrix(matrix_t *A, matrix_t *B, matrix_t *result);
int s21_sub_matrix(matrix_t *A, matrix_t *B, matrix_t *result);
int s21_mult_number(matrix_t *A, double number, matrix_t *result);
int s21_mult_matrix(matrix_t *A, matrix_t *B, matrix_t *result);
int s21_transpose(matrix_t *A, matrix_t *result);
int s21_calc_complements(matrix_t *A, matrix_t *result);
int s21_determinant(matrix_t *A, double *result);
int s21_inverse_matrix(matrix_t *A, matrix_t *result);

//// SUPPORT FUNCTIONS ////

int check_incorrect_matrix(matrix_t *target);
int check_valid_eq_matrix(matrix_t *suspectA, matrix_t *suspectB);
int check_valid_arithmetic_matrix(matrix_t *A, matrix_t *B, matrix_t *result);
int check_valid_multiplication_matrix(matrix_t *A, matrix_t *B,
                                      matrix_t *result);
double find_determinant(matrix_t temp);
void swap_rows(matrix_t *temp, int i, int k);
void create_minor(int erase_row, int erase_col, matrix_t *temp, matrix_t *A);

#endif  // SRC_S21_MATRIX_H_
