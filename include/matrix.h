//
// Created by Anna Dolgopolova on 29.09.2025.
//
#ifndef LINALG_MATRIX_H
#define LINALG_MATRIX_H

#include <cstddef>
#include <initializer_list>
#include <ostream>

namespace linalg {

class Matrix {

private:
  std::size_t m_rows;
  std::size_t m_columns;
  std::size_t m_capacity;
  double *m_ptr;

  class Row {
  private:
    double* m_row_ptr;
    std::size_t m_columns; // количество элементов в строке (то же самое, что кол-во столбцов), для проверки границ
  public:
    Row(double* row_ptr, std::size_t columns) noexcept;
    double& operator[](std::size_t col);
    const double& operator[](std::size_t col) const;
  };

public:
  inline static double EPSILON = 1e-9;

  Matrix();
  explicit Matrix(std::size_t rows);
  explicit Matrix(std::size_t rows, std::size_t columns);

  Matrix(const Matrix &other);
  Matrix(Matrix &&other) noexcept;

  Matrix(std::initializer_list<std::initializer_list<double>> list);
  Matrix(std::initializer_list<double> list);

  ~Matrix() { delete[] m_ptr; }

  Matrix &operator=(const Matrix &other);
  Matrix &operator=(Matrix &&other) noexcept;

  double &operator()(std::size_t row, std::size_t col);
  const double &operator()(std::size_t row, std::size_t col) const;

  Row operator[](std::size_t row);
  const Row operator[](std::size_t row) const;

  Matrix operator+() const;
  Matrix operator-() const;

  Matrix &operator+=(const Matrix &other);
  Matrix &operator-=(const Matrix &other);
  Matrix &operator*=(const Matrix &other);
  Matrix &operator*=(double value);

  bool empty() const noexcept { return size() == 0; }

  std::size_t rows() const noexcept { return m_rows; }
  std::size_t columns() const noexcept { return m_columns; }
  std::size_t capacity() const noexcept { return m_capacity; }
  std::size_t size() const noexcept { return m_rows * m_columns; }

  double *begin() noexcept { return m_ptr; }
  double *end() noexcept { return m_ptr + size(); }
  const double *begin() const noexcept { return m_ptr; }
  const double *end() const noexcept { return m_ptr + size(); }

  void reshape(std::size_t rows, std::size_t columns);
  void reserve(std::size_t number);
  void clear() noexcept;
  void shrink_to_fit();
  void swap(Matrix &other) noexcept;

  double norm() const;
  double trace() const;
  double det() const;
  size_t rank() const;
  void swap_rows(std::size_t i, std::size_t j) noexcept;
  Matrix &gauss_forward();
  Matrix &gauss_backward();

  int swap_count() const;
};

inline void swap(Matrix &left_matrix, Matrix &right_matrix) {
  left_matrix.swap(right_matrix);
}

Matrix operator+(const Matrix &left, const Matrix &right);
Matrix operator-(const Matrix &left, const Matrix &right);
Matrix operator*(const Matrix &left, const Matrix &right);
Matrix operator*(const Matrix &left, double value);
Matrix operator*(double value, const Matrix &right);

bool operator==(const Matrix &left, const Matrix &right);
bool operator!=(const Matrix &left, const Matrix &right);

std::ostream &operator<<(std::ostream &os, const Matrix &matrix);

Matrix transpose(const Matrix &matrix);
Matrix concatenate(const Matrix &left, const Matrix &right);
Matrix solve(const Matrix &matr_A, const Matrix &vec_f);
Matrix power(const Matrix &matr, int power);
Matrix create_identity_matrix(std::size_t length);
Matrix invert(const Matrix &matrix);

} // namespace linalg

#endif // LINALG_MATRIX_H