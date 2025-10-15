//
// Created by Anna Dolgopolova on 29.09.2025.
//
#ifndef LINALG_MATRIX_H
#define LINALG_MATRIX_H

#include <cstddef>
#include <initializer_list>
#include <stdexcept>

namespace linalg {

class Matrix {

private:
  double *m_ptr;
  std::size_t m_rows;
  std::size_t m_columns;
  std::size_t m_capacity;

  static constexpr double GROWTH_FACTOR = 2.0;

public:
  Matrix();
  explicit Matrix(std::size_t rows);
  explicit Matrix(std::size_t rows, std::size_t columns);

  Matrix(const Matrix &other);
  Matrix(Matrix &&other);

  Matrix(std::initializer_list<std::initializer_list<double>> list);
  Matrix(std::initializer_list<double> list);

  ~Matrix() { delete[] this->m_ptr; }

  Matrix &operator=(const Matrix &other);
  Matrix &operator=(Matrix &&other);

  double &operator()(std::size_t row, std::size_t col);
  const double &operator()(std::size_t row, std::size_t col) const;

  Matrix operator+() const;
  Matrix operator-() const;

  Matrix &operator+=(const Matrix &other);
  Matrix &operator-=(const Matrix &other);
  Matrix &operator*=(const Matrix &other);
  Matrix &operator*=(double value);

  bool empty() const { return this->m_rows == 0 || this->m_columns == 0; }

  std::size_t rows() const { return this->m_rows; }
  std::size_t columns() const { return this->m_columns; }
  std::size_t capacity() const { return this->m_capacity; }
  std::size_t size() const { return this->m_rows * this->m_columns; }

  double *begin() { return this->m_ptr; }
  double *end() { return this->m_ptr + this->size(); }
  const double *begin() const { return this->m_ptr; }
  const double *end() const { return this->m_ptr + this->size(); }

  void reshape(std::size_t rows, std::size_t columns);
  void reserve(std::size_t number);
  void clear();
  void shrink_to_fit();
  void swap(Matrix &destination);

  friend Matrix operator+(const Matrix &left, const Matrix &right);
  friend Matrix operator-(const Matrix &left, const Matrix &right);
  friend Matrix operator*(const Matrix &left, const Matrix &right);
  friend Matrix operator*(const Matrix &left, double value);
  friend Matrix operator*(double value, const Matrix &right);

  friend bool operator==(const Matrix &left, const Matrix &right);
  friend bool operator!=(const Matrix &left, const Matrix &right);
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
} // namespace linalg

#endif // LINALG_MATRIX_H
