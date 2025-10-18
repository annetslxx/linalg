//
// Created by Anna on 8.10.25.
//

#include "../include/matrix.h"

#include <cmath>
#include <iostream>
#include <stdexcept>
#include <utility>

namespace linalg {
Matrix::Matrix() : m_ptr(nullptr), m_rows(0), m_columns(0), m_capacity(0) {};

Matrix::Matrix(const std::size_t rows)
    : m_ptr(rows == 0 ? nullptr : new double[rows]()), m_rows(rows), m_columns(rows == 0 ? 0 : 1), m_capacity(rows == 0 ? 0 : rows) {
}

Matrix::Matrix(std::size_t rows, std::size_t columns)
    : m_ptr(new double[rows * columns]()), m_rows(rows), m_columns(columns),
      m_capacity(rows * columns) {
  if ((rows == 0 || columns == 0) && (rows != columns))
    throw std::runtime_error("A matrix cannot contain zero number of columns or rows");
}

Matrix::Matrix(const Matrix &other)
    : m_rows(other.rows()), m_columns(other.columns()),
      m_capacity(other.capacity()), m_ptr(new double[other.capacity()]()) {
  std::copy(other.begin(), other.end(), begin());
}

Matrix::Matrix(Matrix &&other) noexcept : m_ptr(nullptr), m_rows(0), m_columns(0), m_capacity(0){
  swap(other);
}
Matrix::Matrix(std::initializer_list<std::initializer_list<double>> list)
    : m_rows(list.size()), m_columns(m_rows == 0 ? 0 : list.begin()->size()),
      m_capacity(m_rows * m_columns),
      m_ptr(m_capacity == 0 ? nullptr : new double[m_capacity]()) {
  for (const auto &row : list) {
    if (row.size() != m_columns) {
      throw std::runtime_error(
          "All rows must have the same size in initializer list");
    }
  }
  if (!empty()) {
    double *iterator = begin();
    for (const auto &row : list) {
      std::copy(row.begin(), row.end(), iterator);
      iterator += m_columns;
    }
  }
}

Matrix::Matrix(std::initializer_list<double> list)
    : m_rows(list.size()), m_columns(m_rows == 0 ? 0 : 1), m_capacity(m_rows),
      m_ptr(m_capacity == 0 ? nullptr : new double[m_capacity]()) {
  if (!empty())
    std::copy(list.begin(), list.end(), begin());
}

Matrix &Matrix::operator=(const Matrix &other) {
  if (this == &other) {
    return *this;
  }

  reshape(other.rows(), other.columns());
  std::copy(other.begin(), other.end(), begin());

  return *this;
}

Matrix &Matrix::operator=(Matrix &&other) noexcept {
  swap(other);
  return *this;
}

void Matrix::reserve(std::size_t number) {
  if (number > capacity()) {
    auto new_ptr = new double[number]();
    std::copy(begin(), end(), new_ptr);

    delete[] m_ptr;
    m_ptr = new_ptr;
    m_capacity = number;
  }
}

void Matrix::reshape(std::size_t rows, std::size_t columns) {
  if ((rows == 0 || columns == 0) && (rows != columns))
    throw std::runtime_error("A matrix cannot contain zero number of columns or rows");
  const std::size_t new_capacity = rows * columns;
  if (new_capacity > capacity())
    reserve(new_capacity);

  m_rows = rows;
  m_columns = columns;
}

void Matrix::clear() noexcept { m_rows = m_columns = 0; }

void Matrix::shrink_to_fit() {
  if (capacity() == size())
    return;

  if (size() == 0) {
    delete[] m_ptr;
    m_ptr = nullptr;
    m_capacity = 0;
    return;
  }

  const size_t new_capacity = size();
  auto new_ptr = new double[new_capacity]();
  std::copy(begin(), end(), new_ptr);

  delete[] m_ptr;
  m_ptr = new_ptr;
  m_capacity = new_capacity;
}

void Matrix::swap(Matrix &other) noexcept {
  std::swap(m_ptr, other.m_ptr);
  std::swap(m_rows, other.m_rows);
  std::swap(m_columns, other.m_columns);
  std::swap(m_capacity, other.m_capacity);
}

double &Matrix::operator()(std::size_t row, std::size_t col) {
  if (row >= rows() || col >= columns()) {
    throw std::runtime_error(
        "Matrix indices out of range"); // Возможно стоит использовать
    // out_of_range
  }
  return begin()[row * columns() + col];
}

const double &Matrix::operator()(std::size_t row, std::size_t col) const {
  if (row >= rows() || col >= columns()) {
    throw std::runtime_error(
        "Matrix indices out of range"); // Возможно стоит использовать
    // out_of_range
  }
  return begin()[row * columns() + col];
}

int get_double_width(double value, std::streamsize precision) {
  int width = 0;
  if (value < 0) {
    width += 1;
    value = -value;
  }

  if (value - floor(value) > Matrix::EPSILON)
    ++width;

  value *= pow(10, precision);

  auto int_part = (long long)value;
  bool has_digit_on_end = false;
  if (int_part == 0) {
    width += 1;
  } else {
    while (int_part > 0) {
      has_digit_on_end = (has_digit_on_end || int_part % 10 != 0);
      width = has_digit_on_end ? width + 1 : width;
      int_part /= 10;
    }
  }

  return width;
}

std::ostream &operator<<(std::ostream &os, const Matrix &matrix) {
  if (matrix.empty()) {
    os << "||";
    return os;
  }

  int max_width = 0;
  int first_col_width = 0;
  int precision = (int)os.precision();

  for (std::size_t j = 0; j < matrix.columns(); ++j) {
    for (std::size_t i = 0; i < matrix.rows(); ++i) {
      int current_width = get_double_width(matrix(i, j), precision);

      if (j == 0) {
        if (current_width > first_col_width) {
          first_col_width = current_width;
        }
      }

      if (current_width > max_width) {
        max_width = current_width;
      }
    }
  }

  os.precision(precision + 1);

  for (std::size_t i = 0; i < matrix.rows(); ++i) {
    os << "|";
    for (std::size_t j = 0; j < matrix.columns(); ++j) {
      if (j == 0)
        os.width(first_col_width);
      else
        os.width(max_width);

      os << matrix(i, j);
      if (j < matrix.columns() - 1) {
        os << " ";
      }
    }
    os << "|";
    if (i < matrix.rows() - 1) {
      os << "\n";
    }
  }

  os.precision(precision);

  return os;
}

Matrix Matrix::operator+() const { return *this; }

Matrix Matrix::operator-() const {
  Matrix result(rows(), columns());
  for (std::size_t i = 0; i < size(); ++i) {
    result.begin()[i] = -begin()[i];
  }
  return result;
}

Matrix &Matrix::operator+=(const Matrix &other) {
  // Проверка размерностей
  if (rows() != other.rows() || columns() != other.columns()) {
    throw std::runtime_error("Matrix dimensions must match for addition");
  }

  for (std::size_t i = 0; i < size(); ++i) {
    begin()[i] += other.begin()[i];
  }

  return *this;
}

Matrix &Matrix::operator-=(const Matrix &other) {
  if (rows() != other.rows() || columns() != other.columns()) {
    throw std::runtime_error("Matrix dimensions must match for subtraction");
  }

  for (std::size_t i = 0; i < size(); ++i) {
    begin()[i] -= other.begin()[i];
  }

  return *this;
}

Matrix &Matrix::operator*=(double value) {
  for (std::size_t i = 0; i < size(); ++i) {
    begin()[i] *= value;
  }
  return *this;
}

Matrix &Matrix::operator*=(const Matrix &other) {
  // A(m×n) * B(n×k) = C(m×k)
  if (columns() != other.rows()) {
    throw std::runtime_error(
        "Matrix dimensions incompatible for multiplication");
  }

  Matrix result(rows(), other.columns());

  for (std::size_t i = 0; i < rows(); ++i) {
    for (std::size_t j = 0; j < other.columns(); ++j) {
      double sum = 0.0;
      for (std::size_t k = 0; k < columns(); ++k) {
        sum += (*this)(i, k) * other(k, j);
      }
      result(i, j) = sum;
    }
  }

  *this = std::move(result);
  return *this;
}

Matrix operator+(const Matrix &left, const Matrix &right) {
  Matrix result(left);
  result += right;
  return result;
}

Matrix operator-(const Matrix &left, const Matrix &right) {
  Matrix result(left);
  result -= right;
  return result;
}

Matrix operator*(const Matrix &left, const Matrix &right) {
  Matrix result(left);
  result *= right;
  return result;
}

Matrix operator*(const Matrix &matrix, double value) {
  Matrix result(matrix);
  result *= value;
  return result;
}

Matrix operator*(double value, const Matrix &matrix) { return matrix * value; }

bool operator==(const Matrix &left, const Matrix &right) {
  // Разные размерности - не равны
  if (left.rows() != right.rows() || left.columns() != right.columns()) {
    return false;
  }

  for (std::size_t i = 0; i < left.size(); ++i) {
    if (!(std::abs(left.begin()[i] - right.begin()[i]) <= Matrix::EPSILON)) {
      return false;
    }
  }

  return true;
}

bool operator!=(const Matrix &left, const Matrix &right) {
  return !(left == right); // Переиспользуем ==
}

double Matrix::norm() const {
  if (empty())
    throw std::runtime_error("Matrix is empty");
  double sum = 0.0;
  for (const auto &x : *this)
    sum += x * x;
  return std::sqrt(sum);
}

double Matrix::trace() const {
  if (empty())
    throw std::runtime_error("Matrix is empty");
  if (rows() != columns())
    throw std::runtime_error("Trace is defined only for square matrices");
  double sum = 0.0;
  for (std::size_t i = 0; i < rows(); ++i)
    sum += (*this)(i, i);
  return sum;
}

Matrix concatenate(const Matrix &left, const Matrix &right) {
  if (left.rows() != right.rows())
    throw std::runtime_error("Matrix must have the same number of rows");
  Matrix result(left.rows(), left.columns() + right.columns());

  auto result_iterator = result.begin();
  for (size_t i = 0; i < result.rows(); ++i) {
    std::copy(left.begin() + i * left.columns(),
              left.begin() + (i + 1) * left.columns(), result_iterator);
    result_iterator += left.columns();
    std::copy(right.begin() + i * right.columns(),
              right.begin() + (i + 1) * right.columns(), result_iterator);
    result_iterator += right.columns();
  }

  return result;
}

void Matrix::swap_rows(std::size_t i, std::size_t j) noexcept {
  if (i == j)
    return;
  for (std::size_t col = 0; col < columns(); ++col)
    std::swap((*this)(i, col), (*this)(j, col));
}

int Matrix::swap_count() const {

  Matrix matrix(*this);
  int count = 0;
  std::size_t pivot_col = 0;

  for (std::size_t pivot_row = 0;
       pivot_row < matrix.rows() && pivot_col < matrix.columns(); ++pivot_row) {
    // Находим максимальный элемент в текущем столбце
    double max_val = std::abs(matrix(pivot_row, pivot_col));
    std::size_t max_row = pivot_row;

    for (std::size_t row = pivot_row + 1; row < matrix.rows(); ++row) {
      if (std::abs(matrix(row, pivot_col)) > max_val) {
        max_val = std::abs(matrix(row, pivot_col));
        max_row = row;
      }
    }

    // Если весь столбец нулевой - переходим к следующему
    if (max_val < Matrix::EPSILON) {
      ++pivot_col;
      --pivot_row; // Компенсируем инкремент в цикле
      continue;
    }
    if (pivot_row != max_row)
      ++count; // увеличиваем счетчик, если макс. элемент не в верхней строке

    ++pivot_col;
  }

  return count;
}

Matrix &Matrix::gauss_forward() {
  if (empty())
    throw std::runtime_error("Matrix must not be empty");
  std::size_t pivot_col = 0;

  for (std::size_t pivot_row = 0; pivot_row < rows() && pivot_col < columns();
       ++pivot_row) {
    // Находим максимальный элемент в текущем столбце
    double max_val = std::abs((*this)(pivot_row, pivot_col));
    std::size_t max_row = pivot_row;

    for (std::size_t row = pivot_row + 1; row < rows(); ++row) {
      if (std::abs((*this)(row, pivot_col)) > max_val) {
        max_val = std::abs((*this)(row, pivot_col));
        max_row = row;
      }
    }

    // Если весь столбец нулевой - переходим к следующему
    if (max_val < EPSILON) {
      ++pivot_col;
      --pivot_row; // Компенсируем инкремент в цикле
      continue;
    }
    swap_rows(pivot_row, max_row);

    // Обнуляем элементы ниже
    for (std::size_t row = pivot_row + 1; row < rows(); ++row) {
      // коэф = текущий эл-т / эл-т в строке, которая поднялась вверх
      double coeff = (*this)(row, pivot_col) / (*this)(pivot_row, pivot_col);
      for (std::size_t col = pivot_col; col < columns(); ++col) {
        (*this)(row, col) -= coeff * (*this)(pivot_row, col);
      }
    }

    ++pivot_col;
       }

  return *this;
}

Matrix &Matrix::gauss_backward() {

  if (empty())
    throw std::runtime_error("Matrix must not be empty");

  for (int i = (int)rows() - 1; i >= 0;
       --i) { // чтобы i могло стать отрицательным
    double diag = (*this)(i, i);
    if (std::fabs(diag) < EPSILON)
      throw std::runtime_error("Zero on diagonal");

    // Нормируем текущую строку
    for (std::size_t j = 0; j < columns(); ++j)
      (*this)(i, j) /= diag;

    // Вычитаем текущую строку из всех выше стоящих
    for (int k = i - 1; k >= 0; --k) {
      double factor = (*this)(k, i);
      for (std::size_t j = 0; j < columns();
           ++j) // из каждого элемента строки вычитаем элемент строки,
             // где опорный элемент стал единицей
               (*this)(k, j) -= factor * (*this)(i, j);
    }
       }
  return *this;
}

double Matrix::det() const {
  if (empty())
    throw std::runtime_error("Matrix is empty");
  if (rows() != columns())
    throw std::runtime_error("Matrix must be square");
  Matrix temp = *this;
  temp.gauss_forward();
  double det = 1.0;
  for (std::size_t i = 0; i < temp.rows(); ++i)
    det *= temp(i, i);
  return det * std::pow(-1, swap_count());
}

size_t Matrix::rank() const {
  Matrix temp = *this;
  temp.gauss_forward();
  size_t rank = 0;

  bool is_not_null;
  for (size_t row = 0; row < temp.rows(); ++row) {
    is_not_null = false; // считаем, что строка нулевая
    for (size_t col = 0; col < temp.columns() && !is_not_null; ++col) {
      is_not_null = std::fabs(temp(row, col)) > EPSILON;
    }
    rank = is_not_null ? rank + 1 : rank;
  }
  return rank;
}

Matrix solve(const Matrix &matr_A, const Matrix &vec_f) {
  if (matr_A.rank() != vec_f.rows())
    throw std::runtime_error("The system has either no solution or infinitely many solutions");
  Matrix full = concatenate(matr_A, vec_f).gauss_forward().gauss_backward();
  Matrix solution(vec_f.rows(), vec_f.columns());
  for (std::size_t i = 0; i < full.rows(); ++i) {
    solution(i, 0) = full(i, full.columns() - 1);
  }

  return solution;
}

Matrix transpose(const Matrix &matrix) {
  Matrix ret(matrix.columns(), matrix.rows());
  for (std::size_t row = 0; row < matrix.rows(); ++row) {
    for (std::size_t column = 0; column < matrix.columns(); ++column) {
      ret(column, row) = matrix(row, column);
    }
  }
  return ret;
}

Matrix create_identity_matrix(std::size_t length) {
  Matrix ret(length, length);
  for (std::size_t i = 0; i < length; i++) {
    ret(i, i) = 1;
  }
  return ret;
}

Matrix invert(const Matrix &matrix) {
  if (matrix.rows() != matrix.columns())
    throw std::runtime_error("Matrix must be square");
  Matrix full = concatenate(matrix, create_identity_matrix(matrix.rows()))
                    .gauss_forward()
                    .gauss_backward();
  Matrix ret(matrix.rows(), matrix.columns());

  for (std::size_t row = 0; row < ret.rows(); ++row) {
    for (std::size_t column = 0; column < ret.columns(); ++column) {
      ret(row, column) = full(row, column + matrix.columns());
    }
  }

  return ret;
}

Matrix power(const Matrix &matrix, int power) {
  if (matrix.rows() != matrix.columns())
    throw std::runtime_error("Matrix must be square");

  if (power == 0)
    return create_identity_matrix(matrix.rows());

  Matrix ret = power < 0 ? invert(matrix) : matrix;
  power = std::abs(power);
  Matrix mlt(ret); // для того, чтобы если матрица стала обратной, умножать на
                   // нее, а не на исходную
  while (--power)
    ret *= mlt;
  return ret;
}

} // namespace linalg