//
// Created by Anna on 8.10.25.
//

#include "matrix.h"
namespace linalg {

Matrix::Matrix() : m_ptr(nullptr), m_rows(0), m_columns(0), m_capacity(0) {};

Matrix::Matrix(const std::size_t rows)
    : m_ptr(new double[rows]()), m_rows(rows), m_columns(1), m_capacity(rows) {}

Matrix::Matrix(std::size_t rows, std::size_t columns)
    : m_ptr(new double[rows * columns]()), m_rows(rows), m_columns(columns),
      m_capacity(rows * columns) {}

Matrix::Matrix(const Matrix &other)
    : m_rows(other.rows()), m_columns(other.columns()),
      m_capacity(other.capacity()), m_ptr(new double[other.capacity()]) {

  std::copy(other.begin(), other.end(), begin());
}

Matrix::Matrix(Matrix &&other) noexcept
    : m_ptr(other.begin()), m_rows(other.rows()), m_columns(other.columns()),
      m_capacity(other.capacity()) {

  other.m_ptr = nullptr;
  other.m_rows = 0;
  other.m_columns = 0;
  other.m_capacity = 0;
}

Matrix::Matrix(std::initializer_list<std::initializer_list<double>> list)
    : m_rows(list.size()), m_columns(m_rows == 0 ? 0 : list.begin()->size()),
      m_capacity(m_rows * m_columns),
      m_ptr(m_capacity == 0 ? nullptr : new double[m_capacity]) {
  for (const auto &row : list) {
    if (row.size() != m_columns) {
      throw std::runtime_error(
          "All rows must have the same size in initializer list");
    }
  }
  if (!empty()) {
    double *dest = begin();
    for (const auto &row : list) {
      std::copy(row.begin(), row.end(), dest);
      dest += m_columns;
    }
  }
}

Matrix::Matrix(std::initializer_list<double> list)
    : m_rows(list.size()), m_columns(m_rows == 0 ? 0 : 1), m_capacity(m_rows),
      m_ptr(m_capacity == 0 ? nullptr : new double[m_capacity]) {
  if (empty())
    return;

  std::copy(list.begin(), list.end(), begin());
}

} // namespace linalg
