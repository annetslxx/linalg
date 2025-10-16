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
      m_ptr(m_capacity == 0 ? nullptr : new double[m_capacity]()) {
  for (const auto &row : list) {
    if (row.size() != m_columns) {
      throw std::runtime_error(
          "All rows must have the same size in initializer list");
    }
  }
  if (!empty()) {
    double *iterator = this->begin();
    for (const auto &row : list) {
      std::copy(row.begin(), row.end(), iterator);
      iterator += m_columns;
    }
  }
}

Matrix::Matrix(std::initializer_list<double> list)
    : m_rows(list.size()), m_columns(m_rows == 0 ? 0 : 1), m_capacity(m_rows),
      m_ptr(m_capacity == 0 ? nullptr : new double[m_capacity]()) {
  if (!(this->empty()))
    std::copy(list.begin(), list.end(), this->begin());
}

// TODO add equality checking after == implemented
Matrix &Matrix::operator=(const Matrix &other) {
  if (this == &other) {
    return *this;
  }

  if (this->capacity() < other.size()) {
    auto new_ptr = new double[other.size()];

    std::copy(other.begin(), other.end(), new_ptr);

    if (!(this->empty())) {
      delete[] this->m_ptr;
    }

    this->m_ptr = new_ptr;
    this->m_capacity = other.size();
  } else {
    std::copy(other.begin(), other.end(), this->begin());
  }

  this->m_rows = other.rows();
  this->m_columns = other.columns();

  return *this;
}

Matrix &Matrix::operator=(Matrix &&other) noexcept {

  delete[] m_ptr;

  this->m_ptr = other.m_ptr;
  this->m_rows = other.rows();
  this->m_columns = other.columns();
  this->m_capacity = other.capacity();

  other.m_ptr = nullptr;
  other.m_rows = 0;
  other.m_columns = 0;
  other.m_capacity = 0;

  return *this;
}

void Matrix::reshape(std::size_t rows, std::size_t columns) {

  this->m_rows = rows;
  this->m_columns = columns;

  const std::size_t new_capacity = this->size();

  if (new_capacity > this->m_capacity) {
    auto new_ptr = new double[rows * columns]();

    delete[] this->m_ptr;
    this->m_ptr = new_ptr;
    this->m_capacity = new_capacity;
  }
}

void Matrix::reserve(std::size_t number) {
  if (number > this->capacity()) {
    auto new_ptr = new double[number]();

    delete[] this->m_ptr;
    this->m_ptr = new_ptr;
    this->m_capacity = number;
  }
}

void Matrix::clear() { m_rows = m_columns = 0; }

void Matrix::shrink_to_fit() {
  if (this->capacity() == this->size()) {
    return;
  }

  if (this->size() == 0) {
    delete[] this->m_ptr;
    this->m_ptr = nullptr;
    this->m_capacity = 0;
    return;
  }

  auto new_ptr = new double[this->size()]();
  std::copy(this->begin(), this->end(), new_ptr);

  delete[] this->m_ptr;
  this->m_ptr = new_ptr;
  this->m_capacity = this->size();
}

void Matrix::swap(Matrix &other) {
  std::swap(this->m_ptr, other.m_ptr);
  std::swap(this->m_rows, other.m_rows);
  std::swap(this->m_columns, other.m_columns);
  std::swap(this->m_capacity, other.m_capacity);
}
} // namespace linalg
