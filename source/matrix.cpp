//
// Created by Anna on 8.10.25.
//

#include "matrix.h"
namespace linalg {

Matrix::Matrix() : m_ptr(nullptr), m_rows(0), m_columns(0), m_capacity(0) {};
Matrix::Matrix(std::size_t rows)
    : m_ptr(new double[rows]()), m_rows(rows), m_columns(1), m_capacity(rows) {};
Matrix::Matrix(std::size_t rows, std::size_t columns)
    : m_ptr(new double[rows * columns]()), m_rows(rows), m_columns(columns),
      m_capacity(rows * columns) {};

} // namespace linalg
