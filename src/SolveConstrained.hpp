#pragma once

#include <vector>
#include <Eigen/Sparse>
#include <Eigen/Dense>


void solveConstrainedSymmetric(const Eigen::SparseMatrix<double>& A, const Eigen::MatrixXd& B,
                               const std::vector<int>& constr, const Eigen::MatrixXd& constrValues,
                               Eigen::MatrixXd& X);

