#pragma once

#include <cstddef>
#include <cstdint>
#include <cwctype>
#include <algorithm>
#include <array>
#include <chrono>
#include <iostream>
#include <cstdio>
#include <random>
#include <thread>
#include <vector>

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

const size_t WIDTH = 700;
const size_t HEIGHT = 600;
const size_t DIM = 4;
const size_t MAX_PRECISION = 25;

typedef bool PSDChecker(const double*);

// auxiliary functions
double* read_matrices();
void build_eigen_matrix(const double* M, Eigen::Matrix4cd& H);

// available PSD checkers
bool check_with_charpoly_noh(const double* M);
bool check_with_sylvester_noh(const double* M);
bool check_eignevalues(const double* M);
bool check_choletsky(const double* M);
void compute_elementary_symmetric(const double* M, double* result); // helper function for the Pauli checker
bool check_pauli(const double* M);

const std::vector<std::pair<std::string, PSDChecker*> > checkers{
    std::make_pair(std::string("Charpoly avoiding H"), &check_with_charpoly_noh),
    std::make_pair(std::string("Sylvester criterion avoiding H"), &check_with_sylvester_noh),
    std::make_pair(std::string("Eigen's self-adjoint eigensolver"), &check_eignevalues),
    std::make_pair(std::string("Eigen's Choletsky"), &check_choletsky),
    std::make_pair(std::string("Formula via Pauli matrices"), &check_pauli)
};

// generic functions to run checkers on data
void run_checker_seq(PSDChecker* checker, double* matrices, size_t num_matrices, bool* result);
void run_checker_parallel(PSDChecker* checker, double* matrices, size_t num_matrices, bool* result, size_t num_threads);
