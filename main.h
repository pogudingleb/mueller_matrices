#pragma once

#include <cstddef>
#include <cstdint>
#include <cwctype>
#include <algorithm>
#include <array>
#include <chrono>
#include <iostream>
#include <cstdio>

const size_t WIDTH = 700;
const size_t HEIGHT = 600;
const size_t DIM = 4;
const size_t MAX_PRECISION = 25;

double* read_matrices();
