#include <chrono>
#include <fstream>
#include <sstream>
#include <iostream>

const size_t N = 600;
const size_t M = 700;
const size_t DIM = 4;

void read_entry(size_t i, size_t j, double* dest) {
    std::string fname = "data/MM_" + std::to_string(i) + std::to_string(j) + ".dat";
    std::ifstream infile(fname);

    std::string line;
    while (getline(infile, line)) {
        std::istringstream ss(line);
        double x;
        size_t counter = i * DIM + j;
        size_t step = DIM * DIM;
        while (ss >> x) {
            dest[counter] = x;
            counter += step;
        }
    }
}

double* read_matrices() {
    size_t total_length = N * M * DIM;
    double* result = (double *)malloc(total_length * sizeof(double));
    for (size_t i = 0; i < DIM; ++i) {
        for (size_t j = 0; j < DIM; ++j) {
            read_entry(i, j, result);
        }
    }
    return result;
}

int main() {
    auto start = std::chrono::system_clock::now();
    double* matrices = read_matrices();
    auto end = std::chrono::system_clock::now();

    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    std::cout << "Matrices read from files in " << elapsed << " milliseconds" << std::endl;
}
