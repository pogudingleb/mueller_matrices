#include "main.h"

void run_checker_seq(PSDChecker* checker, double* matrices, size_t num_matrices, bool* result) {
    size_t step = DIM * DIM;
    for (size_t i = 0; i < num_matrices; ++i) {
        result[i] = (*checker)(matrices + step * i);
    }
}

void run_checker_parallel(PSDChecker* checker, double* matrices, size_t num_matrices, bool* result, size_t num_threads) {
    size_t chunk_size = num_matrices / num_threads;
    std::vector<std::thread> workers(num_threads - 1); 
    double* matrices_local = matrices;
    bool* result_local = result;

    for (size_t i = 0; i < num_threads - 1; ++i) {
        workers[i] = std::thread(&run_checker_seq, checker, matrices_local, chunk_size, result_local);
        matrices_local += chunk_size * DIM * DIM;
        result_local += chunk_size;
    }
    run_checker_seq(checker, matrices_local, num_matrices - (num_threads - 1) * chunk_size, result_local);

    for (size_t i = 0; i < num_threads - 1; ++i) {
        workers[i].join();
    }
}
