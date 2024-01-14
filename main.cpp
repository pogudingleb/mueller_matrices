#include "main.h"

int main(int argc, char** argv) {

    if (argc != 3) {
        std::cout << "Wrong command line arguments. The function should be called as:" << std::endl;
        std::cout << "\t ./matrices method num_threads" << std::endl;
        std::cout << "Where num_threads is the number of threads to be used and method is the index of the method. "
                  << "For the latter, the available options are:" << std::endl;
        for (size_t i = 0; i < checkers.size(); ++i) {
            std::cout << "\t" << i << " : " << checkers[i].first << std::endl;
        }
        return 0;
    }

    size_t checker_index = std::atoi(argv[1]);
    PSDChecker* checker_fun = checkers[checker_index].second;
    size_t num_threads = std::atoi(argv[2]);

    auto start = std::chrono::system_clock::now();
    double* matrices = read_matrices();
    auto end = std::chrono::system_clock::now();

    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    std::cout << "Matrices read from files in " << elapsed << " milliseconds" << std::endl;

    const size_t NUM_MATRICES = HEIGHT * WIDTH;
    bool* result = (bool*)malloc(sizeof(bool) * NUM_MATRICES);
    std::cout << "Checking PSD using " << checkers[checker_index].first << std::endl; 
    start = std::chrono::system_clock::now();
    if (num_threads == 1) {
        run_checker_seq(checker_fun, matrices, NUM_MATRICES, result);
    } else {
        run_checker_parallel(checker_fun, matrices, NUM_MATRICES, result, num_threads);
    }
    end = std::chrono::system_clock::now();
    elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    std::cout << "Checking performed in " << elapsed << " milliseconds" << std::endl;

    size_t nonphysical_count = 0;
    for (size_t i = 0; i < NUM_MATRICES; ++i) {
        if (!result[i]) {
            nonphysical_count++;
        }
    }
    std::cout << "The number of non-physical matrices is " << nonphysical_count << std::endl;

    free(matrices);
    free(result);
    return 0;
}
