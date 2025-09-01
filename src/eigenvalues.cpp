#include "main.h"

void compute_charpoly_coefficients(double* matrices, size_t num_matrices, double* result) {
    const size_t STEP = DIM * DIM;
    for (size_t i = 0; i < num_matrices; ++i) {
        compute_elementary_symmetric(matrices + STEP * i, result + DIM * i);
        result[DIM * i] *= -1.;
        result[DIM * i + 2] *= -1;
    }
}

void ferrari_massive(double* coefficients, size_t num_matrices, double* result) {
    for (size_t i = 0; i < num_matrices; ++i) {
        ferrari_method(coefficients + DIM * i, result + DIM * i);
    }
}

int main(int argc, char** argv) {

    if (argc != 2) {
        std::cout << "Wrong command line arguments. The function should be called as:" << std::endl;
        std::cout << "\t ./eigenvalues num_threads" << std::endl;
        std::cout << "Where num_threads is the number of threads to be used." << std::endl;
        return 0;
    }

    size_t num_threads = std::atoi(argv[1]);

    auto start = std::chrono::system_clock::now();
    double* matrices = read_matrices();
    auto end = std::chrono::system_clock::now();

    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    std::cout << "Matrices read from files in " << elapsed << " milliseconds" << std::endl;

    const size_t NUM_MATRICES = HEIGHT * WIDTH;
    const size_t STEP = DIM * DIM;
    double* coefficients = (double*)malloc(DIM * sizeof(double) * NUM_MATRICES);

    //----------------------

    std::cout << "Computing coefficients of the charpoly via the Pauli matrices method" << std::endl; 
    start = std::chrono::system_clock::now();

    std::vector<std::thread> workers(num_threads - 1);
    size_t chunk_size = NUM_MATRICES / num_threads;
    size_t offset = 0;
    for (size_t i = 0; i < num_threads - 1; ++i) {
        workers[i] = std::thread(&compute_charpoly_coefficients, matrices + offset * STEP, chunk_size, coefficients + offset * DIM);
        offset += chunk_size;
    }
    compute_charpoly_coefficients(matrices + offset * STEP, NUM_MATRICES - chunk_size * (num_threads - 1), coefficients + offset * DIM);
    for (auto &t: workers) {
        t.join();
    }

    end = std::chrono::system_clock::now();
    std::cout << "Coefficients of the charpolys computed in " 
        << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << " milliseconds" << std::endl;

    //----------------------

    double* eigenvalues = (double*)malloc(DIM * sizeof(double) * NUM_MATRICES);
    start = std::chrono::system_clock::now();

    std::vector<std::thread> eigenworkers(num_threads - 1);
    offset = 0;
    for (size_t i = 0; i < num_threads - 1; ++i) {
        eigenworkers[i] = std::thread(&ferrari_massive, coefficients + offset * DIM, chunk_size, eigenvalues + offset * DIM);
        offset += chunk_size;
    }
    ferrari_massive(coefficients + offset * DIM, NUM_MATRICES - chunk_size * (num_threads - 1), eigenvalues + offset * DIM);
    for (auto &t: eigenworkers) {
        t.join();
    }

    end = std::chrono::system_clock::now();
    std::cout << "Eigenvalues computed in " 
        << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << " milliseconds" << std::endl;

    //----------------------

    size_t count_negative = 0;
    for (size_t i = 0; i < NUM_MATRICES; ++i) {
        bool neg_value = (eigenvalues[i * DIM] < 0.);
        bool nonrealizable = !((coefficients[i * DIM] <= 0.) && (coefficients[i * DIM + 1] >= 0.) && (coefficients[i * DIM + 2] <= 0.) && (coefficients[i * DIM + 3] >= 0.));
        if (neg_value != nonrealizable) {
            std::cout << "Problem at " << i << std::endl;
            std::cout << "Coefficients ";
            for (size_t j = 0; j < DIM; ++j) {
                std::cout << coefficients[i * DIM + j] << " ";
            }
            std::cout << std::endl;
            std::cout << "Real parts of roots ";
            for (size_t j = 0; j < DIM; ++j) {
                std::cout << eigenvalues[i * DIM + j] << " ";
            }
            std::cout << std::endl << std::endl;
        }
        if (eigenvalues[i * DIM] <= 0.) {
            ++count_negative;
        }
    }
    std::cout << "There are " << count_negative << " physically nonrealizable matrices" << std::endl;

    //----------------------

    std::cout << "Comparing against Eigen's solver" << std::endl;
    double diff = 0.;
    for (size_t i = 0; i < NUM_MATRICES; ++i) {
        Eigen::Matrix4cd H;
        build_eigen_matrix(matrices + STEP * i, H);
        Eigen::SelfAdjointEigenSolver<Eigen::Matrix4cd> es;
        es.compute(H);
        // eigenvalues are sorted
        for (size_t j = 0; j < DIM; ++j) {
            // the eigenvalue we compute are twice smaller that the ones of the H-matrix !
            diff = std::max(diff, std::abs(es.eigenvalues()(j, 0) - 2.0 * eigenvalues[i * DIM + j]));
        }
    }
    std::cout << "The maximal discrepancy is " << diff << std::endl;


    free(matrices);
    free(coefficients);
    free(eigenvalues);
    return 0;
}
