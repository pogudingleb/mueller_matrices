#include "main.h"

bool test_pauli(size_t num_tests) {
    bool correct = true;
    double elem_symmeric[4];
    double M[16];
    for (size_t i = 0; i < num_tests; ++i) {
        std::uniform_real_distribution<double> unif(-1., 1.);
        std::default_random_engine re;
        for (size_t j = 0; j < DIM * DIM; ++j) {
            M[j] = unif(re);
        }
        compute_elementary_symmetric(M, elem_symmeric); 
        Eigen::Matrix4cd H;
        build_eigen_matrix(M, H);
        H = 0.5 * H;
        Eigen::SelfAdjointEigenSolver<Eigen::Matrix4cd> es;
        es.compute(H);
        // ugly code...
        double e1 = es.eigenvalues()(0, 0);
        double e2 = es.eigenvalues()(1, 0);
        double e3 = es.eigenvalues()(2, 0);
        double e4 = es.eigenvalues()(3, 0);
        double elem1 = e1 + e2 + e3 + e4;
        double elem2 = e1 * (e2 + e3 + e4) + e2 * (e3 + e4) + e3 * e4;
        double elem3 = e1 * e2 * e3 + e1 * e2 * e4 + e1 * e3 * e4 + e2 * e3 * e4;
        double elem4 = e1 * e2 * e3 * e4;
        double eps = 1e-9;
        if ((std::abs(elem1 - elem_symmeric[0]) > eps) || (std::abs(elem2 - elem_symmeric[1]) > eps) || 
            (std::abs(elem3 - elem_symmeric[2]) > eps) || (std::abs(elem4 - elem_symmeric[3]) > eps)) {
            std::cout << "Problem with the matrix" << std::endl;
            for (size_t i = 0; i < DIM; ++i) {
                for (size_t j = 0; j < DIM; ++j) {
                    std::cout << M[i * DIM + j] << " ";
                }
                std::cout << std::endl;
            }
            std::cout << "Elementary symmetric polynomials via the Pauli matrices:" << std::endl;
            for (size_t i = 0; i < DIM; ++i) {
                std::cout << elem_symmeric[i] << " " << std::endl;
            }
            std::cout << "The ones computed directly from the eigenvalue:" << std::endl;
            std::cout << elem1 << " " << elem2 << " " << elem3 << " " << elem4 << std::endl;
            correct = false;
        }
    }
    return correct;
}

bool test_correctness(double* matrices, size_t num_matrices, bool* result, PSDChecker* test_checker) {
    bool* correct_result = (bool*)malloc(sizeof(bool) * num_matrices);

    run_checker_seq(test_checker, matrices, num_matrices, correct_result);
    bool correct = true;
    for (size_t i = 0; i < num_matrices; ++i) {
        if (result[i] != correct_result[i]) {
            correct = false;
        }
    }

    free(correct_result);
    return correct;
}

int main(int argc, char** argv) {

    double* matrices = read_matrices();
    const size_t NUM_MATRICES = HEIGHT * WIDTH;
    bool* result = (bool*)malloc(sizeof(bool) * NUM_MATRICES);

    for (size_t i = 0; i < checkers.size(); ++i) {
        std::cout << "Checking " << checkers[i].first << " on the data ... ";
        run_checker_seq(checkers[i].second, matrices, NUM_MATRICES, result);
        bool correct = test_correctness(matrices, NUM_MATRICES, result, checkers[0].second);
        if (correct) {
            std::cout << "OK" << std::endl;
        } else {
            std::cout << "Not OK !!!" << std::endl;
        }
    }

    std::cout << "Additional testing for the method based on Pauli matrices ... ";
    if (test_pauli(10000)) {
        std::cout << "OK" << std::endl;
    } else {
        std::cout << "Not OK !!!" << std::endl;
    }

    free(matrices);
    free(result);
    return 0;
}
