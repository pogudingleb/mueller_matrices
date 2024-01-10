#include "main.h"

// ------------
// PSD Checkers
// ------------

// Some common low-level functions
template <size_t i, size_t j>
inline double& get(double* M) {
    return M[i * DIM + j];
}

inline double pow2(double x) {
    return x * x;
}

inline double pow3(double x) {
    return x * x * x;
}

inline double pow4(double x) {
    double y = x * x;
    return y * y;
}

// Checker via the charpoly

bool check_with_charpoly_noh(double* M) {
    double C1 = get<0, 0>(M);

    double C2 = 3 * pow2(get<0 ,0>(M)) - pow2(get<0 ,1>(M)) - pow2(get<0 ,2>(M)) - pow2(get<0 ,3>(M)) - pow2(get<1 ,0>(M)) - pow2(get<1 ,1>(M)) - pow2(get<1 ,2>(M)) - pow2(get<1 ,3>(M)) - pow2(get<2 ,0>(M)) - pow2(get<2 ,1>(M)) - pow2(get<2 ,2>(M)) - pow2(get<2 ,3>(M)) - pow2(get<3 ,0>(M)) - pow2(get<3 ,1>(M)) - pow2(get<3 ,2>(M)) - pow2(get<3 ,3>(M));

    double C3 = 4*pow3(get<0 ,0>(M))+(-4*pow2(get<0 ,1>(M))-4*pow2(get<0 ,2>(M))-4*pow2(get<0 ,3>(M))-4*pow2(get<1 ,0>(M))-4*pow2(get<1 ,1>(M))-4*pow2(get<1 ,2>(M))-4*pow2(get<1 ,3>(M))-4*pow2(get<2 ,0>(M))-4*pow2(get<2 ,1>(M))-4*pow2(get<2 ,2>(M))-4*pow2(get<2 ,3>(M))-4*pow2(get<3 ,0>(M))-4*pow2(get<3 ,1>(M))-4*pow2(get<3 ,2>(M))-4*pow2(get<3 ,3>(M)))*get<0 ,0>(M)+(8*get<1 ,0>(M)*get<1 ,1>(M)+8*get<2 ,0>(M)*get<2 ,1>(M)+8*get<3 ,0>(M)*get<3 ,1>(M))*get<0 ,1>(M)+(8*get<1 ,0>(M)*get<1 ,2>(M)+8*get<2 ,0>(M)*get<2 ,2>(M)+8*get<3 ,0>(M)*get<3 ,2>(M))*get<0 ,2>(M)+(8*get<1 ,0>(M)*get<1 ,3>(M)+8*get<2 ,0>(M)*get<2 ,3>(M)+8*get<3 ,0>(M)*get<3 ,3>(M))*get<0 ,3>(M)+(8*get<2 ,2>(M)*get<3 ,3>(M)-8*get<2 ,3>(M)*get<3 ,2>(M))*get<1 ,1>(M)+(-8*get<2 ,1>(M)*get<3 ,3>(M)+8*get<2 ,3>(M)*get<3 ,1>(M))*get<1 ,2>(M)+8*get<1 ,3>(M)*(get<2 ,1>(M)*get<3 ,2>(M)-get<2 ,2>(M)*get<3 ,1>(M));

    double C4 = pow4(get<0 ,0>(M))+pow4(get<0 ,1>(M))+pow4(get<0 ,2>(M))+pow4(get<1 ,0>(M))+pow4(get<1 ,1>(M))+pow4(get<1 ,2>(M))+pow4(get<2 ,0>(M))+pow4(get<2 ,1>(M))+pow4(get<2 ,2>(M))+pow4(get<0 ,3>(M))+pow4(get<1 ,3>(M))+pow4(get<2 ,3>(M))+8*get<2 ,2>(M)*get<2 ,3>(M)*get<3 ,2>(M)*get<3 ,3>(M)-8*get<3 ,0>(M)*(get<2 ,1>(M)*get<3 ,1>(M)+get<2 ,2>(M)*get<3 ,2>(M)+get<2 ,3>(M)*get<3 ,3>(M))*get<2 ,0>(M)+8*get<3 ,1>(M)*(get<2 ,2>(M)*get<3 ,2>(M)+get<2 ,3>(M)*get<3 ,3>(M))*get<2 ,1>(M)+8*get<1 ,3>(M)*(get<2 ,2>(M)*get<2 ,3>(M)+get<3 ,2>(M)*get<3 ,3>(M))*get<1 ,2>(M)+pow2(pow2(get<3 ,0>(M))-pow2(get<3 ,1>(M))-pow2(get<3 ,2>(M))-pow2(get<3 ,3>(M)))+(2*pow2(get<3 ,0>(M))-2*pow2(get<3 ,1>(M))-2*pow2(get<3 ,2>(M))+2*pow2(get<3 ,3>(M)))*pow2(get<2 ,3>(M))+(2*pow2(get<2 ,3>(M))+2*pow2(get<3 ,0>(M))-2*pow2(get<3 ,1>(M))+2*pow2(get<3 ,2>(M))-2*pow2(get<3 ,3>(M)))*pow2(get<2 ,2>(M))+(2*pow2(get<2 ,2>(M))+2*pow2(get<2 ,3>(M))+2*pow2(get<3 ,0>(M))+2*pow2(get<3 ,1>(M))-2*pow2(get<3 ,2>(M))-2*pow2(get<3 ,3>(M)))*pow2(get<2 ,1>(M))+(-2*pow2(get<2 ,1>(M))-2*pow2(get<2 ,2>(M))-2*pow2(get<2 ,3>(M))+2*pow2(get<3 ,0>(M))+2*pow2(get<3 ,1>(M))+2*pow2(get<3 ,2>(M))+2*pow2(get<3 ,3>(M)))*pow2(get<2 ,0>(M))+(2*pow2(get<2 ,0>(M))-2*pow2(get<2 ,1>(M))-2*pow2(get<2 ,2>(M))+2*pow2(get<2 ,3>(M))+2*pow2(get<3 ,0>(M))-2*pow2(get<3 ,1>(M))-2*pow2(get<3 ,2>(M))+2*pow2(get<3 ,3>(M)))*pow2(get<1 ,3>(M))+(2*pow2(get<1 ,3>(M))+2*pow2(get<2 ,0>(M))-2*pow2(get<2 ,1>(M))+2*pow2(get<2 ,2>(M))-2*pow2(get<2 ,3>(M))+2*pow2(get<3 ,0>(M))-2*pow2(get<3 ,1>(M))+2*pow2(get<3 ,2>(M))-2*pow2(get<3 ,3>(M)))*pow2(get<1 ,2>(M))+((8*get<2 ,1>(M)*get<2 ,2>(M)+8*get<3 ,1>(M)*get<3 ,2>(M))*get<1 ,2>(M)+8*get<1 ,3>(M)*(get<2 ,1>(M)*get<2 ,3>(M)+get<3 ,1>(M)*get<3 ,3>(M)))*get<1 ,1>(M)+(2*pow2(get<1 ,2>(M))+2*pow2(get<1 ,3>(M))+2*pow2(get<2 ,0>(M))+2*pow2(get<2 ,1>(M))-2*pow2(get<2 ,2>(M))-2*pow2(get<2 ,3>(M))+2*pow2(get<3 ,0>(M))+2*pow2(get<3 ,1>(M))-2*pow2(get<3 ,2>(M))-2*pow2(get<3 ,3>(M)))*pow2(get<1 ,1>(M))+((-8*get<2 ,0>(M)*get<2 ,1>(M)-8*get<3 ,0>(M)*get<3 ,1>(M))*get<1 ,1>(M)+(-8*get<2 ,0>(M)*get<2 ,2>(M)-8*get<3 ,0>(M)*get<3 ,2>(M))*get<1 ,2>(M)-8*get<1 ,3>(M)*(get<2 ,0>(M)*get<2 ,3>(M)+get<3 ,0>(M)*get<3 ,3>(M)))*get<1 ,0>(M)+(-2*pow2(get<1 ,1>(M))-2*pow2(get<1 ,2>(M))-2*pow2(get<1 ,3>(M))+2*pow2(get<2 ,0>(M))+2*pow2(get<2 ,1>(M))+2*pow2(get<2 ,2>(M))+2*pow2(get<2 ,3>(M))+2*pow2(get<3 ,0>(M))+2*pow2(get<3 ,1>(M))+2*pow2(get<3 ,2>(M))+2*pow2(get<3 ,3>(M)))*pow2(get<1 ,0>(M))+((-8*get<2 ,1>(M)*get<3 ,2>(M)+8*get<2 ,2>(M)*get<3 ,1>(M))*get<1 ,0>(M)+(8*get<2 ,0>(M)*get<3 ,2>(M)-8*get<2 ,2>(M)*get<3 ,0>(M))*get<1 ,1>(M)+8*get<1 ,2>(M)*(-get<2 ,0>(M)*get<3 ,1>(M)+get<2 ,1>(M)*get<3 ,0>(M)))*get<0 ,3>(M)+(-2*pow2(get<1 ,0>(M))+2*pow2(get<1 ,1>(M))+2*pow2(get<1 ,2>(M))-2*pow2(get<1 ,3>(M))-2*pow2(get<2 ,0>(M))+2*pow2(get<2 ,1>(M))+2*pow2(get<2 ,2>(M))-2*pow2(get<2 ,3>(M))-2*pow2(get<3 ,0>(M))+2*pow2(get<3 ,1>(M))+2*pow2(get<3 ,2>(M))-2*pow2(get<3 ,3>(M)))*pow2(get<0 ,3>(M))+((-8*get<1 ,2>(M)*get<1 ,3>(M)-8*get<2 ,2>(M)*get<2 ,3>(M)-8*get<3 ,2>(M)*get<3 ,3>(M))*get<0 ,3>(M)+(8*get<2 ,1>(M)*get<3 ,3>(M)-8*get<2 ,3>(M)*get<3 ,1>(M))*get<1 ,0>(M)+(-8*get<2 ,0>(M)*get<3 ,3>(M)+8*get<2 ,3>(M)*get<3 ,0>(M))*get<1 ,1>(M)-8*get<1 ,3>(M)*(-get<2 ,0>(M)*get<3 ,1>(M)+get<2 ,1>(M)*get<3 ,0>(M)))*get<0 ,2>(M)+(2*pow2(get<0 ,3>(M))-2*pow2(get<1 ,0>(M))+2*pow2(get<1 ,1>(M))-2*pow2(get<1 ,2>(M))+2*pow2(get<1 ,3>(M))-2*pow2(get<2 ,0>(M))+2*pow2(get<2 ,1>(M))-2*pow2(get<2 ,2>(M))+2*pow2(get<2 ,3>(M))-2*pow2(get<3 ,0>(M))+2*pow2(get<3 ,1>(M))-2*pow2(get<3 ,2>(M))+2*pow2(get<3 ,3>(M)))*pow2(get<0 ,2>(M))+((-8*get<1 ,1>(M)*get<1 ,2>(M)-8*get<2 ,1>(M)*get<2 ,2>(M)-8*get<3 ,1>(M)*get<3 ,2>(M))*get<0 ,2>(M)+(-8*get<1 ,1>(M)*get<1 ,3>(M)-8*get<2 ,1>(M)*get<2 ,3>(M)-8*get<3 ,1>(M)*get<3 ,3>(M))*get<0 ,3>(M)+(-8*get<2 ,2>(M)*get<3 ,3>(M)+8*get<2 ,3>(M)*get<3 ,2>(M))*get<1 ,0>(M)+(8*get<2 ,0>(M)*get<3 ,3>(M)-8*get<2 ,3>(M)*get<3 ,0>(M))*get<1 ,2>(M)+8*get<1 ,3>(M)*(-get<2 ,0>(M)*get<3 ,2>(M)+get<2 ,2>(M)*get<3 ,0>(M)))*get<0 ,1>(M)+(2*pow2(get<0 ,2>(M))+2*pow2(get<0 ,3>(M))-2*pow2(get<1 ,0>(M))-2*pow2(get<1 ,1>(M))+2*pow2(get<1 ,2>(M))+2*pow2(get<1 ,3>(M))-2*pow2(get<2 ,0>(M))-2*pow2(get<2 ,1>(M))+2*pow2(get<2 ,2>(M))+2*pow2(get<2 ,3>(M))-2*pow2(get<3 ,0>(M))-2*pow2(get<3 ,1>(M))+2*pow2(get<3 ,2>(M))+2*pow2(get<3 ,3>(M)))*pow2(get<0 ,1>(M))+((8*get<1 ,0>(M)*get<1 ,1>(M)+8*get<2 ,0>(M)*get<2 ,1>(M)+8*get<3 ,0>(M)*get<3 ,1>(M))*get<0 ,1>(M)+(8*get<1 ,0>(M)*get<1 ,2>(M)+8*get<2 ,0>(M)*get<2 ,2>(M)+8*get<3 ,0>(M)*get<3 ,2>(M))*get<0 ,2>(M)+(8*get<1 ,0>(M)*get<1 ,3>(M)+8*get<2 ,0>(M)*get<2 ,3>(M)+8*get<3 ,0>(M)*get<3 ,3>(M))*get<0 ,3>(M)+(8*get<2 ,2>(M)*get<3 ,3>(M)-8*get<2 ,3>(M)*get<3 ,2>(M))*get<1 ,1>(M)+(-8*get<2 ,1>(M)*get<3 ,3>(M)+8*get<2 ,3>(M)*get<3 ,1>(M))*get<1 ,2>(M)+8*get<1 ,3>(M)*(get<2 ,1>(M)*get<3 ,2>(M)-get<2 ,2>(M)*get<3 ,1>(M)))*get<0 ,0>(M)+(-2*pow2(get<0 ,1>(M))-2*pow2(get<0 ,2>(M))-2*pow2(get<0 ,3>(M))-2*pow2(get<1 ,0>(M))-2*pow2(get<1 ,1>(M))-2*pow2(get<1 ,2>(M))-2*pow2(get<1 ,3>(M))-2*pow2(get<2 ,0>(M))-2*pow2(get<2 ,1>(M))-2*pow2(get<2 ,2>(M))-2*pow2(get<2 ,3>(M))-2*pow2(get<3 ,0>(M))-2*pow2(get<3 ,1>(M))-2*pow2(get<3 ,2>(M))-2*pow2(get<3 ,3>(M)))*pow2(get<0 ,0>(M));

    return ((C1 >= 0.) && (C2 >= 0.) && (C3 >= 0.) && (C4 >= 0.));
}

// --------------------------------
// Generic runners for PSD checkers
// --------------------------------

void run_checker_seq(bool (*checker)(double*), double* matrices, size_t num_matrices, bool* result) {
    size_t step = DIM * DIM;
    for (size_t i = 0; i < num_matrices; ++i) {
        result[i] = (*checker)(matrices + step * i);
    }
}

// --------------
// Interface code
// --------------

int main() {
    auto start = std::chrono::system_clock::now();
    double* matrices = read_matrices();
    auto end = std::chrono::system_clock::now();

    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    std::cout << "Matrices read from files in " << elapsed << " milliseconds" << std::endl;

    bool* result = (bool*)malloc(sizeof(bool) * WIDTH * HEIGHT);
    start = std::chrono::system_clock::now();
    run_checker_seq(&check_with_charpoly_noh, matrices, HEIGHT * WIDTH, result);
    end = std::chrono::system_clock::now();
    elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    std::cout << "Checking performed in " << elapsed << " milliseconds" << std::endl;

    size_t nonphysical_count = 0;
    for (size_t i = 0; i < HEIGHT * WIDTH; ++i) {
        if (!result[i]) {
            nonphysical_count++;
        }
    }
    std::cout << "The number of non-physical matrices is " << nonphysical_count << std::endl;

    free(matrices);
    free(result);
}
