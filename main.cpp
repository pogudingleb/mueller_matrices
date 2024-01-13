#include "main.h"

// ------------
// PSD Checkers
// ------------

typedef bool PSDChecker(const double*);

// Some common low-level functions
template <size_t i, size_t j>
inline double get(const double* M) {
    return M[i * DIM + j];
}

inline double pow2(const double x) {
    return x * x;
}

inline double pow3(const double x) {
    return x * x * x;
}

inline double pow4(const double x) {
    double y = x * x;
    return y * y;
}

void build_eigen_matrix(const double* M, Eigen::Matrix4cd& H) {
    H(0, 0).real(get<0, 0>(M) + get<0, 1>(M) + get<1, 0>(M) + get<1, 1>(M));
    H(0, 0).imag(0);

    H(0, 1).real(get<0, 2>(M) + get<1, 2>(M));
    H(0, 1).imag(get<0, 3>(M) + get<1, 3>(M));

    H(0, 2).real(get<2, 0>(M) + get<2, 1>(M));
    H(0, 2).imag(-get<3, 0>(M) - get<3, 1>(M));

    H(0, 3).real(get<2, 2>(M) + get<3, 3>(M));
    H(0, 3).imag(get<2, 3>(M) - get<3, 2>(M));

    H(1, 0).real(H(0, 1).real());
    H(1, 0).imag(-H(0, 1).imag());

    H(1, 1).real(get<0, 0>(M) - get<0, 1>(M) + get<1, 0>(M) - get<1, 1>(M));
    H(1, 1).imag(0);

    H(1, 2).real(get<2, 2>(M) - get<3, 3>(M));
    H(1, 2).imag(-get<2, 3>(M) - get<3, 2>(M));

    H(1, 3).real(get<2, 0>(M) - get<2, 1>(M));
    H(1, 3).imag(-get<3, 0>(M) + get<3, 1>(M));

    H(2, 0).real(H(0, 2).real());
    H(2, 0).imag(-H(0, 2).imag());

    H(2, 1).real(H(1, 2).real());
    H(2, 1).imag(-H(1, 2).imag());

    H(2, 2).real(get<0, 0>(M) + get<0, 1>(M) - get<1, 0>(M) - get<1, 1>(M));
    H(2, 2).imag(0);

    H(2, 3).real(get<0, 2>(M) - get<1, 2>(M));
    H(2, 3).imag(get<0, 3>(M) - get<1, 3>(M));

    H(3, 0).real(H(0, 3).real());
    H(3, 0).imag(-H(0, 3).imag());

    H(3, 1).real(H(1, 3).real());
    H(3, 1).imag(-H(1, 3).imag());

    H(3, 2).real(H(2, 3).real());
    H(3, 2).imag(-H(2, 3).imag());

    H(3, 3).real(get<0, 0>(M) - get<0, 1>(M) - get<1, 0>(M) + get<1, 1>(M));
    H(3, 3).imag(0);
}

// Checker via the charpoly (not computing H directly)

bool check_with_charpoly_noh(const double* M) {
    double C1 = get<0, 0>(M);

    double C2 = 3 * pow2(get<0 ,0>(M)) - pow2(get<0 ,1>(M)) - pow2(get<0 ,2>(M)) - pow2(get<0 ,3>(M)) - pow2(get<1 ,0>(M)) - pow2(get<1 ,1>(M)) - pow2(get<1 ,2>(M)) - pow2(get<1 ,3>(M)) - pow2(get<2 ,0>(M)) - pow2(get<2 ,1>(M)) - pow2(get<2 ,2>(M)) - pow2(get<2 ,3>(M)) - pow2(get<3 ,0>(M)) - pow2(get<3 ,1>(M)) - pow2(get<3 ,2>(M)) - pow2(get<3 ,3>(M));

    double C3 = 4*pow3(get<0 ,0>(M))+(-4*pow2(get<0 ,1>(M))-4*pow2(get<0 ,2>(M))-4*pow2(get<0 ,3>(M))-4*pow2(get<1 ,0>(M))-4*pow2(get<1 ,1>(M))-4*pow2(get<1 ,2>(M))-4*pow2(get<1 ,3>(M))-4*pow2(get<2 ,0>(M))-4*pow2(get<2 ,1>(M))-4*pow2(get<2 ,2>(M))-4*pow2(get<2 ,3>(M))-4*pow2(get<3 ,0>(M))-4*pow2(get<3 ,1>(M))-4*pow2(get<3 ,2>(M))-4*pow2(get<3 ,3>(M)))*get<0 ,0>(M)+(8*get<1 ,0>(M)*get<1 ,1>(M)+8*get<2 ,0>(M)*get<2 ,1>(M)+8*get<3 ,0>(M)*get<3 ,1>(M))*get<0 ,1>(M)+(8*get<1 ,0>(M)*get<1 ,2>(M)+8*get<2 ,0>(M)*get<2 ,2>(M)+8*get<3 ,0>(M)*get<3 ,2>(M))*get<0 ,2>(M)+(8*get<1 ,0>(M)*get<1 ,3>(M)+8*get<2 ,0>(M)*get<2 ,3>(M)+8*get<3 ,0>(M)*get<3 ,3>(M))*get<0 ,3>(M)+(8*get<2 ,2>(M)*get<3 ,3>(M)-8*get<2 ,3>(M)*get<3 ,2>(M))*get<1 ,1>(M)+(-8*get<2 ,1>(M)*get<3 ,3>(M)+8*get<2 ,3>(M)*get<3 ,1>(M))*get<1 ,2>(M)+8*get<1 ,3>(M)*(get<2 ,1>(M)*get<3 ,2>(M)-get<2 ,2>(M)*get<3 ,1>(M));

    double C4 = pow4(get<0 ,0>(M))+pow4(get<0 ,1>(M))+pow4(get<0 ,2>(M))+pow4(get<1 ,0>(M))+pow4(get<1 ,1>(M))+pow4(get<1 ,2>(M))+pow4(get<2 ,0>(M))+pow4(get<2 ,1>(M))+pow4(get<2 ,2>(M))+pow4(get<0 ,3>(M))+pow4(get<1 ,3>(M))+pow4(get<2 ,3>(M))+8*get<2 ,2>(M)*get<2 ,3>(M)*get<3 ,2>(M)*get<3 ,3>(M)-8*get<3 ,0>(M)*(get<2 ,1>(M)*get<3 ,1>(M)+get<2 ,2>(M)*get<3 ,2>(M)+get<2 ,3>(M)*get<3 ,3>(M))*get<2 ,0>(M)+8*get<3 ,1>(M)*(get<2 ,2>(M)*get<3 ,2>(M)+get<2 ,3>(M)*get<3 ,3>(M))*get<2 ,1>(M)+8*get<1 ,3>(M)*(get<2 ,2>(M)*get<2 ,3>(M)+get<3 ,2>(M)*get<3 ,3>(M))*get<1 ,2>(M)+pow2(pow2(get<3 ,0>(M))-pow2(get<3 ,1>(M))-pow2(get<3 ,2>(M))-pow2(get<3 ,3>(M)))+(2*pow2(get<3 ,0>(M))-2*pow2(get<3 ,1>(M))-2*pow2(get<3 ,2>(M))+2*pow2(get<3 ,3>(M)))*pow2(get<2 ,3>(M))+(2*pow2(get<2 ,3>(M))+2*pow2(get<3 ,0>(M))-2*pow2(get<3 ,1>(M))+2*pow2(get<3 ,2>(M))-2*pow2(get<3 ,3>(M)))*pow2(get<2 ,2>(M))+(2*pow2(get<2 ,2>(M))+2*pow2(get<2 ,3>(M))+2*pow2(get<3 ,0>(M))+2*pow2(get<3 ,1>(M))-2*pow2(get<3 ,2>(M))-2*pow2(get<3 ,3>(M)))*pow2(get<2 ,1>(M))+(-2*pow2(get<2 ,1>(M))-2*pow2(get<2 ,2>(M))-2*pow2(get<2 ,3>(M))+2*pow2(get<3 ,0>(M))+2*pow2(get<3 ,1>(M))+2*pow2(get<3 ,2>(M))+2*pow2(get<3 ,3>(M)))*pow2(get<2 ,0>(M))+(2*pow2(get<2 ,0>(M))-2*pow2(get<2 ,1>(M))-2*pow2(get<2 ,2>(M))+2*pow2(get<2 ,3>(M))+2*pow2(get<3 ,0>(M))-2*pow2(get<3 ,1>(M))-2*pow2(get<3 ,2>(M))+2*pow2(get<3 ,3>(M)))*pow2(get<1 ,3>(M))+(2*pow2(get<1 ,3>(M))+2*pow2(get<2 ,0>(M))-2*pow2(get<2 ,1>(M))+2*pow2(get<2 ,2>(M))-2*pow2(get<2 ,3>(M))+2*pow2(get<3 ,0>(M))-2*pow2(get<3 ,1>(M))+2*pow2(get<3 ,2>(M))-2*pow2(get<3 ,3>(M)))*pow2(get<1 ,2>(M))+((8*get<2 ,1>(M)*get<2 ,2>(M)+8*get<3 ,1>(M)*get<3 ,2>(M))*get<1 ,2>(M)+8*get<1 ,3>(M)*(get<2 ,1>(M)*get<2 ,3>(M)+get<3 ,1>(M)*get<3 ,3>(M)))*get<1 ,1>(M)+(2*pow2(get<1 ,2>(M))+2*pow2(get<1 ,3>(M))+2*pow2(get<2 ,0>(M))+2*pow2(get<2 ,1>(M))-2*pow2(get<2 ,2>(M))-2*pow2(get<2 ,3>(M))+2*pow2(get<3 ,0>(M))+2*pow2(get<3 ,1>(M))-2*pow2(get<3 ,2>(M))-2*pow2(get<3 ,3>(M)))*pow2(get<1 ,1>(M))+((-8*get<2 ,0>(M)*get<2 ,1>(M)-8*get<3 ,0>(M)*get<3 ,1>(M))*get<1 ,1>(M)+(-8*get<2 ,0>(M)*get<2 ,2>(M)-8*get<3 ,0>(M)*get<3 ,2>(M))*get<1 ,2>(M)-8*get<1 ,3>(M)*(get<2 ,0>(M)*get<2 ,3>(M)+get<3 ,0>(M)*get<3 ,3>(M)))*get<1 ,0>(M)+(-2*pow2(get<1 ,1>(M))-2*pow2(get<1 ,2>(M))-2*pow2(get<1 ,3>(M))+2*pow2(get<2 ,0>(M))+2*pow2(get<2 ,1>(M))+2*pow2(get<2 ,2>(M))+2*pow2(get<2 ,3>(M))+2*pow2(get<3 ,0>(M))+2*pow2(get<3 ,1>(M))+2*pow2(get<3 ,2>(M))+2*pow2(get<3 ,3>(M)))*pow2(get<1 ,0>(M))+((-8*get<2 ,1>(M)*get<3 ,2>(M)+8*get<2 ,2>(M)*get<3 ,1>(M))*get<1 ,0>(M)+(8*get<2 ,0>(M)*get<3 ,2>(M)-8*get<2 ,2>(M)*get<3 ,0>(M))*get<1 ,1>(M)+8*get<1 ,2>(M)*(-get<2 ,0>(M)*get<3 ,1>(M)+get<2 ,1>(M)*get<3 ,0>(M)))*get<0 ,3>(M)+(-2*pow2(get<1 ,0>(M))+2*pow2(get<1 ,1>(M))+2*pow2(get<1 ,2>(M))-2*pow2(get<1 ,3>(M))-2*pow2(get<2 ,0>(M))+2*pow2(get<2 ,1>(M))+2*pow2(get<2 ,2>(M))-2*pow2(get<2 ,3>(M))-2*pow2(get<3 ,0>(M))+2*pow2(get<3 ,1>(M))+2*pow2(get<3 ,2>(M))-2*pow2(get<3 ,3>(M)))*pow2(get<0 ,3>(M))+((-8*get<1 ,2>(M)*get<1 ,3>(M)-8*get<2 ,2>(M)*get<2 ,3>(M)-8*get<3 ,2>(M)*get<3 ,3>(M))*get<0 ,3>(M)+(8*get<2 ,1>(M)*get<3 ,3>(M)-8*get<2 ,3>(M)*get<3 ,1>(M))*get<1 ,0>(M)+(-8*get<2 ,0>(M)*get<3 ,3>(M)+8*get<2 ,3>(M)*get<3 ,0>(M))*get<1 ,1>(M)-8*get<1 ,3>(M)*(-get<2 ,0>(M)*get<3 ,1>(M)+get<2 ,1>(M)*get<3 ,0>(M)))*get<0 ,2>(M)+(2*pow2(get<0 ,3>(M))-2*pow2(get<1 ,0>(M))+2*pow2(get<1 ,1>(M))-2*pow2(get<1 ,2>(M))+2*pow2(get<1 ,3>(M))-2*pow2(get<2 ,0>(M))+2*pow2(get<2 ,1>(M))-2*pow2(get<2 ,2>(M))+2*pow2(get<2 ,3>(M))-2*pow2(get<3 ,0>(M))+2*pow2(get<3 ,1>(M))-2*pow2(get<3 ,2>(M))+2*pow2(get<3 ,3>(M)))*pow2(get<0 ,2>(M))+((-8*get<1 ,1>(M)*get<1 ,2>(M)-8*get<2 ,1>(M)*get<2 ,2>(M)-8*get<3 ,1>(M)*get<3 ,2>(M))*get<0 ,2>(M)+(-8*get<1 ,1>(M)*get<1 ,3>(M)-8*get<2 ,1>(M)*get<2 ,3>(M)-8*get<3 ,1>(M)*get<3 ,3>(M))*get<0 ,3>(M)+(-8*get<2 ,2>(M)*get<3 ,3>(M)+8*get<2 ,3>(M)*get<3 ,2>(M))*get<1 ,0>(M)+(8*get<2 ,0>(M)*get<3 ,3>(M)-8*get<2 ,3>(M)*get<3 ,0>(M))*get<1 ,2>(M)+8*get<1 ,3>(M)*(-get<2 ,0>(M)*get<3 ,2>(M)+get<2 ,2>(M)*get<3 ,0>(M)))*get<0 ,1>(M)+(2*pow2(get<0 ,2>(M))+2*pow2(get<0 ,3>(M))-2*pow2(get<1 ,0>(M))-2*pow2(get<1 ,1>(M))+2*pow2(get<1 ,2>(M))+2*pow2(get<1 ,3>(M))-2*pow2(get<2 ,0>(M))-2*pow2(get<2 ,1>(M))+2*pow2(get<2 ,2>(M))+2*pow2(get<2 ,3>(M))-2*pow2(get<3 ,0>(M))-2*pow2(get<3 ,1>(M))+2*pow2(get<3 ,2>(M))+2*pow2(get<3 ,3>(M)))*pow2(get<0 ,1>(M))+((8*get<1 ,0>(M)*get<1 ,1>(M)+8*get<2 ,0>(M)*get<2 ,1>(M)+8*get<3 ,0>(M)*get<3 ,1>(M))*get<0 ,1>(M)+(8*get<1 ,0>(M)*get<1 ,2>(M)+8*get<2 ,0>(M)*get<2 ,2>(M)+8*get<3 ,0>(M)*get<3 ,2>(M))*get<0 ,2>(M)+(8*get<1 ,0>(M)*get<1 ,3>(M)+8*get<2 ,0>(M)*get<2 ,3>(M)+8*get<3 ,0>(M)*get<3 ,3>(M))*get<0 ,3>(M)+(8*get<2 ,2>(M)*get<3 ,3>(M)-8*get<2 ,3>(M)*get<3 ,2>(M))*get<1 ,1>(M)+(-8*get<2 ,1>(M)*get<3 ,3>(M)+8*get<2 ,3>(M)*get<3 ,1>(M))*get<1 ,2>(M)+8*get<1 ,3>(M)*(get<2 ,1>(M)*get<3 ,2>(M)-get<2 ,2>(M)*get<3 ,1>(M)))*get<0 ,0>(M)+(-2*pow2(get<0 ,1>(M))-2*pow2(get<0 ,2>(M))-2*pow2(get<0 ,3>(M))-2*pow2(get<1 ,0>(M))-2*pow2(get<1 ,1>(M))-2*pow2(get<1 ,2>(M))-2*pow2(get<1 ,3>(M))-2*pow2(get<2 ,0>(M))-2*pow2(get<2 ,1>(M))-2*pow2(get<2 ,2>(M))-2*pow2(get<2 ,3>(M))-2*pow2(get<3 ,0>(M))-2*pow2(get<3 ,1>(M))-2*pow2(get<3 ,2>(M))-2*pow2(get<3 ,3>(M)))*pow2(get<0 ,0>(M));

    return ((C1 >= 0.) && (C2 >= 0.) && (C3 >= 0.) && (C4 >= 0.));
}

// Checker via the Sylvester criterion (not computing H directly)

bool check_with_sylvester_noh(const double* M) {
    double C1 = get<0, 0>(M) + get<1, 1>(M) + get<2, 2>(M) + get<3, 3>(M);

    double C2 = pow2(get<0 ,0>(M))+2*get<0 ,0>(M)*get<1 ,0>(M)-pow2(get<0 ,1>(M))-2*get<0 ,1>(M)*get<1 ,1>(M)-pow2(get<0 ,2>(M))-2*get<0 ,2>(M)*get<1 ,2>(M)-pow2(get<0 ,3>(M))-2*get<0 ,3>(M)*get<1 ,3>(M)+pow2(get<1 ,0>(M))-pow2(get<1 ,1>(M))-pow2(get<1 ,2>(M))-pow2(get<1 ,3>(M));

    double C3 = pow3(get<0 ,0>(M))+(get<0 ,1>(M)+get<1 ,0>(M)-get<1 ,1>(M))*pow2(get<0 ,0>(M))+(-pow2(get<0 ,1>(M))+(2*get<1 ,0>(M)-2*get<1 ,1>(M))*get<0 ,1>(M)-pow2(get<2 ,1>(M))-2*get<2 ,0>(M)*get<2 ,1>(M)-pow2(get<2 ,2>(M))+2*get<2 ,2>(M)*get<3 ,3>(M)-pow2(get<2 ,3>(M))-2*get<2 ,3>(M)*get<3 ,2>(M)-pow2(get<3 ,0>(M))-2*get<3 ,0>(M)*get<3 ,1>(M)-pow2(get<3 ,1>(M))-pow2(get<3 ,2>(M))-pow2(get<3 ,3>(M))-pow2(get<0 ,2>(M))-2*get<0 ,2>(M)*get<1 ,2>(M)-pow2(get<0 ,3>(M))-2*get<0 ,3>(M)*get<1 ,3>(M)-pow2(get<1 ,0>(M))-2*get<1 ,0>(M)*get<1 ,1>(M)-pow2(get<1 ,1>(M))-pow2(get<1 ,2>(M))-pow2(get<1 ,3>(M))-pow2(get<2 ,0>(M)))*get<0 ,0>(M)-pow3(get<0 ,1>(M))+(get<1 ,0>(M)-get<1 ,1>(M))*pow2(get<0 ,1>(M))+(-pow2(get<0 ,2>(M))-2*get<0 ,2>(M)*get<1 ,2>(M)-pow2(get<0 ,3>(M))-2*get<0 ,3>(M)*get<1 ,3>(M)+pow2(get<1 ,0>(M))+2*get<1 ,0>(M)*get<1 ,1>(M)+pow2(get<1 ,1>(M))-pow2(get<1 ,2>(M))-pow2(get<1 ,3>(M))+pow2(get<2 ,0>(M))+2*get<2 ,0>(M)*get<2 ,1>(M)+pow2(get<2 ,1>(M))-pow2(get<2 ,2>(M))+2*get<2 ,2>(M)*get<3 ,3>(M)-pow2(get<2 ,3>(M))-2*get<2 ,3>(M)*get<3 ,2>(M)+pow2(get<3 ,0>(M))+2*get<3 ,0>(M)*get<3 ,1>(M)+pow2(get<3 ,1>(M))-pow2(get<3 ,2>(M))-pow2(get<3 ,3>(M)))*get<0 ,1>(M)-pow3(get<1 ,0>(M))-pow2(get<1 ,0>(M))*get<1 ,1>(M)+(pow2(get<0 ,2>(M))+2*get<0 ,2>(M)*get<1 ,2>(M)+pow2(get<0 ,3>(M))+2*get<0 ,3>(M)*get<1 ,3>(M)+pow2(get<1 ,1>(M))+pow2(get<1 ,2>(M))+pow2(get<1 ,3>(M))-pow2(get<2 ,0>(M))-2*get<2 ,0>(M)*get<2 ,1>(M)-pow2(get<2 ,1>(M))-pow2(get<2 ,2>(M))+2*get<2 ,2>(M)*get<3 ,3>(M)-pow2(get<2 ,3>(M))-2*get<2 ,3>(M)*get<3 ,2>(M)-pow2(get<3 ,0>(M))-2*get<3 ,0>(M)*get<3 ,1>(M)-pow2(get<3 ,1>(M))-pow2(get<3 ,2>(M))-pow2(get<3 ,3>(M)))*get<1 ,0>(M)+pow3(get<1 ,1>(M))+(pow2(get<0 ,2>(M))+2*get<0 ,2>(M)*get<1 ,2>(M)+pow2(get<0 ,3>(M))+2*get<0 ,3>(M)*get<1 ,3>(M)+pow2(get<1 ,2>(M))+pow2(get<1 ,3>(M))+pow2(get<2 ,0>(M))+2*get<2 ,0>(M)*get<2 ,1>(M)+pow2(get<2 ,1>(M))-pow2(get<2 ,2>(M))+2*get<2 ,2>(M)*get<3 ,3>(M)-pow2(get<2 ,3>(M))-2*get<2 ,3>(M)*get<3 ,2>(M)+pow2(get<3 ,0>(M))+2*get<3 ,0>(M)*get<3 ,1>(M)+pow2(get<3 ,1>(M))-pow2(get<3 ,2>(M))-pow2(get<3 ,3>(M)))*get<1 ,1>(M)+((2*get<2 ,2>(M)-2*get<3 ,3>(M))*get<2 ,0>(M)+(2*get<2 ,2>(M)-2*get<3 ,3>(M))*get<2 ,1>(M)+2*(get<3 ,0>(M)+get<3 ,1>(M))*(get<2 ,3>(M)+get<3 ,2>(M)))*get<0 ,2>(M)+((2*get<2 ,3>(M)+2*get<3 ,2>(M))*get<2 ,0>(M)+(2*get<2 ,3>(M)+2*get<3 ,2>(M))*get<2 ,1>(M)-2*(get<3 ,0>(M)+get<3 ,1>(M))*(get<2 ,2>(M)-get<3 ,3>(M)))*get<0 ,3>(M)+((2*get<2 ,2>(M)-2*get<3 ,3>(M))*get<2 ,0>(M)+(2*get<2 ,2>(M)-2*get<3 ,3>(M))*get<2 ,1>(M)+2*(get<3 ,0>(M)+get<3 ,1>(M))*(get<2 ,3>(M)+get<3 ,2>(M)))*get<1 ,2>(M)+2*get<1 ,3>(M)*((get<2 ,3>(M)+get<3 ,2>(M))*get<2 ,0>(M)+(get<2 ,3>(M)+get<3 ,2>(M))*get<2 ,1>(M)-(get<3 ,0>(M)+get<3 ,1>(M))*(get<2 ,2>(M)-get<3 ,3>(M)));

    double C4 = pow4(get<0 ,0>(M))+pow4(get<0 ,1>(M))+pow4(get<0 ,2>(M))+pow4(get<1 ,0>(M))+pow4(get<1 ,1>(M))+pow4(get<1 ,2>(M))+pow4(get<2 ,0>(M))+pow4(get<2 ,1>(M))+pow4(get<2 ,2>(M))+pow4(get<0 ,3>(M))+pow4(get<1 ,3>(M))+pow4(get<2 ,3>(M))+8*get<2 ,2>(M)*get<2 ,3>(M)*get<3 ,2>(M)*get<3 ,3>(M)-8*get<3 ,0>(M)*(get<2 ,1>(M)*get<3 ,1>(M)+get<2 ,2>(M)*get<3 ,2>(M)+get<2 ,3>(M)*get<3 ,3>(M))*get<2 ,0>(M)+8*get<3 ,1>(M)*(get<2 ,2>(M)*get<3 ,2>(M)+get<2 ,3>(M)*get<3 ,3>(M))*get<2 ,1>(M)+8*get<1 ,3>(M)*(get<2 ,2>(M)*get<2 ,3>(M)+get<3 ,2>(M)*get<3 ,3>(M))*get<1 ,2>(M)+pow2(pow2(get<3 ,0>(M))-pow2(get<3 ,1>(M))-pow2(get<3 ,2>(M))-pow2(get<3 ,3>(M)))+(2*pow2(get<3 ,0>(M))-2*pow2(get<3 ,1>(M))-2*pow2(get<3 ,2>(M))+2*pow2(get<3 ,3>(M)))*pow2(get<2 ,3>(M))+(2*pow2(get<2 ,3>(M))+2*pow2(get<3 ,0>(M))-2*pow2(get<3 ,1>(M))+2*pow2(get<3 ,2>(M))-2*pow2(get<3 ,3>(M)))*pow2(get<2 ,2>(M))+(2*pow2(get<2 ,2>(M))+2*pow2(get<2 ,3>(M))+2*pow2(get<3 ,0>(M))+2*pow2(get<3 ,1>(M))-2*pow2(get<3 ,2>(M))-2*pow2(get<3 ,3>(M)))*pow2(get<2 ,1>(M))+(-2*pow2(get<2 ,1>(M))-2*pow2(get<2 ,2>(M))-2*pow2(get<2 ,3>(M))+2*pow2(get<3 ,0>(M))+2*pow2(get<3 ,1>(M))+2*pow2(get<3 ,2>(M))+2*pow2(get<3 ,3>(M)))*pow2(get<2 ,0>(M))+(2*pow2(get<2 ,0>(M))-2*pow2(get<2 ,1>(M))-2*pow2(get<2 ,2>(M))+2*pow2(get<2 ,3>(M))+2*pow2(get<3 ,0>(M))-2*pow2(get<3 ,1>(M))-2*pow2(get<3 ,2>(M))+2*pow2(get<3 ,3>(M)))*pow2(get<1 ,3>(M))+(2*pow2(get<1 ,3>(M))+2*pow2(get<2 ,0>(M))-2*pow2(get<2 ,1>(M))+2*pow2(get<2 ,2>(M))-2*pow2(get<2 ,3>(M))+2*pow2(get<3 ,0>(M))-2*pow2(get<3 ,1>(M))+2*pow2(get<3 ,2>(M))-2*pow2(get<3 ,3>(M)))*pow2(get<1 ,2>(M))+((8*get<2 ,1>(M)*get<2 ,2>(M)+8*get<3 ,1>(M)*get<3 ,2>(M))*get<1 ,2>(M)+8*get<1 ,3>(M)*(get<2 ,1>(M)*get<2 ,3>(M)+get<3 ,1>(M)*get<3 ,3>(M)))*get<1 ,1>(M)+(2*pow2(get<1 ,2>(M))+2*pow2(get<1 ,3>(M))+2*pow2(get<2 ,0>(M))+2*pow2(get<2 ,1>(M))-2*pow2(get<2 ,2>(M))-2*pow2(get<2 ,3>(M))+2*pow2(get<3 ,0>(M))+2*pow2(get<3 ,1>(M))-2*pow2(get<3 ,2>(M))-2*pow2(get<3 ,3>(M)))*pow2(get<1 ,1>(M))+((-8*get<2 ,0>(M)*get<2 ,1>(M)-8*get<3 ,0>(M)*get<3 ,1>(M))*get<1 ,1>(M)+(-8*get<2 ,0>(M)*get<2 ,2>(M)-8*get<3 ,0>(M)*get<3 ,2>(M))*get<1 ,2>(M)-8*get<1 ,3>(M)*(get<2 ,0>(M)*get<2 ,3>(M)+get<3 ,0>(M)*get<3 ,3>(M)))*get<1 ,0>(M)+(-2*pow2(get<1 ,1>(M))-2*pow2(get<1 ,2>(M))-2*pow2(get<1 ,3>(M))+2*pow2(get<2 ,0>(M))+2*pow2(get<2 ,1>(M))+2*pow2(get<2 ,2>(M))+2*pow2(get<2 ,3>(M))+2*pow2(get<3 ,0>(M))+2*pow2(get<3 ,1>(M))+2*pow2(get<3 ,2>(M))+2*pow2(get<3 ,3>(M)))*pow2(get<1 ,0>(M))+((-8*get<2 ,1>(M)*get<3 ,2>(M)+8*get<2 ,2>(M)*get<3 ,1>(M))*get<1 ,0>(M)+(8*get<2 ,0>(M)*get<3 ,2>(M)-8*get<2 ,2>(M)*get<3 ,0>(M))*get<1 ,1>(M)+8*get<1 ,2>(M)*(-get<2 ,0>(M)*get<3 ,1>(M)+get<2 ,1>(M)*get<3 ,0>(M)))*get<0 ,3>(M)+(-2*pow2(get<1 ,0>(M))+2*pow2(get<1 ,1>(M))+2*pow2(get<1 ,2>(M))-2*pow2(get<1 ,3>(M))-2*pow2(get<2 ,0>(M))+2*pow2(get<2 ,1>(M))+2*pow2(get<2 ,2>(M))-2*pow2(get<2 ,3>(M))-2*pow2(get<3 ,0>(M))+2*pow2(get<3 ,1>(M))+2*pow2(get<3 ,2>(M))-2*pow2(get<3 ,3>(M)))*pow2(get<0 ,3>(M))+((-8*get<1 ,2>(M)*get<1 ,3>(M)-8*get<2 ,2>(M)*get<2 ,3>(M)-8*get<3 ,2>(M)*get<3 ,3>(M))*get<0 ,3>(M)+(8*get<2 ,1>(M)*get<3 ,3>(M)-8*get<2 ,3>(M)*get<3 ,1>(M))*get<1 ,0>(M)+(-8*get<2 ,0>(M)*get<3 ,3>(M)+8*get<2 ,3>(M)*get<3 ,0>(M))*get<1 ,1>(M)-8*get<1 ,3>(M)*(-get<2 ,0>(M)*get<3 ,1>(M)+get<2 ,1>(M)*get<3 ,0>(M)))*get<0 ,2>(M)+(2*pow2(get<0 ,3>(M))-2*pow2(get<1 ,0>(M))+2*pow2(get<1 ,1>(M))-2*pow2(get<1 ,2>(M))+2*pow2(get<1 ,3>(M))-2*pow2(get<2 ,0>(M))+2*pow2(get<2 ,1>(M))-2*pow2(get<2 ,2>(M))+2*pow2(get<2 ,3>(M))-2*pow2(get<3 ,0>(M))+2*pow2(get<3 ,1>(M))-2*pow2(get<3 ,2>(M))+2*pow2(get<3 ,3>(M)))*pow2(get<0 ,2>(M))+((-8*get<1 ,1>(M)*get<1 ,2>(M)-8*get<2 ,1>(M)*get<2 ,2>(M)-8*get<3 ,1>(M)*get<3 ,2>(M))*get<0 ,2>(M)+(-8*get<1 ,1>(M)*get<1 ,3>(M)-8*get<2 ,1>(M)*get<2 ,3>(M)-8*get<3 ,1>(M)*get<3 ,3>(M))*get<0 ,3>(M)+(-8*get<2 ,2>(M)*get<3 ,3>(M)+8*get<2 ,3>(M)*get<3 ,2>(M))*get<1 ,0>(M)+(8*get<2 ,0>(M)*get<3 ,3>(M)-8*get<2 ,3>(M)*get<3 ,0>(M))*get<1 ,2>(M)+8*get<1 ,3>(M)*(-get<2 ,0>(M)*get<3 ,2>(M)+get<2 ,2>(M)*get<3 ,0>(M)))*get<0 ,1>(M)+(2*pow2(get<0 ,2>(M))+2*pow2(get<0 ,3>(M))-2*pow2(get<1 ,0>(M))-2*pow2(get<1 ,1>(M))+2*pow2(get<1 ,2>(M))+2*pow2(get<1 ,3>(M))-2*pow2(get<2 ,0>(M))-2*pow2(get<2 ,1>(M))+2*pow2(get<2 ,2>(M))+2*pow2(get<2 ,3>(M))-2*pow2(get<3 ,0>(M))-2*pow2(get<3 ,1>(M))+2*pow2(get<3 ,2>(M))+2*pow2(get<3 ,3>(M)))*pow2(get<0 ,1>(M))+((8*get<1 ,0>(M)*get<1 ,1>(M)+8*get<2 ,0>(M)*get<2 ,1>(M)+8*get<3 ,0>(M)*get<3 ,1>(M))*get<0 ,1>(M)+(8*get<1 ,0>(M)*get<1 ,2>(M)+8*get<2 ,0>(M)*get<2 ,2>(M)+8*get<3 ,0>(M)*get<3 ,2>(M))*get<0 ,2>(M)+(8*get<1 ,0>(M)*get<1 ,3>(M)+8*get<2 ,0>(M)*get<2 ,3>(M)+8*get<3 ,0>(M)*get<3 ,3>(M))*get<0 ,3>(M)+(8*get<2 ,2>(M)*get<3 ,3>(M)-8*get<2 ,3>(M)*get<3 ,2>(M))*get<1 ,1>(M)+(-8*get<2 ,1>(M)*get<3 ,3>(M)+8*get<2 ,3>(M)*get<3 ,1>(M))*get<1 ,2>(M)+8*get<1 ,3>(M)*(get<2 ,1>(M)*get<3 ,2>(M)-get<2 ,2>(M)*get<3 ,1>(M)))*get<0 ,0>(M)+(-2*pow2(get<0 ,1>(M))-2*pow2(get<0 ,2>(M))-2*pow2(get<0 ,3>(M))-2*pow2(get<1 ,0>(M))-2*pow2(get<1 ,1>(M))-2*pow2(get<1 ,2>(M))-2*pow2(get<1 ,3>(M))-2*pow2(get<2 ,0>(M))-2*pow2(get<2 ,1>(M))-2*pow2(get<2 ,2>(M))-2*pow2(get<2 ,3>(M))-2*pow2(get<3 ,0>(M))-2*pow2(get<3 ,1>(M))-2*pow2(get<3 ,2>(M))-2*pow2(get<3 ,3>(M)))*pow2(get<0 ,0>(M));

    return ((C1 > 0) && (C2 > 0) && (C3 > 0) && (C4 > 0));
}

// Checker using Eigen's self-adjoint solver
// https://eigen.tuxfamily.org/dox/classEigen_1_1SelfAdjointEigenSolver.html#aaf4ed4172a517a4b9f0ab222f629e261
bool check_eignevalues(const double* M) {
    Eigen::Matrix4cd H;
    build_eigen_matrix(M, H);
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix4cd> es;
    es.compute(H);
    // eigenvalues are sorted
    return (es.eigenvalues()(0, 0) > 0);
}

// Checker with Eigen's Choletsky decomposition
// https://eigen.tuxfamily.org/dox/classEigen_1_1LLT.html
bool check_choletsky(const double* M) {
    Eigen::Matrix4cd H;
    build_eigen_matrix(M, H);
    Eigen::LLT<Eigen::Matrix4cd> llt;
    llt.compute(H);
    return (llt.info() == Eigen::Success);
}

// Checker with the formulas based on the properties of Pauli matrices
bool check_pauli(const double* M) {
    double traceH = 2 * get<0, 0>(M);
    
    double squares[DIM * DIM];
    // should be a better way with template unrolling
    squares[0] = pow2(get<0, 0>(M));
    squares[1] = pow2(get<0, 1>(M));
    squares[2] = pow2(get<0, 2>(M));
    squares[3] = pow2(get<0, 3>(M));
    squares[4] = pow2(get<1, 0>(M));
    squares[5] = pow2(get<1, 1>(M));
    squares[6] = pow2(get<1, 2>(M));
    squares[7] = pow2(get<1, 3>(M));
    squares[8] = pow2(get<2, 0>(M));
    squares[9] = pow2(get<2, 1>(M));
    squares[10] = pow2(get<2, 2>(M));
    squares[11] = pow2(get<2, 3>(M));
    squares[12] = pow2(get<3, 0>(M));
    squares[13] = pow2(get<3, 1>(M));
    squares[14] = pow2(get<3, 2>(M));
    squares[15] = pow2(get<3, 3>(M));

    // different aggregations of squares
    // rowsums for lower-right 3-by-3 corner
    double row_squaresums[3];
    row_squaresums[0] = get<1, 1>(squares) + get<1, 2>(squares) + get<1, 3>(squares);
    row_squaresums[1] = get<2, 1>(squares) + get<2, 2>(squares) + get<2, 3>(squares);
    row_squaresums[2] = get<3, 1>(squares) + get<3, 2>(squares) + get<3, 3>(squares);

    double squares0x = get<0, 1>(squares) + get<0, 2>(squares) + get<0, 3>(squares);
    double squaresx0 = get<1, 0>(squares) + get<2, 0>(squares) + get<3, 0>(squares);
    double corner_squaresum = row_squaresums[0] + row_squaresums[1] + row_squaresums[2];

    double traceH2 = get<0, 0>(squares) + squares0x + squaresx0 + corner_squaresum;

    // computing the determinant of M
    double det01 = get<0, 2>(M) * get<1, 3>(M) - get<0, 3>(M) * get<1, 2>(M);
    double det02 = get<0, 2>(M) * get<2, 3>(M) - get<0, 3>(M) * get<2, 2>(M);
    double det03 = get<0, 2>(M) * get<3, 3>(M) - get<0, 3>(M) * get<3, 2>(M);
    double det12 = get<1, 2>(M) * get<2, 3>(M) - get<1, 3>(M) * get<2, 2>(M);
    double det13 = get<1, 2>(M) * get<3, 3>(M) - get<1, 3>(M) * get<3, 2>(M);
    double det23 = get<2, 2>(M) * get<3, 3>(M) - get<2, 3>(M) * get<3, 2>(M);

    double det012 = get<0, 1>(M) * det12 - get<1, 1>(M) * det02 + get<2, 1>(M) * det01;
    double det013 = get<0, 1>(M) * det13 - get<1, 1>(M) * det03 + get<3, 1>(M) * det01;
    double det023 = get<0, 1>(M) * det23 - get<2, 1>(M) * det03 + get<3, 1>(M) * det02;
    double det123 = get<1, 1>(M) * det23 - get<2, 1>(M) * det13 + get<3, 1>(M) * det12;

    double det = get<0, 0>(M) * det123 - get<1, 0>(M) * det023 + get<2, 0>(M) * det013 - get<3, 0>(M) * det012;

    // other auxiliary computations for traceH3
    double prod0x1xs = get<0, 1>(M) * get<1, 1>(M) + get<0, 2>(M) * get<1, 2>(M) + get<0, 3>(M) * get<1, 3>(M);
    double prod0x2xs = get<0, 1>(M) * get<2, 1>(M) + get<0, 2>(M) * get<2, 2>(M) + get<0, 3>(M) * get<2, 3>(M);
    double prod0x3xs = get<0, 1>(M) * get<3, 1>(M) + get<0, 2>(M) * get<3, 2>(M) + get<0, 3>(M) * get<3, 3>(M);
    double with_pairs_mixed = get<1, 0>(M) * prod0x1xs + get<2, 0>(M) * prod0x2xs + get<3, 0>(M) * prod0x3xs;

    double traceH3 = get<0, 0>(M) * (1.5 * traceH2 - get<0, 0>(squares)) + 3 * (with_pairs_mixed + det123);

    // auxiliary computation for traceH4
    double row_squaresums_sq = pow2(row_squaresums[0]) + pow2(row_squaresums[1]) + pow2(row_squaresums[2]);
    double anticommuting_with0 = get<1, 0>(squares) * (row_squaresums[1] + row_squaresums[2]) + get<2, 0>(squares) * (row_squaresums[0] + row_squaresums[2]) + get<3, 0>(squares) * (row_squaresums[0] + row_squaresums[1]);
    double prod0x1x = get<0, 0>(M) * get<1, 0>(M) + prod0x1xs;
    double prod0x2x = get<0, 0>(M) * get<2, 0>(M) + prod0x2xs;
    double prod0x3x = get<0, 0>(M) * get<3, 0>(M) + prod0x3xs;
    double prod1x2x = get<1, 0>(M) * get<2, 0>(M) - get<1, 1>(M) * get<2, 1>(M) - get<1, 2>(M) * get<2, 2>(M) - get<1, 3>(M) * get<2, 3>(M);
    double prod1x3x = get<1, 0>(M) * get<3, 0>(M) - get<1, 1>(M) * get<3, 1>(M) - get<1, 2>(M) * get<3, 2>(M) - get<1, 3>(M) * get<3, 3>(M);
    double prod2x3x = get<2, 0>(M) * get<3, 0>(M) - get<2, 1>(M) * get<3, 1>(M) - get<2, 2>(M) * get<3, 2>(M) - get<2, 3>(M) * get<3, 3>(M);
    double halfmixed_products = -pow2(prod0x1x) - pow2(prod0x2x) - pow2(prod0x3x) + pow2(prod1x2x) + pow2(prod1x3x) + pow2(prod2x3x);
  
    double traceH4 = -2 * det + 0.75 * pow2(traceH2) - 0.5 * row_squaresums_sq - anticommuting_with0 - squares0x * (0.5 * squares0x + corner_squaresum) + 4 * get<0, 0>(M) * (2 * det123 + with_pairs_mixed) - halfmixed_products - get<0, 0>(squares) * (squaresx0 + 0.5 * get<0, 0>(squares)) - 0.5 * (pow2(get<1, 0>(squares)) + pow2(get<2, 0>(squares)) + pow2(get<3, 0>(squares)));

    // converting power sums to symmetric polynomials
    double elem1 = traceH;
    double elem2 = 0.5 * (elem1 * traceH - traceH2);
    double elem3 = (elem2 * traceH - elem1 * traceH2 + traceH3) / 3.;
    double elem4 = 0.25 * (elem3 * traceH - elem2 * traceH2 + elem1 * traceH3 - traceH4);

    return ((elem1 >= 0) && (elem2 >= 0) && (elem3 >= 0) && (elem4 >= 0));
}

// --------------------------------
// Generic runners for PSD checkers
// --------------------------------

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

// --------------
// Interface code
// --------------

int main(int argc, char** argv) {

    // Available PSD checkers with their descriptions
    std::vector<std::pair<std::string, PSDChecker*> > checkers{
        std::make_pair(std::string("Charpoly avoiding H"), &check_with_charpoly_noh),
        std::make_pair(std::string("Sylvester criterion avoiding H"), &check_with_sylvester_noh),
        std::make_pair(std::string("Eigen's self-adjoint eigensolver"), &check_eignevalues),
        std::make_pair(std::string("Eigen's Choletsky"), &check_choletsky),
        std::make_pair(std::string("Formula via Pauli matrices"), &check_pauli)
    };

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

    std::cout << "Checking correctness ";
    bool correct = test_correctness(matrices, NUM_MATRICES, result, checkers[0].second);
    if (correct) {
        std::cout << "OK" << std::endl;
    } else {
        std::cout << "Not OK" << std::endl;
    }

    free(matrices);
    free(result);
    return 0;
}
