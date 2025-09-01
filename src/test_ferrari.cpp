// The original version of this testing code was written by Claude.ai
#include "main.h"

// tests cubic solver for a polynomial with the prescribed roots
bool test_cubic(double r1, double r2, double r3) {
    double p = -(r1 + r2 + r3);
    double q = r1 * (r2 + r3) + r2 * r3;
    double r = -r1 * r2 * r3;

    double sol = solve_cubic_real(p, q, r);
    double err = std::min(std::abs(r1 - sol), std::min(std::abs(r2 - sol), std::abs(r3 - sol)));
    return err < 1e-10;   
}

bool test_quartic(const std::vector<double>& roots) {
    double coeffs[4];
    coeffs[0] = -(roots[0] + roots[1] + roots[2] + roots[3]);
    coeffs[1] = roots[0] * (roots[1] + roots[2] + roots[3]) + roots[1] * (roots[2] + roots[3]) + roots[2] * roots[3];
    coeffs[2] = -(roots[0] * roots[1] * (roots[3] + roots[2]) + roots[2] * roots[3] * (roots[1] + roots[0]));
    coeffs[3] = roots[0] * roots[1] * roots[2] * roots[3];

    double result[4];
    ferrari_method(coeffs, result);

    bool res = true;
    for (size_t i = 0; i < 4; ++i) {
            std::cout << roots[i] << " " << result[i] << std::endl;
        if (std::abs(roots[i] - result[i]) > 1e-10) {
            res = false;
        }
    }
    return res;
}

int main() {
    std::cout << "Testing cubic solver" << std::endl;
    std::vector<std::vector<double> > test_cases = {
        {1., 2., 3.},
//        {1., 1., 1.},
        {1., -2., -2.},
        {-0.5632764, 0.83647, 0.83648},
        {-0.3239842, 1.947853982, 1000.32874938},
        {37478.0, 3874.0, -0.0001}
    };
    for (auto c: test_cases) {
        if (test_cubic(c[0], c[1], c[2])) {
            std::cout << "OK" << std::endl;
        } else {
            std::cout << "NOK" << std::endl;
        }
    }

    std::cout << "Testing quartic solver" << std::endl;
    std::vector<std::vector<double> > test_cases_quartic = {
        {1., 2., 3., 4.},
        {1., 1., 1., 1.},
        {-3., -2., -2., -1.},
        {-0.5632764, 0.83647, 0.83648, 0.9},
        {-0.3239842, 1.947853982, 1000.3287493, 2000.}
    };
    for (auto c: test_cases_quartic) {
        if (test_quartic(c)) {
            std::cout << "OK" << std::endl;
        } else {
            std::cout << "NOK" << std::endl;
        }
    }


    return 0;
}
/*
// Test function to verify polynomial with given roots
void test_polynomial(const std::vector<double>& expected_roots, const std::string& description) {
    std::cout << "\n" << description << "\n";
    std::cout << std::string(description.length(), '=') << "\n";
    
    // Construct polynomial coefficients from roots
    // For roots r1, r2, r3, r4: (t-r1)(t-r2)(t-r3)(t-r4) = t^4 + a3*t^3 + a2*t^2 + a1*t + a0
    double a3 = -(expected_roots[0] + expected_roots[1] + expected_roots[2] + expected_roots[3]);
    double a2 = expected_roots[0]*expected_roots[1] + expected_roots[0]*expected_roots[2] + expected_roots[0]*expected_roots[3] +
                expected_roots[1]*expected_roots[2] + expected_roots[1]*expected_roots[3] + expected_roots[2]*expected_roots[3];
    double a1 = -(expected_roots[0]*expected_roots[1]*expected_roots[2] + expected_roots[0]*expected_roots[1]*expected_roots[3] +
                  expected_roots[0]*expected_roots[2]*expected_roots[3] + expected_roots[1]*expected_roots[2]*expected_roots[3]);
    double a0 = expected_roots[0]*expected_roots[1]*expected_roots[2]*expected_roots[3];
    
    std::cout << "Polynomial: t^4 + (" << a3 << ")t^3 + (" << a2 << ")t^2 + (" << a1 << ")t + (" << a0 << ")\n";
    std::cout << "Expected roots: ";
    for (size_t i = 0; i < expected_roots.size(); ++i) {
        std::cout << expected_roots[i];
        if (i < expected_roots.size() - 1) std::cout << ", ";
    }
    std::cout << "\n\n";
    
    // Solve using Ferrari's method
    std::vector<std::complex<double>> computed_roots = ferrari_method(a3, a2, a1, a0);
    
    std::cout << "Computed roots:\n";
    for (size_t i = 0; i < computed_roots.size(); ++i) {
        std::cout << "  Root " << i+1 << ": " << std::fixed << std::setprecision(8) 
             << computed_roots[i].real();
        if (abs(computed_roots[i].imag()) > 1e-12) {
                std::cout << " + " << computed_roots[i].imag() << "i";
        }
        std::cout << "\n";
    }
    
    // Verification: substitute each root back into the polynomial
    std::cout << "\nVerification (should be close to 0):\n";
    for (size_t i = 0; i < computed_roots.size(); ++i) {
        std::complex<double> t = computed_roots[i];
        std::complex<double> result = t*t*t*t + a3*t*t*t + a2*t*t + a1*t + a0;
        std::cout << "  f(root " << i+1 << ") = " << std::fixed << std::setprecision(2) << std::scientific 
             << abs(result) << "\n";
    }
}

int main() {
    std::cout << "Ferrari's Method for Solving Quartic Equations\n";
    std::cout << "==============================================\n";
    
    // Test 1: Simple integer roots
    test_polynomial({1.0, 2.0, 3.0, 4.0}, "Test 1: Simple integer roots (1, 2, 3, 4)");
    
    // Test 2: Roots with some negatives
    test_polynomial({-1.0, 2.0, -3.0, 4.0}, "Test 2: Mixed sign roots (-1, 2, -3, 4)");
    
    // Test 3: Repeated roots
    test_polynomial({1.0, 1.0, 2.0, 2.0}, "Test 3: Repeated roots (1, 1, 2, 2)");
    
    // Test 4: One repeated root
    test_polynomial({0.0, 0.0, 0.0, 5.0}, "Test 4: Triple zero root and one non-zero (0, 0, 0, 5)");
    
    // Test 5: Fractional roots
    test_polynomial({0.5, -0.5, 1.5, -1.5}, "Test 5: Fractional roots (0.5, -0.5, 1.5, -1.5)");

    // Test 6: Failing the first version
    test_polynomial({0.05367587107, 0.06088803772, 0.07241383492, 0.1919622563}, "Test 6: failing the first version");
    
    return 0;
}
*/
