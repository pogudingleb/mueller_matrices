// starting version of this code was generated using Claude.ai

#include "main.h" 

// Utility function to find one real root of a cubic equation using Cardano's method
// under assumption that all the roots are real
// For equation x^3 + p*x^2 + q*x + r = 0
double solve_cubic_real_cardano(double p, double q, double r) {
    // Convert to depressed cubic: t^3 + at + b = 0, where x = t - p/3
    double a = q - p * p / 3.0;
    double b = r - p * q / 3.0 + 2.0 * p * p * p / 27.0;
    
    // Discriminant
    double discriminant = -(4.0 * a*a*a + 27.0 * b * b); 

    if (std::abs(discriminant) < 1e-16) { // Multiple roots
        if (std::abs(a) < 1e-12) { // Triple root
            return -p / 3.0;
        } else { // One single + one double root
            return 3.0 * b/(a) - p / 3.0;
        }
    }
    
    assert(discriminant > 0);
    double amplitude = 2.0 * std::sqrt(-a / 3.0);
    double angle = std::acos(3.0 * b / (a * amplitude)) / 3.0;
    return amplitude * std::cos(angle) - p / 3.0;
}

double solve_cubic_real(double p, double q, double r) {
    // double curr_point = std::abs(p) + std::abs(q) + std::abs(r) + 1.0;
    double curr_point = solve_cubic_real_cardano(p, q, r);
    double prev_point = 2.0 * curr_point;
    while (std::abs(prev_point - curr_point) > 1e-10) {
        prev_point = curr_point;
        curr_point = curr_point - (curr_point * (curr_point * (curr_point + p) + q) + r) / (curr_point * (3.0 * curr_point + 2.0 * p) + q);
    }
    return curr_point;
}




// Ferrari's method for solving quartic equation t^4 + a3*t^3 + a2*t^2 + a1*t + a0 = 0
// assuming that all the roots are real
void ferrari_method(double* coeffs, double* result) {
    // Convert to depressed quartic y^4 + py^2 + qy + r = 0
    // Substitution: t = y - a3/4
    double tmp1 = coeffs[0] * coeffs[0] / 8.0;
    double p = coeffs[1] - 3.0 * tmp1;
    double q = coeffs[2] - coeffs[0] * coeffs[1] / 2.0 + coeffs[0] * tmp1;
    double r = coeffs[3] - coeffs[0] * coeffs[2] / 4.0 + coeffs[1] * tmp1 * 0.5 - 3.0 * tmp1 * tmp1 / 4.0;
   

    if (std::abs(q) < 1e-14) {
        // biquadratic case
        double disc = p * p - 4.0 * r;
        double r1 = 0.5 * (-p + std::sqrt(disc));
        double r2 = -p - r1;
        result[0] = std::sqrt(r1);
        result[1] = -result[0];
        result[2] = std::sqrt(r2);
        result[3] = -result[2];
    } else {
        // using the resolvent cubic
        double alphasq = solve_cubic_real(2.0 * p, p * p - 4.0 * r, -q * q);
        double alpha = std::sqrt(alphasq);
        double beta = 0.5 * (alphasq + p - q / alpha);
        double gamma = alphasq + p - beta;
        
        double disc1 = alpha * alpha - 4.0 * beta;
        double disc2 = alpha * alpha - 4.0 * gamma;
        result[0] = 0.5 * (-alpha - std::sqrt(disc1));
        result[1] = -alpha - result[0];
        result[2] = 0.5 * (alpha - std::sqrt(disc2));
        result[3] = alpha - result[2];
    }
     
    for (size_t i = 0; i < 4; ++i) {
        result[i] -= coeffs[0] / 4.0;
    }

    std::sort(result, result + 4);
}
