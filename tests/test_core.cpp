#include "qtransport.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <numeric>

void require(bool condition, const char* message) {
    if (!condition) {
        std::cerr << "Failed: " << message << '\n';
        std::abort();
    }
}

void near(double a, double b, double tol = 1e-10) {
    if (std::abs(a - b) > tol) {
        std::cerr << "Mismatch: " << a << " vs " << b << '\n';
        std::abort();
    }
}

int main() {
    // No disorder: transmission is unity.
    for (double E : {-1.6, -0.5, 0.0, 0.7, 1.8}) {
        const auto y = qt::log_transmission(std::vector<double>(500, 0.0), E);
        for (auto v : y) near(v, 0.0, 1e-10);
    }
    // One impurity: T(E)=4 sin^2(k)/(4 sin^2(k)+U^2).
    for (double E : {-1.0, 0.0, 0.8}) {
        const double U = 0.4;
        const auto y = qt::log_transmission({U}, E);
        near(y.back(), std::log((4.0 - E * E) / (4.0 - E * E + U * U)));
    }
    // Binary weights and block construction.
    qt::Config c;
    c.model = "binary";
    c.s1 = 2;
    c.s2 = 5;
    c.p = 1.0;
    qt::SpacingDistribution binary(c);
    near(binary.probability(0), 1.0);
    std::mt19937_64 rng(4);
    const auto chain = qt::build_chain(6, 0.3, binary, rng);
    require((chain == std::vector<double>{0, 0.3, 0, 0.3, 0, 0.3}), "Binary chain");

    // Both generalized schemes agree with standard zeta at s0=d=1.
    c.model = "zeta";
    c.max_rank = 250;
    qt::SpacingDistribution standard(c);
    c.model = "zeta-rank";
    qt::SpacingDistribution rank_standard(c);
    c.model = "zeta-value";
    qt::SpacingDistribution value_standard(c);
    for (std::size_t i = 0; i < standard.size(); ++i) {
        near(standard.probability(i), rank_standard.probability(i));
        near(standard.probability(i), value_standard.probability(i));
    }
    c.s0 = 1;
    c.step = 2;
    qt::SpacingDistribution value_odd(c);
    c.model = "zeta-rank";
    qt::SpacingDistribution rank_odd(c);
    require(value_odd.spacing(1) == 3 && rank_odd.spacing(1) == 3, "Odd spacing");
    near(rank_odd.probability(1) / rank_odd.probability(0), 1.0 / 8.0);
    near(value_odd.probability(1) / value_odd.probability(0), 1.0 / 27.0);

    // The exponential law has the expected adjacent-weight ratio.
    c.model = "exponential";
    c.mu = 0.5;
    qt::SpacingDistribution exponential(c);
    near(exponential.probability(1) / exponential.probability(0), std::exp(-0.5));
    const auto y = qt::log_transmission(std::vector<double>(300, 0.2), 0.2);
    require(std::all_of(y.begin(), y.end(), [](double x) { return std::isfinite(x); }), "Finite trace");
    std::cout << "All numerical checks passed.\n";
}
