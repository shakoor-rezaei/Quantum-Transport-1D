#include "qtransport.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace qt {

std::vector<double> build_chain(int sites, double potential,
                                const SpacingDistribution& dist,
                                std::mt19937_64& rng) {
    if (sites < 1) throw std::invalid_argument("sites must be positive");
    std::vector<double> chain(static_cast<std::size_t>(sites), 0.0);
    std::int64_t index = -1;
    while (index < sites) {
        const auto s = dist.draw(rng);
        if (s > static_cast<std::int64_t>(sites) - 1 - index) break;
        index += s;
        chain[static_cast<std::size_t>(index)] = potential;
    }
    return chain;
}

std::vector<double> log_transmission(const std::vector<double>& chain, double E) {
    if (!std::isfinite(E) || E <= -2 || E >= 2)
        throw std::invalid_argument("Energy must satisfy -2 < E < 2");
    const double k = std::acos(E / 2.0);
    const double log_numerator = std::log(4.0 * std::pow(std::sin(k), 2));
    const std::complex<double> phase(std::cos(k), -std::sin(k));
    std::complex<double> c0(std::cos(k), std::sin(k));
    std::complex<double> c1(1.0, 0.0);
    double log_scale = 0.0;
    std::vector<double> trace(chain.size() + 1);

    for (std::size_t n = 0; n <= chain.size(); ++n) {
        const double denom = std::norm(phase * c1 - c0);
        if (!(denom > 0) || !std::isfinite(denom))
            throw std::runtime_error("Numerical instability in transmission denominator");
        trace[n] = log_numerator - std::log(denom) - 2.0 * log_scale;
        if (n == chain.size()) break;
        const auto c2 = (E - chain[n]) * c1 - c0;
        c0 = c1;
        c1 = c2;
        const double scale = std::max(std::abs(c0), std::abs(c1));
        if (!(scale > 0) || !std::isfinite(scale))
            throw std::runtime_error("Numerical instability in recurrence");
        if (scale > 1e100 || scale < 1e-100) {
            c0 /= scale;
            c1 /= scale;
            log_scale += std::log(scale);
        }
    }
    return trace;
}

Fit linear_fit(const std::vector<double>& y, int begin, int end) {
    if (begin < 0 || end > static_cast<int>(y.size()) || end - begin < 2)
        throw std::invalid_argument("Invalid fitting interval");
    const long double count = end - begin;
    long double sx = 0, sy = 0, sxx = 0, sxy = 0;
    for (int i = begin; i < end; ++i) {
        const long double x = i, v = y[static_cast<std::size_t>(i)];
        sx += x;
        sy += v;
        sxx += x * x;
        sxy += x * v;
    }
    const long double denominator = count * sxx - sx * sx;
    if (!(denominator > 0)) throw std::runtime_error("Degenerate linear fit");
    const long double m = (count * sxy - sx * sy) / denominator;
    const long double b = (sy - m * sx) / count;
    long double residual = 0, total = 0;
    for (int i = begin; i < end; ++i) {
        const long double v = y[static_cast<std::size_t>(i)];
        residual += (v - (m * i + b)) * (v - (m * i + b));
        total += (v - sy / count) * (v - sy / count);
    }
    return {static_cast<double>(m), static_cast<double>(b),
            static_cast<double>(total > 0 ? 1 - residual / total : 1)};
}

std::uint64_t splitmix64(std::uint64_t x) {
    x += 0x9e3779b97f4a7c15ULL;
    x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
    x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
    return x ^ (x >> 31);
}

} // namespace qt
