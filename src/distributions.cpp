#include "qtransport.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace qt {

SpacingDistribution::SpacingDistribution(const Config& c) {
    if (c.model == "binary") {
        if (c.s1 < 1 || c.s2 < 1 || c.s1 == c.s2 ||
            !std::isfinite(c.p) || c.p < 0 || c.p > 1)
            throw std::invalid_argument("Binary model: require distinct s1,s2 >= 1 and 0 <= p <= 1");
        support_ = {c.s1, c.s2};
        probabilities_ = {c.p, 1.0 - c.p};
    } else if (c.model == "exponential") {
        if (!std::isfinite(c.mu) || c.mu <= 0)
            throw std::invalid_argument("Exponential model: mu must be positive");
        if (c.max_rank < 1) throw std::invalid_argument("max-rank must be positive");
        support_.reserve(c.max_rank);
        probabilities_.reserve(c.max_rank);
        // P(s) = (1-exp(-mu))*exp(-mu*(s-1)), truncated at max-rank.
        const long double q = std::exp(-static_cast<long double>(c.mu));
        const long double norm = -std::expm1(-static_cast<long double>(c.mu) * c.max_rank);
        const long double first = -std::expm1(-static_cast<long double>(c.mu));
        for (int s = 1; s <= c.max_rank; ++s) {
            support_.push_back(s);
            probabilities_.push_back(static_cast<double>(first * std::pow(q, s - 1) / norm));
        }
    } else if (c.model == "zeta" || c.model == "zeta-rank" || c.model == "zeta-value") {
        if (!(c.alpha > 0) || !std::isfinite(c.alpha) || c.max_rank < 1)
            throw std::invalid_argument("Zeta model: alpha > 0 and max-rank >= 1 required");
        if (c.s0 < 1 || c.step < 1)
            throw std::invalid_argument("Zeta model: s0 and d must be positive integers");
        if (c.model == "zeta" && (c.s0 != 1 || c.step != 1))
            throw std::invalid_argument("Use zeta-rank or zeta-value for modified spacing support");
        support_.reserve(c.max_rank);
        probabilities_.reserve(c.max_rank);
        const long double power = 1.0L + c.alpha;
        std::vector<long double> weights;
        weights.reserve(c.max_rank);
        long double total = 0;
        for (int i = 0; i < c.max_rank; ++i) {
            const auto spacing = static_cast<std::int64_t>(c.s0) +
                                 static_cast<std::int64_t>(i) * c.step;
            if (spacing < 1 || spacing > std::numeric_limits<int>::max())
                throw std::invalid_argument("Modified spacing exceeds supported integer range");
            const long double base = (c.model == "zeta-value")
                ? static_cast<long double>(spacing)
                : static_cast<long double>(i + 1);
            const long double w = std::pow(base, -power);
            support_.push_back(spacing);
            weights.push_back(w);
            total += w;
        }
        for (const auto w : weights) probabilities_.push_back(static_cast<double>(w / total));
    } else {
        throw std::invalid_argument("Unknown model: " + c.model);
    }

    long double running = 0;
    cdf_.reserve(probabilities_.size());
    for (double p : probabilities_) {
        running += p;
        cdf_.push_back(static_cast<double>(running));
    }
    cdf_.back() = 1.0;
}

std::int64_t SpacingDistribution::draw(std::mt19937_64& rng) const {
    std::uniform_real_distribution<double> uniform(0.0, 1.0);
    const double u = uniform(rng);
    const auto pos = std::lower_bound(cdf_.begin(), cdf_.end(), u);
    const auto i = static_cast<std::size_t>(pos - cdf_.begin());
    return support_[std::min(i, support_.size() - 1)];
}

double SpacingDistribution::mean_spacing() const {
    long double mean = 0;
    for (std::size_t i = 0; i < support_.size(); ++i)
        mean += static_cast<long double>(support_[i]) * probabilities_[i];
    return static_cast<double>(mean);
}

} // namespace qt
