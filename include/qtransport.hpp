#pragma once

#include <complex>
#include <cstddef>
#include <cstdint>
#include <random>
#include <string>
#include <vector>

namespace qt {

struct Config {
    std::string model = "binary";
    int sites = 1000;
    double potential = 0.1;
    int s1 = 2;
    int s2 = 5;
    double p = 0.5;
    double mu = 0.5;
    double alpha = 2.0;
    int s0 = 1;
    int step = 1;
    int max_rank = 100000;
    int realizations = 1000;
    int threads = 1;
    std::uint64_t seed = 12345;
    double fit_start = 0.0;
    double fit_end = 1.0;
    double emin = -1.9;
    double emax = 1.9;
    double de = 0.1;
    bool single_energy = false;
    double energy = 0.0;
    bool overwrite = false;
    std::string out_dir = "results";
};

class SpacingDistribution {
public:
    explicit SpacingDistribution(const Config& cfg);
    std::int64_t draw(std::mt19937_64& rng) const;
    std::size_t size() const { return support_.size(); }
    std::int64_t spacing(std::size_t i) const { return support_.at(i); }
    double probability(std::size_t i) const { return probabilities_.at(i); }
    double mean_spacing() const;

private:
    std::vector<std::int64_t> support_;
    std::vector<double> probabilities_;
    std::vector<double> cdf_;
};

std::vector<double> build_chain(int sites, double potential,
                                const SpacingDistribution& distribution,
                                std::mt19937_64& rng);
std::vector<double> log_transmission(const std::vector<double>& chain, double energy);

struct Fit {
    double slope;
    double intercept;
    double r2;
};

Fit linear_fit(const std::vector<double>& y, int begin, int end);
std::uint64_t splitmix64(std::uint64_t x);

} // namespace qt
