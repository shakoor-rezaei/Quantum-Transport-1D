#include "qtransport.hpp"

#include <algorithm>
#include <atomic>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <mutex>
#include <sstream>
#include <stdexcept>
#include <thread>

namespace fs = std::filesystem;

namespace {

struct Result {
    double E = 0;
    std::vector<double> mean_lnT;
    qt::Fit fit{0, 0, 0};
    double gamma = 0;
    double xi = 0;
};

void usage() {
    std::cout << R"(Quantum transport in a disordered 1D tight-binding chain

  --model          binary | exponential | zeta | zeta-rank | zeta-value
  --N              Number of sites (default 1000)
  --U              Impurity on-site potential (default 0.1)
  --s1 --s2 --p    Binary spacings and P(s1) (defaults 2, 5, 0.5)
  --mu             Exponential decay parameter (default 0.5)
  --alpha          Zeta tail exponent (default 2.0)
  --s0 --d         Modified zeta: spacing s_r=s0+(r-1)d
  --max-rank       Truncation of exponential/zeta support (default 100000)
  --E              Single energy, or use --emin --emax --dE
  --emin --emax    Energy range, inclusive (defaults -1.9,1.9)
  --dE            Energy step (default 0.1)
  --realizations   Independent chains per energy (default 1000)
  --threads        Concurrent energies (default 1)
  --seed           Random seed (default 12345)
  --fit-start      Fit start fraction of N (default 0)
  --fit-end        Fit end fraction of N (default 1)
  --out            Output directory (default results)
  --overwrite      Replace outputs from an earlier run in that directory
  --help           Display help

The output includes <ln T_n> for n=0,...,N and a fit for xi=-2/slope.
)";
}

qt::Config parse(int argc, char* argv[]) {
    qt::Config c;
    for (int i = 1; i < argc; ++i) {
        const std::string key = argv[i];
        if (key == "--help" || key == "-h") { usage(); std::exit(0); }
        if (key == "--overwrite") { c.overwrite = true; continue; }
        if (i + 1 == argc) throw std::invalid_argument("Missing value for " + key);
        const std::string val = argv[++i];
        if      (key == "--model")        c.model = val;
        else if (key == "--N")            c.sites = std::stoi(val);
        else if (key == "--U")            c.potential = std::stod(val);
        else if (key == "--s1")           c.s1 = std::stoi(val);
        else if (key == "--s2")           c.s2 = std::stoi(val);
        else if (key == "--p")            c.p = std::stod(val);
        else if (key == "--mu")           c.mu = std::stod(val);
        else if (key == "--alpha")        c.alpha = std::stod(val);
        else if (key == "--s0")           c.s0 = std::stoi(val);
        else if (key == "--d")            c.step = std::stoi(val);
        else if (key == "--max-rank")     c.max_rank = std::stoi(val);
        else if (key == "--realizations") c.realizations = std::stoi(val);
        else if (key == "--threads")      c.threads = std::stoi(val);
        else if (key == "--seed")         c.seed = std::stoull(val);
        else if (key == "--fit-start")    c.fit_start = std::stod(val);
        else if (key == "--fit-end")      c.fit_end = std::stod(val);
        else if (key == "--E")            { c.single_energy = true; c.energy = std::stod(val); }
        else if (key == "--emin")         c.emin = std::stod(val);
        else if (key == "--emax")         c.emax = std::stod(val);
        else if (key == "--dE")           c.de = std::stod(val);
        else if (key == "--out")          c.out_dir = val;
        else throw std::invalid_argument("Unknown option: " + key);
    }
    if (c.sites < 2 || c.realizations < 1 || c.threads < 1)
        throw std::invalid_argument("Require N>=2, realizations>=1, threads>=1");
    if (!std::isfinite(c.potential) ||
        !(0 <= c.fit_start && c.fit_start < c.fit_end && c.fit_end <= 1))
        throw std::invalid_argument("Invalid potential or fit interval");
    if (c.single_energy) {
        if (!std::isfinite(c.energy) || c.energy <= -2 || c.energy >= 2)
            throw std::invalid_argument("Require -2 < E < 2");
    } else if (!std::isfinite(c.emin) || !std::isfinite(c.emax) ||
               !std::isfinite(c.de) || c.emin <= -2 || c.emax >= 2 ||
               c.emin > c.emax || c.de <= 0) {
        throw std::invalid_argument("Invalid energy grid: require -2<emin<=emax<2 and dE>0");
    }
    if (c.out_dir.empty()) throw std::invalid_argument("Output directory is empty");
    return c;
}

std::vector<double> energy_grid(const qt::Config& c) {
    if (c.single_energy) return {c.energy};
    std::vector<double> energies;
    const long double span = static_cast<long double>(c.emax) - c.emin;
    const long double steps = span / c.de;
    if (steps > 1000000) throw std::invalid_argument("Energy grid exceeds 1,000,001 points");
    for (std::uint64_t i = 0; i <= static_cast<std::uint64_t>(std::floor(steps + 1e-10L)); ++i) {
        const long double x = static_cast<long double>(c.emin) + i * static_cast<long double>(c.de);
        if (x <= c.emax + 1e-10L) energies.push_back(static_cast<double>(x));
    }
    if (energies.empty()) throw std::runtime_error("Empty energy grid");
    return energies;
}

Result simulate(double E, const qt::Config& c, const qt::SpacingDistribution& dist) {
    Result result;
    result.E = E;
    result.mean_lnT.assign(static_cast<std::size_t>(c.sites + 1), 0.0);
    std::vector<long double> sums(static_cast<std::size_t>(c.sites + 1), 0.0L);
    for (int realization = 0; realization < c.realizations; ++realization) {
        // A realization has the same disorder seed at every energy.
        std::mt19937_64 rng(qt::splitmix64(c.seed + static_cast<std::uint64_t>(realization)));
        const auto chain = qt::build_chain(c.sites, c.potential, dist, rng);
        const auto trace = qt::log_transmission(chain, E);
        for (std::size_t n = 0; n < trace.size(); ++n) sums[n] += trace[n];
    }
    for (std::size_t n = 0; n < sums.size(); ++n)
        result.mean_lnT[n] = static_cast<double>(sums[n] / c.realizations);
    const int begin = std::max(0, std::min(c.sites - 1,
        static_cast<int>(std::floor(c.fit_start * c.sites))));
    const int end = std::max(begin + 2, std::min(c.sites + 1,
        static_cast<int>(std::ceil(c.fit_end * c.sites)) + 1));
    result.fit = qt::linear_fit(result.mean_lnT, begin, end);
    result.gamma = -0.5 * result.fit.slope;
    constexpr double tolerance = 1e-13;
    result.xi = result.gamma > tolerance ? 1.0 / result.gamma :
        (result.gamma >= -tolerance ? std::numeric_limits<double>::infinity()
                                    : std::numeric_limits<double>::quiet_NaN());
    return result;
}

void write_outputs(const qt::Config& c, const qt::SpacingDistribution& dist,
                   const std::vector<Result>& results) {
    const fs::path base = c.out_dir;
    fs::create_directories(base);
    const auto summary_path = base / "summary.csv";
    if (fs::exists(summary_path) && !c.overwrite)
        throw std::runtime_error("Output exists. Choose a new --out directory or use --overwrite");
    auto open = [](const fs::path& path) {
        std::ofstream f(path);
        if (!f) throw std::runtime_error("Cannot write: " + path.string());
        return f;
    };
    auto cfg_file = open(base / "run_config.txt");
    cfg_file << std::setprecision(17)
             << "model=" << c.model << "\nN=" << c.sites << "\nU=" << c.potential
             << "\ns1=" << c.s1 << "\ns2=" << c.s2 << "\np=" << c.p
             << "\nmu=" << c.mu << "\nalpha=" << c.alpha << "\ns0=" << c.s0
             << "\nd=" << c.step << "\nmax_rank=" << c.max_rank
             << "\ntruncated_mean_spacing=" << dist.mean_spacing()
             << "\nrealizations=" << c.realizations << "\nseed=" << c.seed
             << "\nthreads=" << c.threads << "\nfit_start=" << c.fit_start
             << "\nfit_end=" << c.fit_end << "\n";
    auto preview = open(base / "distribution_preview.csv");
    preview << "rank,spacing,probability\n" << std::setprecision(17);
    for (std::size_t i = 0; i < std::min<std::size_t>(dist.size(), 40); ++i)
        preview << i + 1 << ',' << dist.spacing(i) << ',' << dist.probability(i) << '\n';
    std::mt19937_64 sample_rng(qt::splitmix64(c.seed));
    const auto sample = qt::build_chain(c.sites, c.potential, dist, sample_rng);
    auto sample_file = open(base / "sample_chain.csv");
    sample_file << "site,potential\n" << std::setprecision(17);
    for (std::size_t i = 0; i < std::min<std::size_t>(sample.size(), 500); ++i)
        sample_file << i << ',' << sample[i] << '\n';

    auto summary = open(summary_path);
    summary << "energy,file,N,U,realizations,mean_lnT_N,slope,gamma,xi,r_squared,fit_start_site,fit_end_site\n";
    summary << std::setprecision(17);
    const auto start = std::max(0, std::min(c.sites - 1,
        static_cast<int>(std::floor(c.fit_start * c.sites))));
    const auto end = std::max(start + 2, std::min(c.sites + 1,
        static_cast<int>(std::ceil(c.fit_end * c.sites)) + 1));
    for (std::size_t i = 0; i < results.size(); ++i) {
        std::ostringstream name;
        name << "energy_" << std::setfill('0') << std::setw(5) << i << ".dat";
        auto data = open(base / name.str());
        data << std::setprecision(17);
        for (double y : results[i].mean_lnT) data << y << '\n';
        const auto& r = results[i];
        summary << r.E << ',' << name.str() << ',' << c.sites << ',' << c.potential
                << ',' << c.realizations << ',' << r.mean_lnT.back() << ','
                << r.fit.slope << ',' << r.gamma << ',' << r.xi << ',' << r.fit.r2
                << ',' << start << ',' << end - 1 << '\n';
    }
}

} // namespace

int main(int argc, char* argv[]) {
    try {
        if (argc == 1) { usage(); return 0; }
        const auto c = parse(argc, argv);
        const qt::SpacingDistribution dist(c);
        const auto energies = energy_grid(c);
        if (fs::exists(fs::path(c.out_dir) / "summary.csv") && !c.overwrite)
            throw std::runtime_error("Output exists. Use --out or --overwrite");
        std::vector<Result> results(energies.size());
        std::atomic<std::size_t> next{0};
        std::exception_ptr error;
        std::mutex error_mutex;
        auto worker = [&] {
            while (true) {
                const auto i = next.fetch_add(1);
                if (i >= energies.size()) break;
                try { results[i] = simulate(energies[i], c, dist); }
                catch (...) {
                    std::lock_guard<std::mutex> lock(error_mutex);
                    if (!error) error = std::current_exception();
                    next.store(energies.size());
                    break;
                }
            }
        };
        std::vector<std::thread> pool;
        for (int i = 0; i < std::min<int>(c.threads, energies.size()); ++i)
            pool.emplace_back(worker);
        for (auto& thread : pool) thread.join();
        if (error) std::rethrow_exception(error);
        write_outputs(c, dist, results);
        std::cout << "Model: " << c.model << " | energies: " << energies.size()
                  << " | mean spacing (truncated): " << std::setprecision(8)
                  << dist.mean_spacing() << "\nResults: " << c.out_dir << "/summary.csv\n";
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << '\n';
        return 1;
    }
}
