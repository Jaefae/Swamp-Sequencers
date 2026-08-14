// SwampBench — suffix array vs. suffix tree
//
// Measures the three axes the two structures actually trade against each other:
// build time, index memory, and query latency.
//
// Output discipline: the CSV goes to stdout and nothing else ever does, so
//
//     SwampBench --max-size 200000 --steps 4 --queries 200 > results.csv
//
// captures exactly the data. Progress chatter goes to stderr, which is why CI
// can redirect stdout straight into an artifact.
//
// The CSV is "tidy" (one measurement per row) rather than one column per
// metric, so new metrics can be added without breaking anything already
// parsing it:
//
//     experiment,structure,x,metric,value
//
//   experiment – "scaling" (x = text length n) or "query" (x = pattern length m)
//   structure  – "array" or "tree"
//   metric     – build_ms | memory_bytes | query_us | total_hits | node_count
//
// total_hits is not a performance number; it is the sanity check. Every pattern
// is sampled *from the text being searched*, so a run in which some structure
// reports fewer hits than queries has a correctness bug, and a run in which the
// two structures disagree has a bug in one of them.

#include "SwampSeqLib/suffix_array.h"
#include "SwampSeqLib/suffix_tree.h"

#include <chrono>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <random>
#include <string>
#include <vector>

namespace {

using Clock = std::chrono::steady_clock;

// Tunables

// Pattern length used by the scaling sweep. Long enough that a sampled pattern
// is almost always unique in the text, so the scaling numbers measure the
// search itself rather than the cost of reporting a huge occurrence list.
constexpr size_t kScalingPatternLength = 20;

// Pattern lengths for the query sweep. Short patterns are match-list bound,
// long ones are descent bound; the crossover is the interesting part.
constexpr size_t kQueryLengths[] = {4, 8, 12, 16, 24, 32, 48, 64};

struct Config {
  std::string fasta;             // empty => generate a synthetic genome
  size_t maxSize = 1000000;      // largest n in the scaling sweep
  size_t steps = 5;              // number of points in the scaling sweep
  size_t queries = 1000;         // patterns timed per measurement
  uint32_t seed = 20240817;      // drives both text generation and sampling
};

// CSV

void emit(const char *experiment, const char *structure, size_t x,
          const char *metric, double value) {
  std::cout << experiment << ',' << structure << ',' << x << ',' << metric << ','
            << std::fixed << std::setprecision(4) << value << '\n';
}

// Input

// Read a FASTA file, dropping '>' header lines and keeping only ACGT. That also
// silently drops line endings, lowercase soft-masked regions' case, and the
// ambiguity codes (N and friends) that neither index treats specially.
bool loadFasta(const std::string &path, std::string &out) {
  std::ifstream in(path);
  if (!in) {
    std::cerr << "bench: cannot open FASTA file: " << path << '\n';
    return false;
  }
  std::string line;
  while (std::getline(in, line)) {
    if (!line.empty() && line[0] == '>')
      continue;
    for (char c : line) {
      switch (c) {
      case 'a': case 'A': out.push_back('A'); break;
      case 'c': case 'C': out.push_back('C'); break;
      case 'g': case 'G': out.push_back('G'); break;
      case 't': case 'T': out.push_back('T'); break;
      default: break;
      }
    }
  }
  return true;
}

// Uniform random ACGT. Real genomes are not uniform, which is exactly why the
// default is reproducible-synthetic and --fasta exists for the real thing.
std::string generateGenome(size_t n, std::mt19937 &rng) {
  static const char kBases[] = {'A', 'C', 'G', 'T'};
  std::uniform_int_distribution<int> pick(0, 3);
  std::string text;
  text.reserve(n);
  for (size_t i = 0; i < n; ++i)
    text.push_back(kBases[pick(rng)]);
  return text;
}

// Measurement

// Sample `count` patterns of length `m` from `text`. Sampling from the text is
// what makes total_hits meaningful (see the header comment).
std::vector<std::string> samplePatterns(const std::string &text, size_t m,
                                        size_t count, std::mt19937 &rng) {
  std::vector<std::string> pats;
  if (text.size() < m || count == 0)
    return pats;
  std::uniform_int_distribution<size_t> pick(0, text.size() - m);
  pats.reserve(count);
  for (size_t i = 0; i < count; ++i)
    pats.push_back(text.substr(pick(rng), m));
  return pats;
}

struct QueryStats {
  double usPerQuery = 0.0;
  size_t totalHits = 0;
};

// Time `search` across every pattern. The hit count is accumulated inside the
// timed loop on purpose: it consumes the returned vector, so the optimizer
// cannot discard the search as dead code.
template <typename Index>
QueryStats timeQueries(const Index &index, const std::vector<std::string> &pats) {
  QueryStats stats;
  if (pats.empty())
    return stats;

  const auto start = Clock::now();
  for (const std::string &p : pats)
    stats.totalHits += index.search(p).size();
  const auto elapsed = Clock::now() - start;

  const double us =
      std::chrono::duration<double, std::micro>(elapsed).count();
  stats.usPerQuery = us / static_cast<double>(pats.size());
  return stats;
}

double millisSince(Clock::time_point start) {
  return std::chrono::duration<double, std::milli>(Clock::now() - start).count();
}

// Experiments

// Build both indexes at increasing text lengths. This is the experiment behind
// the memory and build-time columns in the README.
void runScaling(const std::string &text, const Config &cfg, std::mt19937 &rng) {
  for (size_t step = 1; step <= cfg.steps; ++step) {
    const size_t n = text.size() * step / cfg.steps;
    if (n < kScalingPatternLength)
      continue;

    std::cerr << "bench: scaling n=" << n << " ...";
    const std::string sub = text.substr(0, n);
    const std::vector<std::string> pats =
        samplePatterns(sub, kScalingPatternLength, cfg.queries, rng);

    // Each index is built in its own scope and timed in place. Both classes
    // hold a raw pointer into their own owned copy of the text, so they are
    // deliberately never moved after construction.
    {
      const auto start = Clock::now();
      const SuffixArray sa{sub};
      emit("scaling", "array", n, "build_ms", millisSince(start));
      emit("scaling", "array", n, "memory_bytes",
           static_cast<double>(sa.memoryBytes()));
      const QueryStats q = timeQueries(sa, pats);
      emit("scaling", "array", n, "query_us", q.usPerQuery);
      emit("scaling", "array", n, "total_hits",
           static_cast<double>(q.totalHits));
    }
    {
      const auto start = Clock::now();
      const SuffixTree st{sub};
      emit("scaling", "tree", n, "build_ms", millisSince(start));
      emit("scaling", "tree", n, "memory_bytes",
           static_cast<double>(st.memoryBytes()));
      const QueryStats q = timeQueries(st, pats);
      emit("scaling", "tree", n, "query_us", q.usPerQuery);
      emit("scaling", "tree", n, "total_hits",
           static_cast<double>(q.totalHits));
      emit("scaling", "tree", n, "node_count",
           static_cast<double>(st.nodeCount()));
    }
    std::cerr << " done\n";
  }
}

// Sweep pattern length against a single pair of indexes built over the whole
// text, isolating query cost from build cost.
void runQuery(const std::string &text, const Config &cfg, std::mt19937 &rng) {
  std::cerr << "bench: building indexes over " << text.size() << " bases ...";
  const SuffixArray sa{text};
  const SuffixTree st{text};
  std::cerr << " done\n";

  for (size_t m : kQueryLengths) {
    if (m > text.size())
      break;

    std::cerr << "bench: query m=" << m << " ...";
    const std::vector<std::string> pats =
        samplePatterns(text, m, cfg.queries, rng);

    const QueryStats a = timeQueries(sa, pats);
    const QueryStats t = timeQueries(st, pats);

    emit("query", "array", m, "query_us", a.usPerQuery);
    emit("query", "tree", m, "query_us", t.usPerQuery);
    emit("query", "array", m, "total_hits", static_cast<double>(a.totalHits));
    emit("query", "tree", m, "total_hits", static_cast<double>(t.totalHits));
    std::cerr << " done\n";
  }
}

// Argument parsing

const char *kUsage =
    "usage: SwampBench [--fasta PATH] [--max-size N] [--steps N] "
    "[--queries N] [--seed N]";

// Distinguishes "print usage and stop successfully" from "bad arguments", so
// --help exits 0 and a typo exits 1.
enum class ParseResult { Ok, HelpRequested, Error };

ParseResult parseArgs(int argc, char **argv, Config &cfg) {
  for (int i = 1; i < argc; ++i) {
    const std::string arg = argv[i];

    auto next = [&](size_t &dest) {
      if (i + 1 >= argc) {
        std::cerr << "bench: " << arg << " requires a value\n";
        return false;
      }
      dest = static_cast<size_t>(std::strtoull(argv[++i], nullptr, 10));
      return true;
    };

    if (arg == "--help" || arg == "-h") {
      std::cerr << kUsage << '\n';
      return ParseResult::HelpRequested;
    } else if (arg == "--fasta") {
      if (i + 1 >= argc) {
        std::cerr << "bench: --fasta requires a value\n";
        return ParseResult::Error;
      }
      cfg.fasta = argv[++i];
    } else if (arg == "--max-size") {
      if (!next(cfg.maxSize)) return ParseResult::Error;
    } else if (arg == "--steps") {
      if (!next(cfg.steps)) return ParseResult::Error;
    } else if (arg == "--queries") {
      if (!next(cfg.queries)) return ParseResult::Error;
    } else if (arg == "--seed") {
      size_t s = 0;
      if (!next(s)) return ParseResult::Error;
      cfg.seed = static_cast<uint32_t>(s);
    } else {
      std::cerr << "bench: unknown argument: " << arg << '\n' << kUsage << '\n';
      return ParseResult::Error;
    }
  }

  if (cfg.steps < 1) {
    std::cerr << "bench: --steps must be at least 1\n";
    return ParseResult::Error;
  }
  return ParseResult::Ok;
}

} // namespace

int main(int argc, char **argv) {
  Config cfg;
  switch (parseArgs(argc, argv, cfg)) {
  case ParseResult::HelpRequested: return 0;
  case ParseResult::Error:         return 1;
  case ParseResult::Ok:            break;
  }

  std::mt19937 rng(cfg.seed);

  std::string text;
  if (!cfg.fasta.empty()) {
    if (!loadFasta(cfg.fasta, text))
      return 1;
    std::cerr << "bench: loaded " << text.size() << " bases from " << cfg.fasta
              << '\n';
    // --max-size truncates a real genome; it cannot extend one.
    if (cfg.maxSize > text.size()) {
      std::cerr << "bench: clamping --max-size to " << text.size() << '\n';
      cfg.maxSize = text.size();
    }
    text.resize(cfg.maxSize);
  } else {
    text = generateGenome(cfg.maxSize, rng);
    std::cerr << "bench: generated " << text.size() << " synthetic bases (seed "
              << cfg.seed << ")\n";
  }

  if (text.empty()) {
    std::cerr << "bench: empty input, nothing to measure\n";
    return 1;
  }

  std::cout << "experiment,structure,x,metric,value\n";

  runScaling(text, cfg, rng);
  runQuery(text, cfg, rng);

  std::cout.flush();
  return 0;
}
