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

// Long enough that a sampled pattern is almost always unique, so the scaling
// numbers measure the search rather than the cost of reporting every hit.
constexpr size_t kScalingPatternLength = 20;

constexpr size_t kQueryLengths[] = {4, 8, 12, 16, 24, 32, 48, 64};

struct Config {
  std::string fasta;
  size_t maxSize = 1000000;
  size_t steps = 5;
  size_t queries = 1000;
  uint32_t seed = 20240817;
};

void emit(const char *experiment, const char *structure, size_t x,
          const char *metric, double value) {
  std::cout << experiment << ',' << structure << ',' << x << ',' << metric << ','
            << std::fixed << std::setprecision(4) << value << '\n';
}

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

std::string generateGenome(size_t n, std::mt19937 &rng) {
  static const char kBases[] = {'A', 'C', 'G', 'T'};
  std::uniform_int_distribution<int> pick(0, 3);
  std::string text;
  text.reserve(n);
  for (size_t i = 0; i < n; ++i)
    text.push_back(kBases[pick(rng)]);
  return text;
}

// Sampling from the text is what makes total_hits a usable correctness check.
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

template <typename Index>
QueryStats timeQueries(const Index &index, const std::vector<std::string> &pats) {
  QueryStats stats;
  if (pats.empty())
    return stats;

  const auto start = Clock::now();
  // Accumulating inside the timed loop consumes the result, so the optimizer
  // cannot discard the search as dead code.
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

void runScaling(const std::string &text, const Config &cfg, std::mt19937 &rng) {
  for (size_t step = 1; step <= cfg.steps; ++step) {
    const size_t n = text.size() * step / cfg.steps;
    if (n < kScalingPatternLength)
      continue;

    std::cerr << "bench: scaling n=" << n << " ...";
    const std::string sub = text.substr(0, n);
    const std::vector<std::string> pats =
        samplePatterns(sub, kScalingPatternLength, cfg.queries, rng);

    // Both classes hold a raw pointer into their own copy of the text, so they
    // are built in place and never moved.
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

const char *kUsage =
    "usage: SwampBench [--fasta PATH] [--max-size N] [--steps N] "
    "[--queries N] [--seed N]";

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
