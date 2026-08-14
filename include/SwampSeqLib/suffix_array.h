#pragma once
#include "SwampSeqLib/genome_mapper.h"
#include <cstdint>
#include <string>
#include <vector>

// offset is a 0-based position in the original text; length equals pattern.size().
struct SearchResult {
  size_t offset;
  size_t length;
  bool operator==(const SearchResult &other) const;
};

// Suffix array built with SA-IS (Nong, Zhang & Chan, 2009) in O(n) time and
// O(n) auxiliary space. Search runs in O(m log n).
class SuffixArray {
public:
  // The mapper must be valid and must outlive this object; only a raw pointer
  // to its buffer is stored.
  explicit SuffixArray(const GenomeMapper &mapper);

  // Copies the text, so lifetime is managed automatically.
  explicit SuffixArray(const std::string &text);

  SuffixArray(const SuffixArray &) = delete;
  SuffixArray &operator=(const SuffixArray &) = delete;
  SuffixArray(SuffixArray &&) = default;
  SuffixArray &operator=(SuffixArray &&) = default;

  // Every position the pattern occurs at, ascending. Empty if the pattern is
  // absent or the array has not been built.
  std::vector<SearchResult>
  search(const std::string &pattern) const;

  size_t size() const noexcept;
  bool ready() const noexcept;
  const std::vector<size_t> &sa() const noexcept;

  // Heap bytes held by the index, excluding the text both structures share.
  // The benchmark plots this against SuffixTree::memoryBytes().
  size_t memoryBytes() const noexcept;

private:
  const char *_data = nullptr; // not owned when built from a GenomeMapper
  std::string _owned;
  size_t _num = 0;
  std::vector<size_t> _sa;
  bool _ready = false;

  void buildSuffixArray();

  // Recursive worker over an integer alphabet [0, alphabetSize), where the
  // sentinel 0 must be a unique minimum.
  static void sa_is(const std::vector<uint32_t> &s, std::vector<size_t> &sa, uint32_t alphabetSize);

  size_t lowerBound(const std::string &pattern) const;
  size_t upperBound(const std::string &pattern) const;
};
