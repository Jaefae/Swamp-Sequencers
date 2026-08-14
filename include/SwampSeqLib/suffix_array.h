#pragma once
#include "SwampSeqLib/genome_mapper.h" // GenomeMapper
#include <cstdint>
#include <string>
#include <vector>

// ─────────────────────────────────────────────────────────────────────────────
// SearchResult
//   Returned by SuffixArray::search().
//   offset – 0-based position of the match in the original text.
//   length – length of the matched pattern (equals pattern.size()).
// ─────────────────────────────────────────────────────────────────────────────

struct SearchResult {
  size_t offset;
  size_t length;
  bool operator==(const SearchResult &other) const;
};

// ─────────────────────────────────────────────────────────────────────────────
// SuffixArray

// Builds a suffix array over a read-only text using the SA-IS algorithm
// (Nong, Zhang & Chan, 2009) in O(n) time and O(n) auxiliary space.

// The text must remain alive for the lifetime of the object when constructed
// from a GenomeMapper (the object stores a raw pointer to the mapper's
// buffer). When constructed from a std::string the string is copied
// internally, so lifetime is managed automatically.

// Search runs in O(m log n) where m = pattern length and n = text length.
// ─────────────────────────────────────────────────────────────────────────────

class SuffixArray {
public:

  // ── Construction ──────────────────────────────────────────────────────────

  // Build from a GenomeMapper.  The mapper must be valid and must outlive this object.
  explicit SuffixArray(const GenomeMapper &mapper);

  // Build from an arbitrary string (text is copied internally).
  explicit SuffixArray(const std::string &text);

  // Non-copyable; moving is fine.
  SuffixArray(const SuffixArray &) = delete;
  SuffixArray &operator=(const SuffixArray &) = delete;
  SuffixArray(SuffixArray &&) = default;
  SuffixArray &operator=(SuffixArray &&) = default;

  // ── Query ─────────────────────────────────────────────────────────────────

  // Return all positions at which a pattern occurs in the text, sorted by
  // ascending offset. Returns an empty vector if the pattern is absent or
  // the array has not been built yet.
  std::vector<SearchResult>
  search(const std::string &pattern) const;

  // Number of characters in the indexed text.
  size_t size() const noexcept;

  // True once the suffix array has been successfully built.
  bool ready() const noexcept;

  // Direct read-only access to the underlying suffix array.
  const std::vector<size_t> &sa() const noexcept;

  // Heap bytes held by the index itself, excluding the indexed text (which
  // both structures share). This is the number the benchmark plots against
  // SuffixTree::memoryBytes(); for a suffix array it is simply the array.
  size_t memoryBytes() const noexcept;

private:
  // ── Internal data ─────────────────────────────────────────────────────────

  const char *_data = nullptr; // Pointer to the text (not owned when constructed from GenomeMapper).
  std::string _owned;  // Storage when constructed from string.
  size_t _num = 0; // Text length.
  std::vector<size_t> _sa; // The suffix array.
  bool _ready = false;

  // ── SA-IS implementation ──────────────────────────────────────────────────

  // Top-level driver: populate _sa using SA-IS.
  void buildSuffixArray();

  // Recursive SA-IS worker operating on an integer alphabet [0, alphabetSize).
  // s – input string as integers (sentinel = 0, must be unique minimum)
  // sa – output suffix array (same length as s)
  // alphabetSize – one past the maximum symbol value
  static void sa_is(const std::vector<uint32_t> &s, std::vector<size_t> &sa, uint32_t alphabetSize);

  // ── Binary-search helpers for pattern search ──────────────────────────────

  size_t lowerBound(const std::string &pattern) const;
  size_t upperBound(const std::string &pattern) const;
};