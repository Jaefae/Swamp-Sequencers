#pragma once
#include "SwampSeqLib/genome_mapper.h"
#include <cstdint>
#include <string>
#include <unordered_map>
#include <vector>

// offset is a 0-based position in the original text; length equals pattern.size().
struct STSearchResult {
  int64_t offset;
  int64_t length;
  bool operator==(const STSearchResult &other) const;
};

// Compressed trie of all suffixes, built with Ukkonen's online algorithm in
// O(n) time and O(n) space. Search runs in O(m).
class SuffixTree {
public:
  // The mapper must outlive this object; only a raw pointer to it is stored.
  explicit SuffixTree(const GenomeMapper &mapper);

  // Copies the text, so lifetime is managed automatically.
  explicit SuffixTree(const std::string &text);

  SuffixTree(const SuffixTree &) = delete;
  SuffixTree &operator=(const SuffixTree &) = delete;
  SuffixTree(SuffixTree &&) = default;
  SuffixTree &operator=(SuffixTree &&) = default;

  ~SuffixTree() = default;

  // Every position the pattern occurs at, ascending. Empty if the pattern is
  // absent or the tree has not been built.
  std::vector<STSearchResult> search(const std::string &pattern) const;

  int64_t size() const noexcept { return _num; }
  bool ready() const noexcept { return _ready; }

  // Heap bytes held by the index, excluding the text both structures share:
  // node pool, per-internal-node end values, and every node's child map. The
  // maps are estimated from bucket and element counts, so this is close to but
  // not exactly what the allocator reports.
  size_t memoryBytes() const noexcept;

  // Ukkonen allocates at most 2n+1 nodes.
  size_t nodeCount() const noexcept { return _nodes.size(); }

private:
  struct Node {
    int64_t  start = -1;
    int64_t *end  = nullptr; // points into _nodeEnds or _globalEnd
    int64_t suffixLink = -1;
    int64_t suffixIndex= -1; // -1 for internal nodes
    std::unordered_map<int64_t, int64_t> children;

    // Leaves use *end, which equals _globalEnd at query time.
    int64_t edgeLength() const { return *end - start; }
  };

  const char  *_data = nullptr; // not owned when built from a GenomeMapper
  std::string _owned;
  int64_t _num = 0; // excludes the sentinel
  bool _ready = false;

  std::vector<Node> _nodes; // pool; avoids pointer invalidation on reallocation
  std::vector<int64_t> _nodeEnds; // per-internal-node; leaves share _globalEnd
  int64_t _globalEnd = 0; // shared open end for all leaves, bumped each phase

  // Ukkonen active point.
  int64_t _activeNode = 0; // index into _nodes
  int64_t _activeEdge = -1; // first character of the active edge, as an index
  int64_t _activeLength = 0;
  int64_t _remaining = 0; // suffixes yet to be inserted

  static constexpr int64_t ROOT = 0;
  static constexpr int64_t NO_NODE = -1;

  void buildSuffixTree();
  int64_t newLeaf(int64_t start);
  int64_t newInternal(int64_t start, int64_t end);
  void extendTree(int64_t pos);

  // DFS stamping suffix indices onto every leaf.
  void annotateSuffixIndices(int64_t nodeIdx, int64_t labelHeight);

  void collectLeaves(int64_t nodeIdx, std::vector<STSearchResult> &out, int64_t patternLength) const;
};
