#include "SwampSeqLib/suffix_tree.h"
#include <algorithm>
#include <cstdint>
#include <stdexcept>

// ─────────────────────────────────────────────────────────────────────────────
// STSearchResult
// ─────────────────────────────────────────────────────────────────────────────

bool STSearchResult::operator==(const STSearchResult &other) const {
  return (this->offset == other.offset) && (this->length == other.length);
}

// ─────────────────────────────────────────────────────────────────────────────
// Construction
// ─────────────────────────────────────────────────────────────────────────────

SuffixTree::SuffixTree(const GenomeMapper &mapper) {
  // The mapper must be valid before we can safely read its underlying data.
  if (!mapper.isValid())
    throw std::runtime_error("SuffixTree: GenomeMapper is not valid.");

  // The suffix tree references external storage provided by GenomeMapper.
  // The mapper must outlive this SuffixTree instance.
  _data = mapper.data();
  _num = static_cast<int64_t>(mapper.size());

  // Build the tree immediately so the object is ready for search queries.
  buildSuffixTree();
  _ready = true;
}

SuffixTree::SuffixTree(const std::string &text) {
  // Keep an internal copy so the tree can safely reference stable storage.
  _owned = text;
  _data = _owned.c_str();
  _num = static_cast<int64_t>(_owned.size());

  // Build the tree immediately so the object is ready for search queries.
  buildSuffixTree();
  _ready = true;
}

// ─────────────────────────────────────────────────────────────────────────────
// Ukkonen's Algorithm — O(n) time, O(n) space

// Overview
// The tree is built left to right, one character at a time. Each new character either extends existing
// leaves for free, creates new ones where the text diverges, or stops early when the character already
// exists in the tree — carrying unresolved suffixes forward to the next character. An active point
// tracks where to resume, so nothing is rescanned, and suffix links between internal nodes keep each
// step O(1), giving the whole build O(n).
// ─────────────────────────────────────────────────────────────────────────────

// ── Node allocation
// ───────────────────────────────────────────────────────────

int64_t SuffixTree::newLeaf(int64_t start) {
  Node leaf;
  leaf.start = start;
  leaf.end = &_globalEnd; // open — shared with all leaves
  leaf.suffixIndex = -1;  // filled in by annotateSuffixIndices()
  leaf.suffixLink = NO_NODE;

  _nodes.push_back(std::move(leaf));
  return static_cast<int64_t>(_nodes.size()) - 1;
}

int64_t SuffixTree::newInternal(int64_t start, int64_t end) {
  // Internal nodes need their own private end value.
  _nodeEnds.push_back(end);
  int64_t *endPtr = &_nodeEnds.back();

  Node node;
  node.start = start;
  node.end = endPtr;
  node.suffixIndex = -1;
  node.suffixLink = ROOT; // default suffix link to root

  _nodes.push_back(std::move(node));
  return static_cast<int64_t>(_nodes.size()) - 1;
}

// ── Single-phase extension
// ────────────────────────────────────────────────────

void SuffixTree::extendTree(int64_t pos) {
  // globalEnd advances to pos+1, extending all open leaves.
  _globalEnd = pos + 1;
  ++_remaining;

  int64_t lastNewInternal = NO_NODE; //for suffix-link chaining

  while (_remaining > 0) {
    // Determine the character we need to insert at the active point.
    if (_activeLength == 0) {
      _activeEdge = pos; // edge key = current character index
    }

    // Read the character at the active edge from the text.
    // Use _num as the sentinel index (value = -1, never in real text).
    auto charAt = [&](int64_t idx) -> int64_t {
      if (idx == _num)
        return -1; // sentinel
      return static_cast<unsigned char>(_data[idx]);
    };

    const int64_t activeChar = charAt(_activeEdge);

    // Does the active node have a child edge starting with activeChar?
    auto &activeChildren = _nodes[_activeNode].children;
    auto childIt = activeChildren.find(activeChar);

    if (childIt == activeChildren.end()) {
      // ── no edge for this character, create a new leaf. ──────────
      int64_t leaf = newLeaf(pos);
      // Re-fetch reference: _nodes may have reallocated.
      _nodes[_activeNode].children[activeChar] = leaf;

      // Chain suffix link from the previously created internal node.
      if (lastNewInternal != NO_NODE) {
        _nodes[lastNewInternal].suffixLink = _activeNode;
        lastNewInternal = NO_NODE;
      }
    } else {
      // There is an existing child edge. Walk down if activeLength ≥ edge len.
      int64_t childIdx = childIt->second;

      // (skip/count): if activeLength spans an entire edge, descend.
      const int64_t edgeLen = _nodes[childIdx].edgeLength();
      if (_activeLength >= edgeLen) {
        _activeEdge += edgeLen;
        _activeLength -= edgeLen;
        _activeNode = childIdx;
        continue; // re-examine with updated active point
      }

      // Check if the next character on the edge already matches pos character.
      const int64_t nextOnEdge = _nodes[childIdx].start + _activeLength;
      if (charAt(nextOnEdge) == charAt(pos)) {
        // ── character already present...stop this phase. ─────────────
        ++_activeLength;
        if (lastNewInternal != NO_NODE) {
          _nodes[lastNewInternal].suffixLink = _activeNode;
        }
        break; // showstopper: remaining increments will be handled later
      }

      int64_t splitNode = newInternal(_nodes[childIdx].start, nextOnEdge);

      // The old child's edge now starts at nextOnEdge.
      _nodes[childIdx].start = nextOnEdge;

      // Attach an old child and new leaf to splitNode.
      // After newInternal/_nodes may have reallocated — use indices.
      _nodes[splitNode].children[charAt(nextOnEdge)] = childIdx;

      int64_t newLeafIdx = newLeaf(pos);
      _nodes[splitNode].children[charAt(pos)] = newLeafIdx;

      // Attach splitNode to activeNode, replacing old child.
      _nodes[_activeNode].children[activeChar] = splitNode;

      // Suffix-link the previous internal node to this one.
      if (lastNewInternal != NO_NODE) {
        _nodes[lastNewInternal].suffixLink = splitNode;
      }
      lastNewInternal = splitNode;
    }

    // One more suffix has been explicitly inserted.
    --_remaining;

    // Follow a suffix link or step toward root.
    if (_activeNode == ROOT && _activeLength > 0) {
      --_activeLength;
      _activeEdge = pos - _remaining + 1;
    } else if (_nodes[_activeNode].suffixLink != NO_NODE &&
               _nodes[_activeNode].suffixLink != NO_NODE) {
      _activeNode = _nodes[_activeNode].suffixLink;
    } else {
      _activeNode = ROOT;
    }
  }
}

// ── DFS suffix-index annotation
// After the tree is built, each leaf needs to know which suffix it represents.
// This is calculated by walking the tree and subtracting the root-to-leaf path length from the total text length.
// Done iteratively with an explicit stack rather than recursion to avoid stack overflow on deep genomic inputs.

void SuffixTree::annotateSuffixIndices(int64_t rootIdx, int64_t /*unused*/) {
  // Stack entries: (nodeIndex, accumulatedLabelHeight)
  struct Frame {
    int64_t nodeIdx;
    int64_t height;
  };

  std::vector<Frame> stack;
  stack.reserve(1024);
  stack.push_back({rootIdx, 0});

  while (!stack.empty()) {
    auto [nodeIdx, height] = stack.back();
    stack.pop_back();

    Node &node = _nodes[nodeIdx];

    if (node.children.empty()) {
      // Leaf — stamp suffix index.
      node.suffixIndex = (_num + 1) - height;
      continue;
    }

    // Push all children with their accumulated height.
    for (auto &[edgeChar, childIdx] : node.children) {
      const int64_t childHeight = height + _nodes[childIdx].edgeLength();
      stack.push_back({childIdx, childHeight});
    }
  }
}

// ── Top-level build driver
// ────────────────────────────────────────────────────

void SuffixTree::buildSuffixTree() {
  if (_num == 0)
    return;

  // Total characters processed = text + 1 sentinel.
  const int64_t total = _num + 1;

  // Pre-reserve to avoid repeated reallocation (Ukkonen creates at most 2n nodes).
  // _nodeEnds must also be reserved.
  _nodes.reserve(static_cast<size_t>(2 * total + 2));
  _nodeEnds.reserve(static_cast<size_t>(total + 2));

  // Create the root (node 0). Root has no edge label.
  // Give it a dummy private end so the end pointer is never null.
  _nodeEnds.push_back(0);
  Node root;
  root.start = -1;
  root.end = &_nodeEnds.back();
  root.suffixLink = NO_NODE;
  root.suffixIndex = -1;
  _nodes.push_back(std::move(root));

  _activeNode = ROOT;
  _activeEdge = -1;
  _activeLength = 0;
  _remaining = 0;
  _globalEnd = 0;

  // Process each real character then the sentinel.
  for (int64_t i = 0; i < total; ++i) {
    extendTree(i);
  }

  // Stamp suffix indices onto all leaves via DFS.
  annotateSuffixIndices(ROOT, 0);
}

// ─────────────────────────────────────────────────────────────────────────────
// Memory accounting
// ─────────────────────────────────────────────────────────────────────────────

size_t SuffixTree::memoryBytes() const noexcept {
  // Node pool and the private end values for internal nodes. capacity()
  // rather than size() because buildSuffixTree() reserves the 2n+2 worst case
  // up front, and that reservation is genuinely resident.
  size_t total = _nodes.capacity() * sizeof(Node) +
                 _nodeEnds.capacity() * sizeof(int64_t);

  // Each node carries its own child hash map, which allocates outside the
  // node itself: a bucket array plus one heap node per child edge. This is an
  // estimate — the exact layout is implementation defined — but it dominates
  // the footprint and cannot be ignored.
  for (const Node &node : _nodes) {
    total += node.children.bucket_count() * sizeof(void *);
    total += node.children.size() *
             (sizeof(std::pair<const int64_t, int64_t>) + sizeof(void *));
  }

  return total;
}

// ─────────────────────────────────────────────────────────────────────────────
// Search — O(m) descent, then DFS to collect all leaf offsets
// ─────────────────────────────────────────────────────────────────────────────

// Iterative DFS — avoids stack overflow on deep subtrees (same root cause
// as annotateSuffixIndices: genome-scale trees can be millions of nodes deep).
void SuffixTree::collectLeaves(int64_t subtreeRoot, std::vector<STSearchResult> &out, int64_t patternLength) const {
  std::vector<int64_t> stack;
  stack.reserve(1024); // increases the vector’s capacity to hold 1024 items without needing to grow immediately.

  stack.push_back(subtreeRoot);

  while (!stack.empty()) {
    const int64_t nodeIdx = stack.back();
    stack.pop_back();

    const Node &node = _nodes[nodeIdx];

    if (node.children.empty()) {
      // Leaf — record if it is a real suffix.
      if (node.suffixIndex >= 0 && node.suffixIndex < _num) {
        out.push_back({node.suffixIndex, patternLength});
      }
      continue;
    }

    for (const auto &[edgeChar, childIdx] : node.children) {
      stack.push_back(childIdx);
    }
  }
}

std::vector<STSearchResult>
SuffixTree::search(const std::string &pattern) const {
  if (!_ready || pattern.empty() || _nodes.empty())
    return {};

  const int64_t patternLength = static_cast<int64_t>(pattern.size());

  // ── Descend the tree matching the pattern ────────────────────────────────
  int64_t currentNode = ROOT;
  int64_t patPos = 0; // how many pattern characters matched so far

  while (patPos < patternLength) {
    const Node &node = _nodes[currentNode];

    // Look for a child edge whose first character matches pattern.
    const int64_t edgeChar =
        static_cast<unsigned char>(pattern[static_cast<size_t>(patPos)]);

    auto childIt = node.children.find(edgeChar);
    if (childIt == node.children.end()) {
      return {}; // no match
    }

    int64_t childIdx = childIt->second;
    const Node &childNode = _nodes[childIdx];

    // Walk along this edge, comparing pattern characters to text characters.
    const int64_t edgeStart = childNode.start;
    const int64_t edgeEnd = *childNode.end; // exclusive

    for (int64_t edgePos = edgeStart;
         edgePos < edgeEnd && patPos < patternLength; ++edgePos, ++patPos) {
      if (static_cast<unsigned char>(_data[edgePos]) !=
          static_cast<unsigned char>(pattern[static_cast<size_t>(patPos)])) {
        return {}; // mismatch
      }
    }

    // If we consumed the whole edge but still have a pattern left, descend.
    if (patPos < patternLength) {
      currentNode = childIdx;
    } else {
      // Pattern fully matched — collect all leaves in this subtree.
      std::vector<STSearchResult> results;
      collectLeaves(childIdx, results, patternLength);

      std::sort(results.begin(), results.end(),
                [](const STSearchResult &a, const STSearchResult &b) {
                  return a.offset < b.offset;
                });

      return results;
    }
  }

  // Pattern matched exactly at an internal node — collect its whole subtree.
  std::vector<STSearchResult> results;
  collectLeaves(currentNode, results, patternLength);

  std::sort(results.begin(), results.end(), [](const STSearchResult &a, const STSearchResult &b) {
              return a.offset < b.offset;
            });

  return results;
}
