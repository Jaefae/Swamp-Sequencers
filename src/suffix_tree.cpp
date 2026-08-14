#include "SwampSeqLib/suffix_tree.h"
#include <algorithm>
#include <cstdint>
#include <stdexcept>

bool STSearchResult::operator==(const STSearchResult &other) const {
  return (this->offset == other.offset) && (this->length == other.length);
}

SuffixTree::SuffixTree(const GenomeMapper &mapper) {
  if (!mapper.isValid())
    throw std::runtime_error("SuffixTree: GenomeMapper is not valid.");

  // References the mapper's storage, so the mapper must outlive this tree.
  _data = mapper.data();
  _num = static_cast<int64_t>(mapper.size());

  buildSuffixTree();
  _ready = true;
}

SuffixTree::SuffixTree(const std::string &text) {
  _owned = text;
  _data = _owned.c_str();
  _num = static_cast<int64_t>(_owned.size());

  buildSuffixTree();
  _ready = true;
}

// Ukkonen's algorithm, O(n) time and space.
//
// The tree is built left to right, one character at a time. Each new character
// either extends existing leaves for free, creates new ones where the text
// diverges, or stops early when the character is already present — carrying
// unresolved suffixes forward. An active point tracks where to resume so
// nothing is rescanned, and suffix links between internal nodes keep each step
// O(1).

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
  node.suffixLink = ROOT;

  _nodes.push_back(std::move(node));
  return static_cast<int64_t>(_nodes.size()) - 1;
}

void SuffixTree::extendTree(int64_t pos) {
  // Advancing globalEnd extends every open leaf at once.
  _globalEnd = pos + 1;
  ++_remaining;

  int64_t lastNewInternal = NO_NODE; // for suffix-link chaining

  while (_remaining > 0) {
    if (_activeLength == 0) {
      _activeEdge = pos;
    }

    // _num indexes the sentinel, which is never a real character.
    auto charAt = [&](int64_t idx) -> int64_t {
      if (idx == _num)
        return -1;
      return static_cast<unsigned char>(_data[idx]);
    };

    const int64_t activeChar = charAt(_activeEdge);

    auto &activeChildren = _nodes[_activeNode].children;
    auto childIt = activeChildren.find(activeChar);

    if (childIt == activeChildren.end()) {
      int64_t leaf = newLeaf(pos);
      // Re-index rather than reuse activeChildren: _nodes may have reallocated.
      _nodes[_activeNode].children[activeChar] = leaf;

      if (lastNewInternal != NO_NODE) {
        _nodes[lastNewInternal].suffixLink = _activeNode;
        lastNewInternal = NO_NODE;
      }
    } else {
      int64_t childIdx = childIt->second;

      // Skip/count: if activeLength spans a whole edge, descend past it.
      const int64_t edgeLen = _nodes[childIdx].edgeLength();
      if (_activeLength >= edgeLen) {
        _activeEdge += edgeLen;
        _activeLength -= edgeLen;
        _activeNode = childIdx;
        continue;
      }

      const int64_t nextOnEdge = _nodes[childIdx].start + _activeLength;
      if (charAt(nextOnEdge) == charAt(pos)) {
        // Already present — stop here; the remaining suffixes are carried to
        // the next phase.
        ++_activeLength;
        if (lastNewInternal != NO_NODE) {
          _nodes[lastNewInternal].suffixLink = _activeNode;
        }
        break;
      }

      int64_t splitNode = newInternal(_nodes[childIdx].start, nextOnEdge);

      _nodes[childIdx].start = nextOnEdge;

      _nodes[splitNode].children[charAt(nextOnEdge)] = childIdx;

      int64_t newLeafIdx = newLeaf(pos);
      _nodes[splitNode].children[charAt(pos)] = newLeafIdx;

      _nodes[_activeNode].children[activeChar] = splitNode;

      if (lastNewInternal != NO_NODE) {
        _nodes[lastNewInternal].suffixLink = splitNode;
      }
      lastNewInternal = splitNode;
    }

    --_remaining;

    // Follow a suffix link, or step back toward the root.
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

// Each leaf learns which suffix it represents by subtracting its root-to-leaf
// path length from the text length. Iterative rather than recursive: genomic
// inputs go deep enough to overflow the call stack.
void SuffixTree::annotateSuffixIndices(int64_t rootIdx, int64_t /*unused*/) {
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
      node.suffixIndex = (_num + 1) - height;
      continue;
    }

    for (auto &[edgeChar, childIdx] : node.children) {
      const int64_t childHeight = height + _nodes[childIdx].edgeLength();
      stack.push_back({childIdx, childHeight});
    }
  }
}

void SuffixTree::buildSuffixTree() {
  if (_num == 0)
    return;

  const int64_t total = _num + 1; // text plus sentinel

  // Ukkonen creates at most 2n nodes; reserving up front avoids reallocation.
  _nodes.reserve(static_cast<size_t>(2 * total + 2));
  _nodeEnds.reserve(static_cast<size_t>(total + 2));

  // Root has no edge label, but gets a dummy end so the pointer is never null.
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

  for (int64_t i = 0; i < total; ++i) {
    extendTree(i);
  }

  annotateSuffixIndices(ROOT, 0);
}

size_t SuffixTree::memoryBytes() const noexcept {
  // capacity() rather than size(): buildSuffixTree() reserves the 2n+2 worst
  // case up front and that reservation is genuinely resident.
  size_t total = _nodes.capacity() * sizeof(Node) +
                 _nodeEnds.capacity() * sizeof(int64_t);

  // Each node's child map allocates outside the node: a bucket array plus one
  // heap node per edge. Estimated, since the layout is implementation defined,
  // but it dominates the footprint.
  for (const Node &node : _nodes) {
    total += node.children.bucket_count() * sizeof(void *);
    total += node.children.size() *
             (sizeof(std::pair<const int64_t, int64_t>) + sizeof(void *));
  }

  return total;
}

// Iterative for the same reason as annotateSuffixIndices: genome-scale
// subtrees are deep enough to overflow the call stack.
void SuffixTree::collectLeaves(int64_t subtreeRoot, std::vector<STSearchResult> &out, int64_t patternLength) const {
  std::vector<int64_t> stack;
  stack.reserve(1024);

  stack.push_back(subtreeRoot);

  while (!stack.empty()) {
    const int64_t nodeIdx = stack.back();
    stack.pop_back();

    const Node &node = _nodes[nodeIdx];

    if (node.children.empty()) {
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

  int64_t currentNode = ROOT;
  int64_t patPos = 0;

  while (patPos < patternLength) {
    const Node &node = _nodes[currentNode];

    const int64_t edgeChar =
        static_cast<unsigned char>(pattern[static_cast<size_t>(patPos)]);

    auto childIt = node.children.find(edgeChar);
    if (childIt == node.children.end()) {
      return {};
    }

    int64_t childIdx = childIt->second;
    const Node &childNode = _nodes[childIdx];

    const int64_t edgeStart = childNode.start;
    const int64_t edgeEnd = *childNode.end; // exclusive

    for (int64_t edgePos = edgeStart;
         edgePos < edgeEnd && patPos < patternLength; ++edgePos, ++patPos) {
      if (static_cast<unsigned char>(_data[edgePos]) !=
          static_cast<unsigned char>(pattern[static_cast<size_t>(patPos)])) {
        return {};
      }
    }

    if (patPos < patternLength) {
      currentNode = childIdx;
    } else {
      std::vector<STSearchResult> results;
      collectLeaves(childIdx, results, patternLength);

      std::sort(results.begin(), results.end(),
                [](const STSearchResult &a, const STSearchResult &b) {
                  return a.offset < b.offset;
                });

      return results;
    }
  }

  // Pattern ended exactly at an internal node — collect its whole subtree.
  std::vector<STSearchResult> results;
  collectLeaves(currentNode, results, patternLength);

  std::sort(results.begin(), results.end(), [](const STSearchResult &a, const STSearchResult &b) {
              return a.offset < b.offset;
            });

  return results;
}
