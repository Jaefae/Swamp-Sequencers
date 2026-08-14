#include "SwampSeqLib/suffix_array.h"

#include <algorithm>
#include <cstdint>
#include <iostream>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <vector>

// SearchResult comparison for testing
bool SearchResult::operator==(const SearchResult& other) const {
  return ((this->length == other.length) && (this->offset == other.offset));
}

// Construction

SuffixArray::SuffixArray(const GenomeMapper& mapper) {
  if (!mapper.isValid())
    throw std::runtime_error("SuffixArray: GenomeMapper is not valid.");

  _data = mapper.data();
  _num  = mapper.size();

  buildSuffixArray();
  _ready = true;
}

SuffixArray::SuffixArray(const std::string& text) {
  _owned = text;
  _data  = _owned.c_str();
  _num   = _owned.size();

  buildSuffixArray();
  _ready = true;
}

// SA-IS Helpers

// Helpers (file-scope)

// Return bucket sizes: buckets[c] = number of occurrences of symbol c.
std::vector<size_t> getBuckets(const std::vector<uint32_t>& symbols,
                               uint32_t                     alphabetSize) {
  // Safety check: Ensure bucket size covers the requested alphabet size
  std::vector<size_t> buckets(alphabetSize, 0);

  for (uint32_t symbol : symbols) {
    buckets[static_cast<size_t>(symbol)]++;
  }

  return buckets;
}

// Fill heads with the start index of each bucket (for L-type induced sort).
void getBucketHeads(const std::vector<size_t>& buckets,
                    std::vector<size_t>&       heads) {
  size_t offset = 0;
  for (size_t i = 0; i < buckets.size(); ++i) {
    heads[i] = offset;
    offset += buckets[i];
  }
}

// Fill tails with the last index (inclusive) of each bucket (for S-type).
void getBucketTails(const std::vector<size_t>& buckets,
                    std::vector<size_t>&       tails) {
  size_t offset = 0;
  for (size_t i = 0; i < buckets.size(); ++i) {
    offset += buckets[i];
    tails[i] = offset - 1;
  }
}

// Induced sorting pass

void induceSortL(const std::vector<uint32_t>& s, std::vector<size_t>& sa,
                 const std::vector<bool>&   isSType,
                 const std::vector<size_t>& buckets, uint32_t alphabetSize) {
  // Text length.
  const size_t textSize = s.size();

  // Bucket head positions: the next free slot at the start of each bucket.
  std::vector<size_t> heads(alphabetSize);
  getBucketHeads(buckets, heads);

  // Scan the suffix array from left to right.
  // Every already-placed suffix may reveal a predecessor that should be
  // inserted.
  for (size_t i = 0; i < textSize; ++i) {
    // Skip empty slots and the sentinel suffix at position 0.
    if (sa[i] == SIZE_MAX || sa[i] == 0) continue;

    // Consider the suffix that starts one character before the current suffix.
    const size_t predecessor = sa[i] - 1;

    // Only L-type suffixes are induced in this pass.
    if (!isSType[predecessor]) {
      // Place the predecessor at the head of the bucket for its first symbol,
      // then advance that bucket head.
      sa[heads[s[predecessor]]++] = predecessor;
    }
  }
}

void induceSortS(const std::vector<uint32_t>& s, std::vector<size_t>& sa,
                 const std::vector<bool>&   isSType,
                 const std::vector<size_t>& buckets, uint32_t alphabetSize) {
  // Length of the input text.
  const size_t n = s.size();

  // tails[c] will track the current last free position in the bucket for symbol
  // c.
  std::vector<size_t> tails(alphabetSize);
  getBucketTails(buckets, tails);

  // Scan the suffix array from right to left.
  // This lets us induce S-type suffixes while preserving the needed order.
  for (size_t i = n; i-- > 0;) {
    // Skip empty slots and the suffix starting at position 0
    // because it has no predecessor.
    if (sa[i] == SIZE_MAX || sa[i] == 0) continue;

    // The predecessor suffix starts one character earlier.
    const size_t predecessor = sa[i] - 1;

    // Only S-type suffixes are inserted in this pass.
    if (isSType[predecessor]) {
      // Place the predecessor at the tail of its bucket,
      // then move the tail one position left.
      sa[tails[s[predecessor]]--] = predecessor;
    }
  }
}

void SuffixArray::sa_is(const std::vector<uint32_t>& s, std::vector<size_t>& sa,
                        uint32_t alphabetSize) {
  // Special marker for empty slots in the suffix array.
  const size_t kUnusedSlot = SIZE_MAX;

  // Sentinel value for ranks that have not been assigned yet.
  const uint32_t kEmptyRank = std::numeric_limits<uint32_t>::max();

  // Length of the input string, including the sentinel.
  const size_t n = s.size();

  // Step 1: Classify every position as S-type or L-type.
  // S-type means the suffix starting here is lexicographically
  // smaller than the one to its right; L-type means larger.
  std::vector<bool> isSType(n, false);
  isSType[n - 1] = true;  // The last position is always S-type.
  for (size_t i = n - 1; i-- > 0;) {
    isSType[i] = (s[i] < s[i + 1]) || (s[i] == s[i + 1] && isSType[i + 1]);
  }

  // LMS = Leftmost S-type position.
  // A position is LMS if it is S-type and the previous position is L-type.
  auto isLMSPosition = [&](size_t i) -> bool {
    return i > 0 && isSType[i] && !isSType[i - 1];
  };

  // Bucket sizes for each alphabet symbol.
  const std::vector<size_t> buckets = getBuckets(s, alphabetSize);

  // Step 2: Place LMS positions into the ends of their buckets.
  // The ordering used here is important for induced sorting.
  auto placeLMSIntoBuckets =
      [&](const std::vector<size_t>& orderedLMSPositions) {
        sa.assign(n, kUnusedSlot);

        // Tail positions for each bucket.
        std::vector<size_t> tails(alphabetSize);
        getBucketTails(buckets, tails);

        // Insert LMS positions from right to left for stability.
        for (size_t i = orderedLMSPositions.size(); i-- > 0;) {
          const size_t pos    = orderedLMSPositions[i];
          sa[tails[s[pos]]--] = pos;
        }

        // Put the sentinel suffix first.
        sa[0] = n - 1;
      };

  // Collect all LMS positions in text order.
  std::vector<size_t> initialLMSPositions;
  for (size_t i = 0; i < n; ++i) {
    if (isLMSPosition(i)) initialLMSPositions.push_back(i);
  }

  // Place LMS suffixes, then induce the remaining L- and S-type suffixes.
  placeLMSIntoBuckets(initialLMSPositions);
  induceSortL(s, sa, isSType, buckets, alphabetSize);
  induceSortS(s, sa, isSType, buckets, alphabetSize);

  // Step 3: Extract LMS positions in sorted order.
  std::vector<size_t> sortedLMSPositions;
  for (size_t i = 0; i < n; ++i) {
    if (isLMSPosition(sa[i])) sortedLMSPositions.push_back(sa[i]);
  }

  // If there are no LMS positions, the job is done.
  if (sortedLMSPositions.empty()) return;

  // Step 4: Assign ranks to LMS substrings.
  // Equal LMS substrings get the same rank.
  std::vector<uint32_t> rank(n, kEmptyRank);
  uint32_t              currentRank = 0;
  rank[sortedLMSPositions[0]]       = 0;

  for (size_t i = 1; i < sortedLMSPositions.size(); ++i) {
    const size_t previous = sortedLMSPositions[i - 1];
    const size_t current  = sortedLMSPositions[i];
    bool         differ   = false;

    // Compare LMS substrings character by character until they differ
    // or until both reach the next LMS boundary.
    for (size_t d = 0;; ++d) {
      const size_t pIdx = previous + d;
      const size_t cIdx = current + d;

      const bool pLMS = (d > 0) && isLMSPosition(pIdx);
      const bool cLMS = (d > 0) && isLMSPosition(cIdx);

      if (s[pIdx] != s[cIdx] || isSType[pIdx] != isSType[cIdx] ||
          pLMS != cLMS) {
        differ = true;
        break;
      }
      if (d > 0 && (pLMS || cLMS)) break;
    }

    if (differ) ++currentRank;
    rank[current] = currentRank;
  }

  // Step 5: Build the reduced string from LMS ranks.
  // Each LMS substring becomes one symbol in the smaller problem.
  std::vector<size_t>   lmsPositionsInTextOrder;
  std::vector<uint32_t> reducedString;
  for (size_t i = 0; i < n; ++i) {
    if (rank[i] != kEmptyRank) {
      lmsPositionsInTextOrder.push_back(i);
      reducedString.push_back(rank[i]);
    }
  }

  const uint32_t      reducedAlphabetSize = currentRank + 1;
  std::vector<size_t> saReduced(reducedString.size());

  // Step 6: Recursively solve the reduced problem if needed.
  // If ranks are already unique enough, build the reduced SA directly.
  if (reducedAlphabetSize < static_cast<uint32_t>(reducedString.size())) {
    sa_is(reducedString, saReduced, reducedAlphabetSize);
  } else {
    for (size_t i = 0; i < reducedString.size(); ++i) {
      saReduced[reducedString[i]] = i;
    }
  }

  // Map reduced suffix-array order back to original LMS positions.
  std::vector<size_t> orderedLMSPositions(saReduced.size());
  for (size_t i = 0; i < saReduced.size(); ++i) {
    orderedLMSPositions[i] = lmsPositionsInTextOrder[saReduced[i]];
  }

  // Insert the sorted LMS positions back into their bucket tails,
  // then complete the suffix array by inducing L-type and S-type suffixes.
  placeLMSIntoBuckets(orderedLMSPositions);
  induceSortL(s, sa, isSType, buckets, alphabetSize);
  induceSortS(s, sa, isSType, buckets, alphabetSize);
}

// SuffixArray Core Logic

// Build the suffix array for the currently stored text.
void SuffixArray::buildSuffixArray() {
  // If the input text is empty, there are no suffixes to build.
  if (_num == 0) {
    _sa.clear();
    return;
  }

  // Assume an 8-bit alphabet.
  constexpr uint32_t alphabetSize = 256;

  // Add one extra slot for the sentinel character.
  const size_t          totalSize = _num + 1;
  std::vector<uint32_t> textWithSentinel(totalSize);

  // Copy each input character into the integer-based working buffer.
  for (size_t i = 0; i < _num; ++i) {
    textWithSentinel[i] = static_cast<unsigned char>(_data[i]);
  }

  // Append a sentinel symbol that is smaller than all real characters.
  // This guarantees that the empty suffix is sorted first.
  textWithSentinel[_num] = 0;  // Sentinel

  // Allocate space for the full suffix array, including the sentinel suffix.
  std::vector<size_t> suffixArray(totalSize);

  // Run the SA-IS suffix array construction algorithm.
  sa_is(textWithSentinel, suffixArray, alphabetSize);

  // Drop the sentinel suffix from the final result.
  // The stored suffix array contains only suffixes of the original text.
  _sa.assign(suffixArray.begin() + 1, suffixArray.end());
}

std::vector<SearchResult> SuffixArray::search(
    const std::string& pattern) const {
  // If the suffix array is not built yet or the pattern is empty, there can be
  // no matches.
  if (!_ready || pattern.empty()) return {};

  // Find the first suffix that could match the pattern.
  const size_t firstMatch = lowerBound(pattern);

  // Find one past the last suffix that could match the pattern.
  const size_t pastLastMatch = upperBound(pattern);

  // If the range is empty, the pattern does not occur in the text.
  if (firstMatch >= pastLastMatch) return {};

  // Collect all matching suffix positions.
  std::vector<SearchResult> results;
  results.reserve(pastLastMatch - firstMatch);

  for (size_t index = firstMatch; index < pastLastMatch; ++index) {
    // Each match starts at the suffix-array position _sa[index],
    // and its length is the full pattern length.
    results.emplace_back(SearchResult{_sa[index], (uint32_t)pattern.size()});
  }

  // Sort matches by their offset in the original text for easier consumption.
  std::sort(results.begin(), results.end(),
            [](const SearchResult& lhs, const SearchResult& rhs) {
              return lhs.offset < rhs.offset;
            });

  return results;
}

size_t SuffixArray::size() const noexcept { return _num; }

bool SuffixArray::ready() const noexcept { return _ready; }

const std::vector<size_t>& SuffixArray::sa() const noexcept { return _sa; }

size_t SuffixArray::memoryBytes() const noexcept {
  // The whole index is one flat array of suffix offsets: 8 bytes per
  // character on a 64-bit build, with no per-element allocator overhead.
  // capacity() rather than size() because that is what is actually held.
  return _sa.capacity() * sizeof(size_t);
}

size_t SuffixArray::lowerBound(const std::string& pattern) const {
  // Length of the search pattern.
  const size_t patternLength = pattern.size();

  // Binary search range over the suffix array: [low, high).
  size_t low = 0, high = _num;

  while (low < high) {
    // Midpoint of the current search range.
    const size_t mid = low + (high - low) / 2;

    // Starting position
    const size_t suffixStart = _sa[mid];

    // Length of the suffix remaining
    const size_t suffixLength = _num - suffixStart;

    // Only compare as many characters as both strings have available.
    const size_t compareLength = std::min(patternLength, suffixLength);

    // Lexicographically compare the pattern with the current suffix prefix.
    const int cmp =
        pattern.compare(0, patternLength, _data + suffixStart, compareLength);

    // If suffix >= pattern, the lower bound is at mid or to the left of it.
    if (cmp <= 0) high = mid;
    // Otherwise, the lower bound must be to the right of mid.
    else
      low = mid + 1;
  }

  // low is now the first suffix that is not less than pattern.
  return low;
}

size_t SuffixArray::upperBound(const std::string& pattern) const {
  // Cache the pattern length once so we don't recompute it in the loop.
  const size_t patternLength = pattern.size();

  // Search range in the suffix array.
  // We are looking for the first suffix that is strictly greater than pattern.
  size_t low = 0, high = _num;

  while (low < high) {
    // Midpoint of the current search range.
    const size_t mid = low + (high - low) / 2;

    // Starting position of the suffix.
    const size_t suffixPos = _sa[mid];

    // Number of characters remaining.
    const size_t suffixLength = _num - suffixPos;

    // Compare only the portion that exists in both strings.
    const size_t compareLength = std::min(patternLength, suffixLength);

    // Lexicographically compare the full pattern against the suffix prefix.
    const int cmp =
        pattern.compare(0, patternLength, _data + suffixPos, compareLength);

    // For upperBound, move left only when the suffix is greater than the
    // pattern. Otherwise, move right to find the first strictly greater suffix.
    if (cmp < 0)
      high = mid;
    else
      low = mid + 1;
  }

  // low is now the index of the first suffix greater than pattern.
  return low;
}
