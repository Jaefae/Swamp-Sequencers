#include "../include/SwampSeqLib/suffix_array.h"
#include "SwampSeqLib/genome_mapper.h"
#include "SwampSeqLib/suffix_tree.h"
#include <algorithm>
#include <cstddef>
#include <cstdio>
#include <fstream>
#include <gtest/gtest.h>
#include <random>

void sortResults(std::vector<STSearchResult> &res) {
  std::sort(res.begin(), res.end(),
            [](const STSearchResult &a, const STSearchResult &b) {
              return a.offset < b.offset;
            });
}

const std::vector<char> possibleChars{'A', 'C', 'T', 'G', 'N'};

TEST(AlwaysTrue, AlwaysTrue) { ASSERT_EQ(0, 0); }

TEST(SuffixArray, BasicString) {
  const std::string text = "banana";
  SuffixArray SA(text);

  std::vector<size_t> expected{5, 3, 1, 0, 4, 2};
  EXPECT_EQ(SA.sa(), expected);
}

TEST(SuffixArray, GenomeString) {
  const std::string genome = "ACGGTCGATCACGAT";
  SuffixArray SA(genome);

  std::vector<SearchResult> expected{{0, 3}, {10, 3}};
  EXPECT_EQ(SA.search("ACG"), expected);
}

TEST(SuffixArray, Repetition) {
  std::string genome = "";
  for (size_t i = 0; i < 1000; i++)
    genome += 'A';
  SuffixArray SA(genome);

  std::vector<SearchResult> expected;
  for (size_t i = 0; i < 999; i += 1)
    expected.push_back({i, 2});
  auto res = SA.search("AA");
  for (size_t i = 0; i < 500; i++)
    EXPECT_EQ(res[i].offset, expected[i].offset);
}

TEST(SuffixArray, Overlapping) {
  std::string text = "ba";
  for (size_t i = 0; i < 100; i++)
    text += "na";

  SuffixArray SA(text);

  std::vector<SearchResult> expected;
  for (size_t i = 0; i < (text.size() - 1) / 3; i++) {
    size_t offset = i * 2 + 1;
    expected.push_back({offset, 3});
  }

  auto res = SA.search("ana");
  for (size_t i = 0; i < (text.size() - 1) / 3; i++) {
    EXPECT_EQ(res[i].offset, expected[i].offset);
  }
}

TEST(SuffixArray, NoMatch) {
  std::string genome = "";
  for (size_t i = 0; i < 1000; i++)
    genome += "A";

  SuffixArray SA(genome);
  auto res = SA.search("NO MATCH");

  EXPECT_EQ(res.size(), 0);
}

TEST(SuffixArray, Boundaries) {
  const std::string text = "GATTACA";
  SuffixArray SA(text);

  auto resStart = SA.search("GAT");
  ASSERT_EQ(resStart.size(), 1);
  EXPECT_EQ(resStart[0].offset, 0);

  auto resEnd = SA.search("ACA");
  ASSERT_EQ(resEnd.size(), 1);
  EXPECT_EQ(resEnd[0].offset, 4);
}

TEST(SuffixArray, AllMatch) {
  std::string genome(100, 'A');
  SuffixArray SA(genome);

  auto res = SA.search("A");
  EXPECT_EQ(res.size(), 100);
}

// Constructor shouldn't segfault on empty input
TEST(SuffixArray, EmptyInput) {
  EXPECT_NO_THROW({
    SuffixArray SA("");
    auto res = SA.search("A");
    EXPECT_EQ(res.size(), 0);
  });
}

TEST(SuffixArray, SingleChar) {
  SuffixArray SA("G");
  EXPECT_EQ(SA.search("G").size(), 1);
  EXPECT_EQ(SA.search("C").size(), 0);
}

TEST(SuffixTree, BasicString) {
  const std::string text = "banana";
  SuffixTree ST(text);

  // Searching for "a" to verify multiple leaf discovery
  std::vector<STSearchResult> res = ST.search("a");
  sortResults(res);

  ASSERT_EQ(res.size(), 3);
  EXPECT_EQ(res[0].offset, 1);
  EXPECT_EQ(res[1].offset, 3);
  EXPECT_EQ(res[2].offset, 5);
}

TEST(SuffixTree, GenomeString) {
  const std::string genome = "ACGGTCGATCACGAT";
  SuffixTree ST(genome);

  std::vector<STSearchResult> expected{{0, 3}, {10, 3}};
  auto res = ST.search("ACG");
  sortResults(res);

  ASSERT_EQ(res.size(), expected.size());
  for (size_t i = 0; i < res.size(); ++i) {
    EXPECT_EQ(res[i].offset, expected[i].offset);
    EXPECT_EQ(res[i].length, expected[i].length);
  }
}

TEST(SuffixTree, Repetition) {
  std::string genome(1000, 'A');
  SuffixTree ST(genome);

  auto res = ST.search("AA");
  sortResults(res);

  // "AA" appears at every index from 0 to 998
  ASSERT_EQ(res.size(), 999);
  for (int64_t i = 0; i < 999; i++) {
    EXPECT_EQ(res[i].offset, i);
    EXPECT_EQ(res[i].length, 2);
  }
}

TEST(SuffixTree, Overlapping) {
  std::string text = "ba";
  for (size_t i = 0; i < 100; i++)
    text += "na";

  SuffixTree ST(text);

  // Pattern "ana" occurs starting at indices 1, 3, 5...
  auto res = ST.search("ana");
  sortResults(res);

  size_t expected_count = (text.size() - 1) / 2; // "ana" in "bananana..."

  for (size_t i = 0; i < res.size(); i++) {
    EXPECT_EQ(res[i].offset, (int64_t)(i * 2 + 1));
    EXPECT_EQ(res[i].length, 3);
  }
}

TEST(SuffixTree, NoMatch) {
  std::string genome(1000, 'A');
  SuffixTree ST(genome);

  auto res = ST.search("NO MATCH");
  EXPECT_EQ(res.size(), 0);
}

TEST(SuffixTree, Boundaries) {
  const std::string text = "GATTACA";
  SuffixTree ST(text);

  auto resStart = ST.search("GAT");
  ASSERT_EQ(resStart.size(), 1);
  EXPECT_EQ(resStart[0].offset, 0);

  auto resEnd = ST.search("ACA");
  ASSERT_EQ(resEnd.size(), 1);
  EXPECT_EQ(resEnd[0].offset, 4);
}

TEST(SuffixTree, AllMatch) {
  std::string genome(100, 'A');
  SuffixTree ST(genome);

  auto res = ST.search("A");
  EXPECT_EQ(res.size(), 100);
}

TEST(SuffixTree, EmptyInput) {
  EXPECT_NO_THROW({
    SuffixTree ST("");
    auto res = ST.search("A");
    EXPECT_EQ(res.size(), 0);
    EXPECT_EQ(ST.size(), 0);
  });
}

TEST(SuffixTree, SingleChar) {
  SuffixTree ST("G");
  auto resG = ST.search("G");
  EXPECT_EQ(resG.size(), 1);
  if (!resG.empty())
    EXPECT_EQ(resG[0].offset, 0);

  EXPECT_EQ(ST.search("C").size(), 0);
}

TEST(SuffixTree, InternalState) {
  SuffixTree ST("banana");
  EXPECT_TRUE(ST.ready());
  EXPECT_EQ(ST.size(), 6);
}

// Differential tests
//
// The two structures are independent implementations of the same contract, so
// any query where they disagree is a bug in one of them. This catches the
// class of defect the per-structure tests above miss: both are checked against
// hand-written expectations on tiny inputs, while the bugs in this code have
// historically been out-of-bounds indexing that only appears at scale.

// Build a random genome over the standard nucleotide alphabet. Seeded
// explicitly so a CI failure reproduces locally from the test name alone.
static std::string randomGenome(size_t length, uint32_t seed) {
  std::mt19937 rng(seed);
  std::uniform_int_distribution<size_t> pick(0, possibleChars.size() - 1);

  std::string genome;
  genome.reserve(length);
  for (size_t i = 0; i < length; ++i)
    genome.push_back(possibleChars[pick(rng)]);

  return genome;
}

// Compare the two indexes over the same text for one pattern.
static void expectSameMatches(const SuffixArray &SA, const SuffixTree &ST,
                              const std::string &pattern) {
  auto arrayRes = SA.search(pattern);
  auto treeRes = ST.search(pattern);
  sortResults(treeRes);

  ASSERT_EQ(arrayRes.size(), treeRes.size())
      << "match count differs for pattern \"" << pattern << "\"";

  for (size_t i = 0; i < arrayRes.size(); ++i) {
    EXPECT_EQ(static_cast<int64_t>(arrayRes[i].offset), treeRes[i].offset)
        << "offset " << i << " differs for pattern \"" << pattern << "\"";
    EXPECT_EQ(static_cast<int64_t>(arrayRes[i].length), treeRes[i].length);
  }
}

TEST(Equivalence, RandomGenomePresentPatterns) {
  const std::string genome = randomGenome(20000, 1234);

  SuffixArray SA(genome);
  SuffixTree ST(genome);

  // Patterns drawn out of the text itself, so every one has at least one hit.
  // Short patterns exercise result collection over large subtrees; long ones
  // exercise deep edge traversal.
  std::mt19937 rng(5678);
  for (size_t patternLength : {1, 2, 3, 5, 8, 13, 21, 34}) {
    std::uniform_int_distribution<size_t> pick(0, genome.size() - patternLength);
    for (int trial = 0; trial < 25; ++trial) {
      const std::string pattern = genome.substr(pick(rng), patternLength);
      expectSameMatches(SA, ST, pattern);
    }
  }
}

TEST(Equivalence, RandomGenomeAbsentPatterns) {
  // 'Q' is outside the nucleotide alphabet, so every pattern containing it
  // must miss in both structures.
  const std::string genome = randomGenome(5000, 99);

  SuffixArray SA(genome);
  SuffixTree ST(genome);

  for (const std::string &pattern :
       {std::string("Q"), std::string("ACGQ"), std::string("QQQQQQ"),
        std::string("ACGTACGTQACGT")}) {
    expectSameMatches(SA, ST, pattern);
    EXPECT_EQ(SA.search(pattern).size(), 0u);
  }
}

TEST(Equivalence, HighlyRepetitiveText) {
  // Repetitive input is the worst case for the suffix tree's edge splitting
  // and produces the deepest subtrees for leaf collection.
  std::string genome;
  for (size_t i = 0; i < 500; ++i)
    genome += "ACGT";

  SuffixArray SA(genome);
  SuffixTree ST(genome);

  for (const std::string &pattern :
       {std::string("A"), std::string("ACGT"), std::string("ACGTACGT"),
        std::string("TACG"), std::string("GTA")}) {
    expectSameMatches(SA, ST, pattern);
  }
}

// Memory footprint
//
// Guards the property the benchmark reports on: the suffix array is a flat
// 8-bytes-per-character array, while the suffix tree costs orders of magnitude
// more. A regression that silently inflated the array would invalidate every
// tradeoff conclusion drawn from the benchmark.

TEST(MemoryFootprint, SuffixArrayIsOnePointerPerCharacter) {
  const std::string genome = randomGenome(4096, 7);
  SuffixArray SA(genome);

  EXPECT_EQ(SA.memoryBytes(), genome.size() * sizeof(size_t));
}

TEST(MemoryFootprint, SuffixTreeCostsMoreThanSuffixArray) {
  const std::string genome = randomGenome(4096, 7);

  SuffixArray SA(genome);
  SuffixTree ST(genome);

  EXPECT_GT(ST.memoryBytes(), SA.memoryBytes());

  // Ukkonen allocates at most 2n+1 nodes; anything beyond that means the
  // construction is creating nodes it should be reusing.
  EXPECT_LE(ST.nodeCount(), 2 * genome.size() + 1);
}
