# SwampSequencer
[![CI](https://github.com/Jaefae/Swamp-Sequencers/actions/workflows/ci.yml/badge.svg)](https://github.com/Jaefae/Swamp-Sequencers/actions/workflows/ci.yml)
[![sanitizers](https://github.com/Jaefae/Swamp-Sequencers/actions/workflows/sanitizers.yml/badge.svg)](https://github.com/Jaefae/Swamp-Sequencers/actions/workflows/sanitizers.yml)

## Overview
Given a massive DNA genome database, finding where a specific gene sequence pattern appears is the core operation behind disease detection, ancestry matching, and drug research.

This project benchmarks the performance, memory efficiency, and build constraints of two O(n) string-matching data structures: the Suffix Array (using the SA-IS algorithm) and the Suffix Tree (using Ukkonen's algorithm). We test these structures against real genome data, specifically NCBI GenBank's E. Coli strain Sakai sample chromosome found [here](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000008865.2/).
## Getting Started
### 1. Clone the repository
```Bash
git clone https://github.com/Jaefae/Swamp-Sequencers
cd Swamp-Sequencers
```

### 2. Build the project
The project uses an out-of-source build. The first time the project is built it will also download the E. Coli reference genome and place it next to the binary.
```Bash
# Swamp-Sequencers/
mkdir build
cd build
# Swamp-Sequencers/build
cmake ..
make
```
### Usage
Run the TUI interface to query the dataset.
```Bash
# Swamp-Sequencers/build
bin/SwampSequencer 
```
Or run the automated Google Test suite
```Bash
# Swamp-Sequencers/build
bin/SwampTests
```

## Results
See [`bench/README.md`](bench/README.md) for a full writeup.
Reference run in [`results/ecoli-50.csv`](results/ecoli-50.csv): the full E. coli
reference genome (5,594,605 bases after filtering to ACGT), 50-step sweep,
1,000 queries per point, clang 22.1.4 `-O3`, i7-9700, 32GB RAM.

| | Suffix array | Suffix tree | |
|---|---|---|---|
| Memory | 44.8 MB | 1,931 MB | **43.1× larger** |
| Build time | 848 ms | 8,373 ms | **9.9× slower** |
| Query (m = 64) | 1.18 µs | 3.03 µs | **2.6× slower** |

### Memory

<picture>
  <source media="(prefers-color-scheme: dark)" srcset="bench/results/memory-dark.svg">
  <img alt="Index memory against text length: the suffix array grows to 45 MB while the suffix tree reaches 1,931 MB, 43.1 times larger." src="bench/results/memory-light.svg" width="100%">
</picture>

Both are linear in n, as O(n) construction requires — the constant factor is the
whole story. The array is exactly 8 bytes per base; the tree needs ~345.

### Build time

<picture>
  <source media="(prefers-color-scheme: dark)" srcset="bench/results/build-dark.svg">
  <img alt="Build time against text length: the suffix array reaches 848 ms while the suffix tree reaches 8,373 ms, 9.9 times slower." src="bench/results/build-light.svg" width="100%">
</picture>

### Query latency

<picture>
  <source media="(prefers-color-scheme: dark)" srcset="bench/results/query-dark.svg">
  <img alt="Query latency against pattern length: the suffix array stays near 1.2 microseconds while the suffix tree stays near 3 microseconds." src="bench/results/query-light.svg" width="100%">
</picture>

## Build options

| Option | Default | Effect |
|---|---|---|
| `SWAMP_BUILD_TUI` | `ON` | Build the FTXUI front-end and download the reference genome |
| `SWAMP_BUILD_TESTS` | `ON` | Build the GoogleTest suite |
| `SWAMP_BUILD_BENCH` | `OFF` | Build the `SwampBench` benchmark |
| `SWAMP_ENABLE_ASAN` | `OFF` | Build with AddressSanitizer + UndefinedBehaviorSanitizer |

Turning the TUI off skips both FTXUI and the ~5 MB NCBI genome download, which
makes a test-only build much faster:

```Bash
cmake -B build -DSWAMP_BUILD_TUI=OFF
cmake --build build --parallel
ctest --test-dir build --output-on-failure
```

## Continuous integration

Two workflows run on every push and pull request to `main`, each reporting its
own badge above.

[`.github/workflows/ci.yml`](.github/workflows/ci.yml):

- **test**: builds and runs the suite on Linux, macOS, and Windows in both
  Debug and Release.
- **benchmark**: runs a small sweep and publishes the CSV as a build artifact
  so the performance tradeoff is tracked across commits.

[`.github/workflows/sanitizers.yml`](.github/workflows/sanitizers.yml):

- **asan + ubsan**: reruns the suite with `SWAMP_ENABLE_ASAN=ON`. Both indexes
  do heavy raw-pointer and index arithmetic, which is exactly the class of bug a
  plain pass/fail run will not surface. It is a separate workflow because GitHub
  badges report per workflow, combined into `ci.yml` a sanitizer
  failure would be indistinguishable from a plain test failure.

Alongside the per-structure unit tests, the suite includes **differential
tests**: the suffix array and suffix tree are independent implementations of the
same contract, so any query where they disagree is a bug in one of them.

## Attribution
    - Michael Amiot
    - Jack Hendrix
    - Sebastian Mejia
