#include "SwampSeqLib/genome_mapper.h"
#include "SwampSeqLib/suffix_array.h"
#include "SwampSeqLib/suffix_tree.h"
#include <algorithm>
#include <chrono>
#include <ftxui/component/component.hpp>
#include <ftxui/component/screen_interactive.hpp>
#include <ftxui/dom/elements.hpp>
#include <iomanip>
#include <locale>
#include <sstream>
#include <string>
#include <thread>
#include <vector>

using namespace ftxui;

int main() {
  std::string searchPattern = "ATGC";
  std::string filePath = "bin/ecoli.fna";
  int algoChoice = 0;
  std::vector<std::string> algoEntries = {"Suffix Array", "Suffix Tree",
                                          "Compare Both"};

  std::string statusText = "Idle. Ready to benchmark.";
  std::string resultText = "";
  std::vector<std::string> previewLines;
  int selectedResult = 0;
  bool isRunning = false;
  auto screen = ScreenInteractive::TerminalOutput();

  Component inputFilepath = Input(&filePath, "e.g., ./data/genome.fna");
  Component inputPattern = Input(&searchPattern, "e.g., ATGC");
  Component algoRadiobox = Radiobox(&algoEntries, &algoChoice);
  Component resultList = Menu(&previewLines, &selectedResult);

  Component runButton = Button("Run Benchmark", [&] {
    if (isRunning)
      return;

    isRunning = true;
    statusText = "Mapping " + filePath + "...";
    resultText = "";
    previewLines.clear();

    std::thread([&, filePath, searchPattern, algoChoice] {
      GenomeMapper map;
      try {
        map.fromFile(filePath);
        if (!map.isValid()) {
          screen.Post([&] {
            statusText = "Could not open file.";
            isRunning = false;
          });
          screen.Post(Event::Custom);
          return;
        }
      } catch (const std::exception &e) {
        std::string errMsg = e.what();
        screen.Post([&, errMsg] {
          statusText = "Error: " + errMsg;
          isRunning = false;
        });
        screen.Post(Event::Custom);
        return;
      }

      std::vector<SearchResult> saRes;
      double saTimeMs = 0.0;
      int saMatchCount = 0;

      std::vector<STSearchResult> stRes;
      double stTimeMs = 0.0;
      int stMatchCount = 0;

      if (algoChoice == 0 || algoChoice == 2) {
        statusText = "Building Suffix Array...";
        screen.Post(Event::Custom);
        SuffixArray sa(map);
        statusText = "Searching Suffix Array...";
        screen.Post(Event::Custom);
        auto saStart = std::chrono::steady_clock::now();
        saRes = sa.search(searchPattern);
        auto saEnd = std::chrono::steady_clock::now();
        saMatchCount = static_cast<int>(saRes.size());
        saTimeMs =
            std::chrono::duration<double>(saEnd - saStart).count() * 1000.0;
      }

      if (algoChoice == 1 || algoChoice == 2) {
        statusText = "Building Suffix Tree...";
        screen.Post(Event::Custom);
        SuffixTree st(map);
        statusText = "Searching Suffix Tree...";
        screen.Post(Event::Custom);
        auto stStart = std::chrono::steady_clock::now();
        stRes = st.search(searchPattern);
        auto stEnd = std::chrono::steady_clock::now();
        stMatchCount = static_cast<int>(stRes.size());
        stTimeMs =
            std::chrono::duration<double>(stEnd - stStart).count() * 1000.0;
      }

      std::vector<std::string> localStrings;
      const char *data = map.data();

      int totalToDisplay = (algoChoice == 1) ? stMatchCount : saMatchCount;
      int previewLimit = std::min<int>(totalToDisplay, 2000);
      int contextLen = 25;

      for (int i = 0; i < previewLimit; i++) {
        size_t idx = (algoChoice == 1) ? stRes[i].offset : saRes[i].offset;

        size_t startIdx = (idx > (size_t)contextLen) ? (idx - contextLen) : 0;
        size_t endIdx = std::min<size_t>(
            map.size(), idx + searchPattern.size() + contextLen);

        std::string snippet(data + startIdx, endIdx - startIdx);
        localStrings.push_back("Offset(" + std::to_string(idx) + "): ..." +
                               snippet + "...");
      }

      std::stringstream ss;
      ss.imbue(std::locale::classic());
      ss << std::fixed << std::setprecision(1) << saTimeMs;
      std::string saFormattedTime = ss.str();

      ss.str("");
      ss.clear();
      ss << std::fixed << std::setprecision(1) << stTimeMs;
      std::string stFormattedTime = ss.str();

      screen.Post([&, saMatchCount, saFormattedTime, stMatchCount,
                   stFormattedTime, localStrings, algoChoice] {
        previewLines = std::move(localStrings);
        statusText = "Benchmark Complete!";

        std::string out;
        if (algoChoice == 0 || algoChoice == 2) {
          out += "Array: " + saFormattedTime + "ms (" +
                 std::to_string(saMatchCount) + ")";
        }
        if (algoChoice == 2)
          out += " | ";
        if (algoChoice == 1 || algoChoice == 2) {
          out += "Tree: " + stFormattedTime + "ms (" +
                 std::to_string(stMatchCount) + ")";
        }
        resultText = out;
        isRunning = false;
      });

      screen.Post(Event::Custom);
    }).detach();
  });

  auto layout = Container::Vertical(
      {inputFilepath, inputPattern, algoRadiobox, runButton, resultList});

  auto renderer = Renderer(layout, [&] {
    return vbox({text("SWAMP SEQUENCER") | bold | color(Color::Green) | center,
                 separator(),
                 hbox({text("Genome File Path : "), inputFilepath->Render()}),
                 hbox({text("Search Pattern   : "), inputPattern->Render()}),
                 separator(), text("Select Data Structure:"),
                 algoRadiobox->Render(), separator(),
                 runButton->Render() | center, separator(),
                 text(statusText) | color(Color::Yellow),
                 text(resultText) | color(Color::Cyan) | bold,
                 vbox({text("Results (Arrows to scroll)") | dim, separator(),
                       resultList->Render() | vscroll_indicator | frame |
                           size(HEIGHT, LESS_THAN, 15)}) |
                     border | flex}) |
           borderHeavy;
  });

  screen.Loop(renderer);
  return 0;
}
