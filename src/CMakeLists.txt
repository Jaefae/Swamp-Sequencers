# NOTE: no project()/cmake_minimum_required here — this is a subdirectory of the
# top-level SwampSequencer project and inherits its C++20 settings. A nested
# project() call previously reset the standard back to C++17 for this target.

add_library(SwampSeqLib STATIC "")
target_sources(SwampSeqLib PRIVATE
  suffix_array.cpp
  suffix_tree.cpp
  genome_mapper.cpp
)

target_include_directories(SwampSeqLib PUBLIC ${CMAKE_SOURCE_DIR}/include)

# -Wall/-Wextra/-Wpedantic are GCC/Clang spellings; MSVC needs /W4.
if(MSVC)
  target_compile_options(SwampSeqLib PRIVATE /W4)
else()
  target_compile_options(SwampSeqLib PRIVATE -Wall -Wextra -Wpedantic)
endif()

# Skipped when SWAMP_BUILD_TUI=OFF so CI can build the library and tests
# without FTXUI or the genome download.
if(SWAMP_BUILD_TUI)
  add_executable(SwampSequencer "")
  target_sources(SwampSequencer PRIVATE main.cpp)

  FetchContent_Declare(
      ecoli_genome
      URL "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/GCF_000008865.2/download?include_annotation_type=GENOME_FASTA"
      DOWNLOAD_NO_EXTRACT FALSE  # NCBI sends a ZIP, set to TRUE if you want the raw zip
  )

  FetchContent_MakeAvailable(ecoli_genome)

  set(GENOME_SOURCE_FILE "${ecoli_genome_SOURCE_DIR}/ncbi_dataset/data/GCF_000008865.2/GCF_000008865.2_ASM886v2_genomic.fna")

  # This ensures that 'ecoli.fna' is always sitting next to SwampSequencer.exe
  add_custom_command(
      TARGET SwampSequencer POST_BUILD
      COMMAND ${CMAKE_COMMAND} -E copy
              ${GENOME_SOURCE_FILE}
              ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/ecoli.fna
      COMMENT "Copying E. coli genome to binary directory..."
  )

  target_link_libraries(SwampSequencer PRIVATE SwampSeqLib)
  target_link_libraries(SwampSequencer PRIVATE
    ftxui::screen
    ftxui::dom
    ftxui::component
  )
endif()
