#pragma once
#include <string>

// Define memory mapper based on OS
#if defined(_WIN32) || defined(_WIN64)
#define OS_WINDOWS
#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <windows.h>
#else
#define OS_POSIX
#include <fcntl.h>    // open()
#include <sys/mman.h> // mmap(), munmap()
#include <sys/stat.h> // fstat()
#include <unistd.h>
#endif

class GenomeMapper {
#ifdef OS_WINDOWS
  HANDLE _fileHandle = INVALID_HANDLE_VALUE;
#elif defined(OS_POSIX)
  int _fileDescriptor = -1;
#endif
  size_t _mappingSize = 0;
  size_t _fileSize = 0;
  void *_mappingBase = nullptr;
  void *_mappedData = nullptr;
  bool _isValid = false;

public:
  explicit GenomeMapper(const std::string &filePath);
  GenomeMapper();
  ~GenomeMapper();

  GenomeMapper(const GenomeMapper &) = delete;
  GenomeMapper &operator=(const GenomeMapper &) = delete;

  const char *data() const;
  char *data();
  size_t size() const;
  bool isValid() const;
  void fromFile(const std::string &filePath);
};
