#include "SwampSeqLib/genome_mapper.h"

#include <cctype>
#include <cstddef>
#include <cstring>
#include <iostream>
#include <stdexcept>
#ifdef OS_WINDOWS
#include <errhandlingapi.h>
#include <fileapi.h>
#include <handleapi.h>
#include <memoryapi.h>
#include <winnt.h>
#endif
#include <stdexcept>
#include <string>

void GenomeMapper::fromFile(const std::string& filePath) {
#ifdef OS_WINDOWS
  // Open the source file, then copy it into an anonymous writable mapping.
  _fileHandle = CreateFileA(filePath.c_str(), GENERIC_READ, FILE_SHARE_READ,
                            NULL, OPEN_EXISTING, FILE_ATTRIBUTE_NORMAL, NULL);
  if (_fileHandle == INVALID_HANDLE_VALUE) return;

  LARGE_INTEGER size;
  if (!GetFileSizeEx(_fileHandle, &size)) return;
  _mappingSize = static_cast<size_t>(size.QuadPart) + 1;

  HANDLE memoryHandle =
      CreateFileMappingA(INVALID_HANDLE_VALUE, NULL, PAGE_READWRITE, 0,
                         static_cast<DWORD>(_mappingSize), NULL);
  if (memoryHandle == NULL) return;

  _mappingBase =
      MapViewOfFile(memoryHandle, FILE_MAP_WRITE, 0, 0, _mappingSize);
  CloseHandle(memoryHandle);
  if (_mappingBase == nullptr) return;

  auto*  buffer    = static_cast<char*>(_mappingBase);
  size_t totalRead = 0;
  while (totalRead < _mappingSize - 1) {
    unsigned long long kMaxCopyChunk = 1 << 24;
    const size_t       chunkSize =
        std::min(_mappingSize - 1 - totalRead, kMaxCopyChunk);
    DWORD bytesRead = 0;
    if (!ReadFile(_fileHandle, buffer + totalRead,
                  static_cast<DWORD>(chunkSize), &bytesRead, NULL)) {
      UnmapViewOfFile(_mappingBase);
      _mappingBase = nullptr;
      throw std::runtime_error(std::to_string(GetLastError()));
    }
    if (bytesRead == 0) break;
    totalRead += static_cast<size_t>(bytesRead);
  }

  buffer[totalRead] = '\0';
  _mappedData       = _mappingBase;
  _fileSize         = totalRead;
  _isValid          = true;

#elif defined(OS_POSIX)
  // Create a writable copy-on-write view of the source file.
  _fileDescriptor = open(filePath.c_str(), O_RDONLY);
  if (_fileDescriptor == -1) return;

  struct stat sb;
  if (fstat(_fileDescriptor, &sb) == -1) return;
  if (sb.st_size <= 0) return;

  _mappingSize = static_cast<size_t>(sb.st_size);
  _mappingBase = mmap(nullptr, _mappingSize, PROT_READ | PROT_WRITE,
                      MAP_PRIVATE, _fileDescriptor, 0);
  if (_mappingBase == MAP_FAILED) {
    _mappingBase = nullptr;
    return;
  }

  _mappedData = _mappingBase;
  _fileSize   = _mappingSize;
  _isValid    = true;
#endif
  if (!_isValid) return;

  // If we have a header, skip it before compacting sequence bytes.
  char* data    = this->data();
  char* dataEnd = static_cast<char*>(_mappingBase) + _fileSize;
  if (data[0] == '>') {
    char* firstNewLine = static_cast<char*>(
        std::memchr(data, '\n', static_cast<size_t>(dataEnd - data)));
    if (firstNewLine) {
      data = firstNewLine + 1;
    }
  }

  // Remove headers, newlines, and any non-sequence characters in place.
  char* src = data;
  char* dst = data;

  while (src < dataEnd) {
    if (*src == '>') {
      char* nextNewLine = static_cast<char*>(
          std::memchr(src, '\n', static_cast<size_t>(dataEnd - src)));
      if (nextNewLine == nullptr) break;
      src = nextNewLine + 1;
      continue;
    }

    char upper = std::toupper(static_cast<unsigned char>(*src));
    if (upper == 'A' || upper == 'C' || upper == 'T' || upper == 'G' ||
        upper == 'N') {
      *dst++ = upper;
    }
    ++src;
  }
  _mappedData = data;
  _fileSize   = static_cast<size_t>(dst - data);
}

GenomeMapper::GenomeMapper(const std::string& filePath) { fromFile(filePath); }
GenomeMapper::GenomeMapper() {};

GenomeMapper::~GenomeMapper() {
#ifdef OS_WINDOWS
  if (_mappingBase) UnmapViewOfFile(_mappingBase);
  if (_fileHandle != INVALID_HANDLE_VALUE) CloseHandle(_fileHandle);
#elif defined(OS_POSIX)
  if (_mappingBase) munmap(_mappingBase, _mappingSize);
  if (_fileDescriptor != -1) close(_fileDescriptor);
#endif
}

const char* GenomeMapper::data() const {
  return static_cast<const char*>(_mappedData);
}
char*  GenomeMapper::data() { return static_cast<char*>(_mappedData); }
size_t GenomeMapper::size() const { return _fileSize; }
bool   GenomeMapper::isValid() const { return _isValid; }
