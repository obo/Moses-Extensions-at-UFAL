#ifndef UTIL_MMAP__
#define UTIL_MMAP__
// Utilities for mmaped files.  

#include "util/scoped.hh"

#include <cstddef>

#include <sys/types.h>

namespace util {

// (void*)-1 is MAP_FAILED; this is done to avoid including the mmap header here.  
class scoped_mmap {
  public:
    scoped_mmap() : data_((void*)-1), size_(0) {}
    scoped_mmap(void *data, std::size_t size) : data_(data), size_(size) {}
    ~scoped_mmap();

    void *get() const { return data_; }

    const char *begin() const { return reinterpret_cast<char*>(data_); }
    const char *end() const { return reinterpret_cast<char*>(data_) + size_; }
    std::size_t size() const { return size_; }

    void reset(void *data, std::size_t size) {
      scoped_mmap other(data_, size_);
      data_ = data;
      size_ = size;
    }

    void reset() {
      reset((void*)-1, 0);
    }

  private:
    void *data_;
    std::size_t size_;

    scoped_mmap(const scoped_mmap &);
    scoped_mmap &operator=(const scoped_mmap &);
};

/* For when the memory might come from mmap, new char[], or malloc.  Uses NULL
 * and 0 for blanks even though mmap signals errors with (void*)-1).  The reset
 * function checks that blank for mmap.  
 */
class scoped_memory {
  public:
    typedef enum {MMAP_ALLOCATED, ARRAY_ALLOCATED, MALLOC_ALLOCATED, NONE_ALLOCATED} Alloc;

    scoped_memory() : data_(NULL), size_(0), source_(NONE_ALLOCATED) {}

    ~scoped_memory() { reset(); }

    void *get() const { return data_; }
    const char *begin() const { return reinterpret_cast<char*>(data_); }
    const char *end() const { return reinterpret_cast<char*>(data_) + size_; }
    std::size_t size() const { return size_; }

    Alloc source() const { return source_; }

    void reset() { reset(NULL, 0, NONE_ALLOCATED); }

    void reset(void *data, std::size_t size, Alloc from);

    // realloc allows the current data to escape hence the need for this call
    // If realloc fails, destroys the original too and get() returns NULL.
    void call_realloc(std::size_t to);

  private:

    void *data_;
    std::size_t size_;

    Alloc source_;

    scoped_memory(const scoped_memory &);
    scoped_memory &operator=(const scoped_memory &);
};

struct scoped_mapped_file {
  scoped_fd fd;
  scoped_mmap mem;
};

// Wrapper around mmap to check it worked and hide some platform macros.  
void *MapOrThrow(std::size_t size, bool for_write, int flags, bool prefault, int fd, off_t offset = 0);
void *MapForRead(std::size_t size, bool prefault, int fd, off_t offset = 0);

void *MapAnonymous(std::size_t size);

// Open file name with mmap of size bytes, all of which are initially zero.  
void MapZeroedWrite(const char *name, std::size_t size, scoped_fd &file, scoped_mmap &mem);
inline void MapZeroedWrite(const char *name, std::size_t size, scoped_mapped_file &out) {
  MapZeroedWrite(name, size, out.fd, out.mem);
}
 
} // namespace util

#endif // UTIL_SCOPED__
