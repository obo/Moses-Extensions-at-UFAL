#ifndef UTIL_FILE_PIECE__
#define UTIL_FILE_PIECE__

#include "util/ersatz_progress.hh"
#include "util/exception.hh"
#include "util/mmap.hh"
#include "util/scoped.hh"
#include "util/string_piece.hh"

#include <string>

#include <cstddef>

namespace util {

class EndOfFileException : public Exception {
  public:
    EndOfFileException() throw();
    ~EndOfFileException() throw();
};

class ParseNumberException : public Exception {
  public:
    explicit ParseNumberException(StringPiece value) throw();
    ~ParseNumberException() throw() {}
};

int OpenReadOrThrow(const char *name);

// Return value for SizeFile when it can't size properly.  
const off_t kBadSize = -1;
off_t SizeFile(int fd);

class FilePiece {
  public:
    // 32 MB default.
    explicit FilePiece(const char *file, std::ostream *show_progress = NULL, off_t min_buffer = 33554432);
    // Takes ownership of fd.  name is used for messages.  
    explicit FilePiece(const char *name, int fd, std::ostream *show_progress = NULL, off_t min_buffer = 33554432);
     
    char get() throw(EndOfFileException) { 
      if (position_ == position_end_) Shift();
      return *(position_++);
    }

    // Memory backing the returned StringPiece may vanish on the next call.  
    // Leaves the delimiter, if any, to be returned by get().
    StringPiece ReadDelimited() throw(EndOfFileException) {
      SkipSpaces();
      return Consume(FindDelimiterOrEOF());
    }
    // Unlike ReadDelimited, this includes leading spaces and consumes the delimiter.
    // It is similar to getline in that way.  
    StringPiece ReadLine(char delim = '\n') throw(EndOfFileException);

    float ReadFloat() throw(EndOfFileException, ParseNumberException);

    void SkipSpaces() throw (EndOfFileException);

    off_t Offset() const {
      return position_ - data_.begin() + mapped_offset_;
    }

    // Only for testing.  
    void ForceFallbackToRead() {
      fallback_to_read_ = true;
    }
    
  private:
    void Initialize(const char *name, std::ostream *show_progress, off_t min_buffer);

    StringPiece Consume(const char *to) {
      StringPiece ret(position_, to - position_);
      position_ = to;
      return ret;
    }

    const char *FindDelimiterOrEOF() throw(EndOfFileException);

    void Shift() throw (EndOfFileException);
    // Backends to Shift().
    void MMapShift(off_t desired_begin) throw ();
    void ReadShift(off_t desired_begin) throw ();

    const char *position_, *last_space_, *position_end_;

    scoped_fd file_;
    const off_t total_size_;
    const off_t page_;

    size_t default_map_size_;
    off_t mapped_offset_;

    // Order matters: file_ should always be destroyed after this.
    scoped_memory data_;

    bool at_end_;
    bool fallback_to_read_;

    ErsatzProgress progress_;
};

} // namespace util

#endif // UTIL_FILE_PIECE__
