#ifndef UTIL_FILE_PIECE__
#define UTIL_FILE_PIECE__

#include "util/scoped.hh"
#include "util/string_piece.hh"

#include <exception>
#include <string>

#include <cstddef>

namespace util {

class EndOfFileException : public std::exception {
  public:
    EndOfFileException() throw() {}

    ~EndOfFileException() throw() {}

    const char *what() const throw() { return "End of file."; }
};

class ParseNumberException : public std::exception {
  public:
    explicit ParseNumberException(StringPiece value) throw();

    ~ParseNumberException() throw() {}

    const char *what() const throw() { return what_.c_str(); }

  private:
    std::string what_;
};

class FilePiece {
  public:
    // 32 MB default.
    explicit FilePiece(const char *file, off_t min_buffer = 33554432);
     
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
    
  private:
    StringPiece Consume(const char *to) {
      StringPiece ret(position_, to - position_);
      position_ = to;
      return ret;
    }

    const char *FindDelimiterOrEOF() throw(EndOfFileException);

    void Shift() throw (EndOfFileException);

    const char *position_, *last_space_, *position_end_;

    scoped_fd file_;
    const off_t total_size_;
    const off_t page_;

    off_t default_map_size_;
    off_t mapped_offset_;

    // Order matters: file_ should always be destroyed after this.
    scoped_mmap data_;

    bool at_end_;
};

} // namespace util

#endif // UTIL_FILE_PIECE__
