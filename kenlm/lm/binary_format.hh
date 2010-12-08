#ifndef LM_BINARY_FORMAT__
#define LM_BINARY_FORMAT__

#include "lm/config.hh"
#include "lm/read_arpa.hh"

#include "util/file_piece.hh"
#include "util/mmap.hh"
#include "util/scoped.hh"

#include <cstddef>
#include <vector>

#include <inttypes.h>

namespace lm {
namespace ngram {

typedef enum {HASH_PROBING=0, HASH_SORTED=1, TRIE_SORTED=2} ModelType;

struct FixedWidthParameters {
  unsigned char order;
  float probing_multiplier;
  // What type of model is this?  
  ModelType model_type;
  // Does the end of the file have the actual strings in the vocabulary?   
  bool has_vocabulary;
};

struct Parameters {
  FixedWidthParameters fixed;
  std::vector<uint64_t> counts;
};

struct Backing {
  // File behind memory, if any.  
  util::scoped_fd file;
  // Raw block of memory backing the language model data structures
  util::scoped_memory memory;
};

namespace detail {

bool IsBinaryFormat(int fd);

void ReadHeader(int fd, Parameters &params);

void MatchCheck(ModelType model_type, const Parameters &params);

uint8_t *SetupBinary(const Config &config, const Parameters &params, std::size_t memory_size, Backing &backing);

uint8_t *SetupZeroed(const Config &config, ModelType model_type, const std::vector<uint64_t> &counts, std::size_t memory_size, Backing &backing);

void ComplainAboutARPA(const Config &config, ModelType model_type);

} // namespace detail

bool RecognizeBinary(const char *file, ModelType &recognized);

template <class To> void LoadLM(const char *file, const Config &config, To &to) {
  Backing &backing = to.MutableBacking();
  backing.file.reset(util::OpenReadOrThrow(file));

  Parameters params;

  try {
    if (detail::IsBinaryFormat(backing.file.get())) {
      detail::ReadHeader(backing.file.get(), params);
      detail::MatchCheck(To::kModelType, params);
      // Replace the probing_multiplier.  
      Config new_config(config);
      new_config.probing_multiplier = params.fixed.probing_multiplier;
      std::size_t memory_size = To::Size(params.counts, new_config);
      uint8_t *start = detail::SetupBinary(new_config, params, memory_size, backing);
      to.InitializeFromBinary(start, params, new_config, backing.file.get());
    } else {
      detail::ComplainAboutARPA(config, To::kModelType);
      util::FilePiece f(backing.file.release(), file, config.messages);
      ReadARPACounts(f, params.counts);
      std::size_t memory_size = To::Size(params.counts, config);
      uint8_t *start = detail::SetupZeroed(config, To::kModelType, params.counts, memory_size, backing);

      to.InitializeFromARPA(file, f, start, params, config);
    }
  } catch (util::Exception &e) {
    e << " in file " << file;
    throw;
  }

}

} // namespace ngram
} // namespace lm
#endif // LM_BINARY_FORMAT__
