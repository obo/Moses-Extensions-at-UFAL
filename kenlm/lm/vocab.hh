#ifndef LM_VOCAB__
#define LM_VOCAB__

#include "lm/enumerate_vocab.hh"
#include "lm/virtual_interface.hh"
#include "util/key_value_packing.hh"
#include "util/probing_hash_table.hh"
#include "util/sorted_uniform.hh"
#include "util/string_piece.hh"

#include <string>
#include <vector>

namespace lm {
class ProbBackoff;

namespace ngram {
class Config;
class EnumerateVocab;

namespace detail {
uint64_t HashForVocab(const char *str, std::size_t len);
inline uint64_t HashForVocab(const StringPiece &str) {
  return HashForVocab(str.data(), str.length());
}
} // namespace detail

class WriteWordsWrapper : public EnumerateVocab {
  public:
    WriteWordsWrapper(EnumerateVocab *inner);

    ~WriteWordsWrapper();
    
    void Add(WordIndex index, const StringPiece &str);

    void Write(int fd);

  private:
    EnumerateVocab *inner_;

    std::string buffer_;
};

// Vocabulary based on sorted uniform find storing only uint64_t values and using their offsets as indices.  
class SortedVocabulary : public base::Vocabulary {
  private:
    // Sorted uniform requires a GetKey function.  
    struct Entry {
      uint64_t GetKey() const { return key; }
      uint64_t key;
      bool operator<(const Entry &other) const {
        return key < other.key;
      }
    };

  public:
    SortedVocabulary();

    WordIndex Index(const StringPiece &str) const {
      const Entry *found;
      if (util::SortedUniformFind<const Entry *, uint64_t>(begin_, end_, detail::HashForVocab(str), found)) {
        return found - begin_ + 1; // +1 because <unk> is 0 and does not appear in the lookup table.
      } else {
        return 0;
      }
    }

    static size_t Size(std::size_t entries, const Config &config);

    // Everything else is for populating.  I'm too lazy to hide and friend these, but you'll only get a const reference anyway.
    void SetupMemory(void *start, std::size_t allocated, std::size_t entries, const Config &config);

    void ConfigureEnumerate(EnumerateVocab *to, std::size_t max_entries);

    WordIndex Insert(const StringPiece &str);

    // Reorders reorder_vocab so that the IDs are sorted.  
    void FinishedLoading(ProbBackoff *reorder_vocab);

    bool SawUnk() const { return saw_unk_; }

    void LoadedBinary(int fd, EnumerateVocab *to);

  private:
    Entry *begin_, *end_;

    bool saw_unk_;

    EnumerateVocab *enumerate_;

    // Actual strings.  Used only when loading from ARPA and enumerate_ != NULL 
    std::vector<std::string> strings_to_enumerate_;
};

// Vocabulary storing a map from uint64_t to WordIndex. 
class ProbingVocabulary : public base::Vocabulary {
  public:
    ProbingVocabulary();

    WordIndex Index(const StringPiece &str) const {
      Lookup::ConstIterator i;
      return lookup_.Find(detail::HashForVocab(str), i) ? i->GetValue() : 0;
    }

    static size_t Size(std::size_t entries, const Config &config);

    // Everything else is for populating.  I'm too lazy to hide and friend these, but you'll only get a const reference anyway.
    void SetupMemory(void *start, std::size_t allocated, std::size_t entries, const Config &config);

    void ConfigureEnumerate(EnumerateVocab *to, std::size_t max_entries);

    WordIndex Insert(const StringPiece &str);

    void FinishedLoading(ProbBackoff *reorder_vocab);

    bool SawUnk() const { return saw_unk_; }

    void LoadedBinary(int fd, EnumerateVocab *to);

  private:
    // std::identity is an SGI extension :-(
    struct IdentityHash : public std::unary_function<uint64_t, std::size_t> {
      std::size_t operator()(uint64_t arg) const { return static_cast<std::size_t>(arg); }
    };

    typedef util::ProbingHashTable<util::ByteAlignedPacking<uint64_t, WordIndex>, IdentityHash> Lookup;

    Lookup lookup_;

    WordIndex available_;

    bool saw_unk_;

    EnumerateVocab *enumerate_;
};

} // namespace ngram
} // namespace lm

#endif // LM_VOCAB__
