#ifndef LM_BLANK__
#define LM_BLANK__

#include <limits>

#include <inttypes.h>
#include <math.h>

namespace lm {
namespace ngram {

/* Suppose "foo bar" appears with zero backoff but there is no trigram
 * beginning with these words.  Then, when scoring "foo bar", the model could
 * return out_state containing "bar" or even null context if "bar" also has no
 * backoff and is never followed by another word.  Then the backoff is set to
 * kNoExtensionBackoff.  If the n-gram might be extended, then out_state must
 * contain the full n-gram, in which case kExtensionBackoff is set.  In any
 * case, if an n-gram has non-zero backoff, the full state is returned so
 * backoff can be properly charged.  
 * These differ only in sign bit because the backoff is in fact zero in either
 * case.   
 */
const float kNoExtensionBackoff = -0.0;
const float kExtensionBackoff = 0.0;

inline void SetExtension(float &backoff) {
  if (backoff == kNoExtensionBackoff) backoff = kExtensionBackoff;
}

// This compiles down nicely.  
inline bool HasExtension(const float &backoff) {
  typedef union { float f; uint32_t i; } UnionValue;
  UnionValue compare, interpret;
  compare.f = kNoExtensionBackoff;
  interpret.f = backoff;
  return compare.i != interpret.i;
}

/* Suppose "foo bar baz quux" appears in the ARPA but not "bar baz quux" or
 * "baz quux" (because they were pruned).  1.2% of n-grams generated by SRI
 * with default settings on the benchmark data set are like this.  Since search
 * proceeds by finding "quux", "baz quux", "bar baz quux", and finally 
 * "foo bar baz quux" and the trie needs pointer nodes anyway, blanks are
 * inserted.  The blanks have probability kBlankProb and backoff kBlankBackoff.
 * A blank is recognized by kBlankProb in the probability field; kBlankBackoff
 * must be 0 so that inference asseses zero backoff from these blanks.  
 */
const float kBlankProb = -std::numeric_limits<float>::infinity();
const float kBlankBackoff = kNoExtensionBackoff;

} // namespace ngram
} // namespace lm
#endif // LM_BLANK__
