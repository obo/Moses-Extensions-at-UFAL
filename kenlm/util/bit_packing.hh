#ifndef UTIL_BIT_PACKING__
#define UTIL_BIT_PACKING__

/* Bit-level packing routines */

#include <assert.h>
#ifdef __APPLE__
#include <architecture/byte_order.h>
#elif __linux__
#include <endian.h>
#else
#include <arpa/nameser_compat.h>
#endif 

#include <inttypes.h>

namespace util {

/* WARNING WARNING WARNING:
 * The write functions assume that memory is zero initially.  This makes them
 * faster and is the appropriate case for mmapped language model construction.
 * These routines assume that unaligned access to uint64_t is fast and that
 * storage is little endian.  This is the case on x86_64.  I'm not sure how 
 * fast unaligned 64-bit access is on x86 but my target audience is large
 * language models for which 64-bit is necessary.  
 *
 * Call the BitPackingSanity function to sanity check.  Calling once suffices,
 * but it may be called multiple times when that's inconvenient.  
 */

inline uint8_t BitPackShift(uint8_t bit, uint8_t length) {
// Fun fact: __BYTE_ORDER is wrong on Solaris Sparc, but the version without __ is correct.  
#if BYTE_ORDER == LITTLE_ENDIAN
  return bit;
#elif BYTE_ORDER == BIG_ENDIAN
  return 64 - length - bit;
#else
#error "Bit packing code isn't written for your byte order."
#endif
}

/* Pack integers up to 57 bits using their least significant digits. 
 * The length is specified using mask:
 * Assumes mask == (1 << length) - 1 where length <= 57.   
 */
inline uint64_t ReadInt57(const void *base, uint8_t bit, uint8_t length, uint64_t mask) {
  return (*reinterpret_cast<const uint64_t*>(base) >> BitPackShift(bit, length)) & mask;
}
/* Assumes value < (1 << length) and length <= 57.
 * Assumes the memory is zero initially. 
 */
inline void WriteInt57(void *base, uint8_t bit, uint8_t length, uint64_t value) {
  *reinterpret_cast<uint64_t*>(base) |= (value << BitPackShift(bit, length));
}

typedef union { float f; uint32_t i; } FloatEnc;

inline float ReadFloat32(const void *base, uint8_t bit) {
  FloatEnc encoded;
  encoded.i = *reinterpret_cast<const uint64_t*>(base) >> BitPackShift(bit, 32);
  return encoded.f;
}
inline void WriteFloat32(void *base, uint8_t bit, float value) {
  FloatEnc encoded;
  encoded.f = value;
  WriteInt57(base, bit, 32, encoded.i);
}

const uint32_t kSignBit = 0x80000000;

inline float ReadNonPositiveFloat31(const void *base, uint8_t bit) {
  FloatEnc encoded;
  encoded.i = *reinterpret_cast<const uint64_t*>(base) >> BitPackShift(bit, 31);
  // Sign bit set means negative.  
  encoded.i |= kSignBit;
  return encoded.f;
}
inline void WriteNonPositiveFloat31(void *base, uint8_t bit, float value) {
  FloatEnc encoded;
  encoded.f = value;
  encoded.i &= ~kSignBit;
  WriteInt57(base, bit, 31, encoded.i);
}

void BitPackingSanity();

// Return bits required to store integers upto max_value.  Not the most
// efficient implementation, but this is only called a few times to size tries. 
uint8_t RequiredBits(uint64_t max_value);

struct BitsMask {
  void FromMax(uint64_t max_value) {
    bits = RequiredBits(max_value);
    mask = (1 << bits) - 1;
  }
  uint8_t bits;
  uint64_t mask;
};

} // namespace util

#endif // UTIL_BIT_PACKING__
