#ifndef __MURMURHASH3_H__
#define __MURMURHASH3_H__

#include <stdlib.h>
#include <stdint.h>

uint32_t murmur3_hash(const char* data, size_t len);

#endif //__MURMURHASH3_H__
