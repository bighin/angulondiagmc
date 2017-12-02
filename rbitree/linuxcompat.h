#ifndef __LINUXCOMPAT_H__
#define __LINUXCOMPAT_H__

#include <stdio.h>
#include <stdbool.h>
#include <stddef.h>

#define WRITE_ONCE(x, val) 	x=(val)

#define container_of(ptr, type, member) ({                      \
        const typeof( ((type *)0)->member ) *__mptr = (ptr);    \
        (type *)( (char *)__mptr - offsetof(type,member) );})

#endif //__LINUXCOMPAT_H__
