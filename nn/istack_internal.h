#if !defined(_ISTACK_INTERNAL_H)
#define _ISTACK_INTERNAL_H

#include "istack.h"

struct istack {
    int n;
    int nallocated;
    int* v;
};

#endif                          /* _ISTACK_INTERNAL_H */
