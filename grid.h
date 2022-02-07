#include "common.h"
#include <vector>

#ifndef numBins
#define numBins 16
#endif

typedef std::vector<particle_t*> bin_t;

typedef struct
{
    
    bin_t collection[numBins][numBins];

    // Init() initializes this grid_t object by instantiating
    // a bin_t objects at each position (i, j).
    void init();

} grid_t;