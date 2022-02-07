#ifndef GRID_H
#define GRID_H
#include "grid.h"

void grid_t::init()
{

    for(int i = 0; i < numBins; i++)
    {

        for(int j = 0; j < numBins; j++)
        {

            collection[i][j] = bin_t();

        }
    }
}

#endif