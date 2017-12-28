#ifndef BUILDCGRAPH_H_INCLUDED
#define BUILDCGRAPH_H_INCLUDED

extern "C"
{
#include "cgraph.h"
}

CGraph *build_cgraph_osi( const void *_solver );

CGraph* build_cgraph_lp(const void *_lp);

#endif