#include "bron_kerbosch.h"
#include "BKGraph.hpp"
#include "clique.h"
//#include "cut.h"

struct _BronKerbosch
{
    BKGraph *bkg;
    CliqueSet *clqSet;
//    CutPool *cutPool;
//    char usePool;
};

//BronKerbosch *bk_create(const CGraph *cgraph, char usePool)
BronKerbosch *bk_create(const CGraph *cgraph)
{
    BronKerbosch *result = (BronKerbosch *) xmalloc( sizeof(BronKerbosch) );
    result->bkg = new BKGraph(cgraph);
    result->clqSet = NULL;

//    result->usePool = usePool;
//    if(result->usePool)
//        result->cutPool = cut_pool_create()

    return result;
}

int bk_run(BronKerbosch *bk, const int minViol, const double timeLimit)
{
    if (bk->clqSet)
    {
        clq_set_free( &(bk->clqSet) );
        bk->clqSet = NULL;
    }

    int status = bk->bkg->execute( minViol, timeLimit );
    bk->clqSet =bk->bkg->getCliqueSet();

    return status;
}

const CliqueSet *bk_get_clq_set( BronKerbosch *bk )
{
    return bk->clqSet;
}

int bk_get_max_weight( BronKerbosch *bk )
{
    return bk->bkg->getMaxWeight();
}

void bk_free( BronKerbosch *bk )
{
    delete bk->bkg;
	bk->bkg = NULL;
    free( bk );
}
