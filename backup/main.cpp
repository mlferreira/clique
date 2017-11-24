#include <cstdio>
#include <iostream>
#include <cstdlib>
#include <algorithm>
#include <fstream>

using namespace std;
extern "C"
{
    #include "lp.h"
    #include "vint_set.h"
    #include "clique_extender.h"
    #include "build_cgraph.h"
}


int main( int argc, const char **argv )
{
    LinearProgram *lp = lp_create();
    lp_read(lp, argv[1]);

    int nCliques; int* cliques; int sizecliques;

    CGraph *cg = build_cgraph_lp(lp);
    const int nVertices = cgraph_size(cg);
    int costs[nVertices];

    fill(costs, costs + nVertices, 1); /* preenchendo o custo de cada variavel com valor 1 */

    
    ifstream file(argv[2]);
    ofstream saida("saida.txt");
    file >> nCliques;

    for(int k = 0; k < nCliques; k++) {


        file >> sizecliques;
        cliques = new int[sizecliques];
        for(int j=0; j<sizecliques; j++) {
            file >> cliques[j];
        }

        if(sizecliques < 100) {

            IntSet clique;
            vint_set_init(&clique);
            vint_set_add(&clique, cliques, sizecliques);

            CliqueExtender *clqe = clqe_create();

            /* Funcao retorna 0 se o clique nao foi estendido ou o numero de cliques gerados com a extensao.
               Parametros: CliqueExtender, CGraph, Clique, peso do clique, forma de extensao.
               Forma de extensao pode ser: CLQEM_RANDOM, CLQEM_MAX_DEGREE, CLQEM_PRIORITY_GREEDY ou CLQEM_EXACT
               Para extensao do tipo CLQEM_PRIORITY_GREEDY ou CLQEM_EXACT eh necessario informar um custo para cada var
            */
            clqe_set_costs(clqe, costs, cgraph_size(cg));
            int cliqueW = sizecliques;
            int status = clqe_extend(clqe, cg, &clique, cliqueW, CLQEM_EXACT);

            printf("%d - lifting module found %d new cliques:\n", k, status);
            saida << k << " " << status << "\n";

            /* novos cliques sao obtidos por meio da consulta seguinte */
            const CliqueSet *clqSet = clqe_get_cliques(clqe);
            for (int i = 0; i < clq_set_number_of_cliques(clqSet); i++) {
                printf("   %d ", clq_set_number_of_cliques(clqSet));
                saida << clq_set_clique_size(clqSet, i) << " ";
                for (int j = 0; j < clq_set_clique_size(clqSet, i); j++) {
                    saida << clq_set_clique_elements(clqSet, i)[j] << " ";
                    printf("%d ", clq_set_clique_elements(clqSet, i)[j]);
                }
                saida << "\n";
                printf("\n");
            }
            free(cliques);
            vint_set_clean(&clique);
            clqe_free(&clqe);
        }
    }

    file.close();
    saida.close();
    //vint_set_clean(&clique);
    //clqe_free(&clqe);
    cgraph_free(&cg);
    lp_free(&lp);

    return 0;
}