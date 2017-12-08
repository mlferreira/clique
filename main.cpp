#include <cstdio>
#include <iostream>
#include <cstdlib>
#include <algorithm>
#include <fstream>
#include <ctime>
using namespace std;
extern "C" {
#include "lp.h"
#include "vint_set.h"
#include "clique_extender.h"
#include "build_cgraph.h"
}


void clq_preproccess(LinearProgram *mip, string file = "resultado.lp");

int check_dominance(int* var, int* res, int sizeRes);


int main(int argc, const char **argv) {

    string arqSaida = argv[1];
    while( arqSaida.back() != '.') {
        arqSaida.pop_back();
    }
    arqSaida.pop_back();
    arqSaida += "_new";



    int startLeitura=clock();

    LinearProgram *lp = lp_create();
    lp_read(lp, argv[1]);

    int endLeitura=clock();


    fstream tempo("0000_tempos.txt", fstream::app);
    tempo << "\n" << argv[1] << " "
          << "leitura: " << (endLeitura-startLeitura)/double(CLOCKS_PER_SEC)*1000 << " ";
    tempo.close();


    clq_preproccess(lp, argv[1]);


    int end=clock();

    tempo.open("0000_tempos.txt", fstream::app);
    tempo << "total: " << (end-startLeitura)/double(CLOCKS_PER_SEC)*1000 << " ";
    tempo.close();


    lp_write_lp(lp, arqSaida.c_str());


    lp_free(&lp);

    return 0;
}


void clq_preproccess(LinearProgram *mip, string file) {
    int n=0, rowsRemoved=0, rowsAdded=0, contadorNome=0;
    int startRest, endRestr, startClique, endClique, start=clock();
    double timeClique=0.0, timeRestr=0.0;
    char* name = new char[10];

    string arqPrint = file;
    while( arqPrint.back() != '.') {
        arqPrint.pop_back();
    }
    arqPrint.pop_back();
    arqPrint += "_print.txt";



    int startGraph=clock();

    CGraph *cg = build_cgraph_lp(mip);
    const int nVertices = cgraph_size(cg);

    int endGraph=clock();



    int startPreprocess=clock();

    int costs[nVertices];
    fill(costs, costs + nVertices, 1); /* preenchendo o custo de cada variavel com valor 1 */

    int nRes = lp_rows(mip);
    int nVar = lp_cols(mip);

    int* rowsToRemove = new int[nRes];
    int* idx = new int[nVar];
    int* vetor = new int[nVar];
    double* coef = new double[nVar];
    fill(coef, coef+nVar, 1.0);
    int resToCheck[nRes];
    for(int i=0; i<nRes; i++) {
        if(lp_sense(mip, i) == 'L' && lp_rhs(mip, i) == 1) {
            resToCheck[i] = 1;
        } else if(lp_sense(mip, i) == 'E' && lp_rhs(mip, i) == 1) {
            resToCheck[i] = 2;
        } else {
            resToCheck[i] = 0;
        }
    }
    int endPreprocess=clock();



    ofstream saida(arqPrint);


    IntSet clq, clique, res;
    vint_set_init( &clq );
    vint_set_init( &clique );
    vint_set_init( &res );


    int startFor=clock();
    for(int k = 0; k < nRes; k++) {

        if(resToCheck[k] > 0) {

            n = lp_row(mip, k, idx, coef);



            startClique=clock();

            vint_set_clear(&clique);
            vint_set_add(&clique, idx, n);

            CliqueExtender *clqe = clqe_create();

            /* Funcao retorna 0 se o clique nao foi estendido ou o numero de cliques gerados com a extensao.
               Parametros: CliqueExtender, CGraph, Clique, peso do clique, forma de extensao.
               Forma de extensao pode ser: CLQEM_RANDOM, CLQEM_MAX_DEGREE, CLQEM_PRIORITY_GREEDY ou CLQEM_EXACT
               Para extensao do tipo CLQEM_PRIORITY_GREEDY ou CLQEM_EXACT eh necessario informar um custo para cada var
            */
            clqe_set_costs(clqe, costs, cgraph_size(cg));
            int cliqueW = n;
            int status = clqe_extend(clqe, cg, &clique, cliqueW, CLQEM_EXACT);

            endClique=clock();
            timeClique += (endClique-startClique)/double(CLOCKS_PER_SEC)*1000;


            if (status > 0) {

                printf("\n%s - lifting module found %d new cliques:\n", lp_row_name(mip, k, name), status);
                saida << "Restricao " << k << ", " << status << " novos cliques\n";

                printf("   Original (%d): ", lp_row(mip, k, idx, coef));
                saida << "    Original" << " (" << lp_row(mip, k, idx, coef) << "): ";
                for (int j = 0; j < n; j++) {
                    saida << idx[j] << " ";
                    printf("%d ", idx[j]);
                }
                saida << "\n";
                printf("\n");


                /* novos cliques sao obtidos por meio da consulta seguinte */
                const CliqueSet *clqSet = clqe_get_cliques(clqe);


                startRest=clock();
                for (int i=0; i<clq_set_number_of_cliques(clqSet); i++) {

                    // nao tem complemento
                    if(clq_set_clique_elements(clqSet, i)[clq_set_clique_size(clqSet, i)-1] < nVar) {


                        fill(vetor, vetor+nVar, 0);

                        printf("   %d (%d): ", i, clq_set_clique_size(clqSet, i));
                        saida << "    " << i << " (" << clq_set_clique_size(clqSet, i) << "): ";

                        for (int j=0; j<clq_set_clique_size(clqSet, i); j++) {

                            saida << clq_set_clique_elements(clqSet, i)[j] << " ";
                            printf("%d ", clq_set_clique_elements(clqSet, i)[j]);

                            //monta vetor de 0/1
                            vetor[clq_set_clique_elements(clqSet, i)[j]] = 1;
                        }

                        saida << "\n";
                        printf("\n");


                        // transforma o clique estendido em IntSet
                        //vint_set_clear( &clq );
                        //vint_set_add(&clq, clq_set_clique_elements(clqSet, i), clq_set_clique_size(clqSet, i));

                        for (int j=0; j<nRes; j++) {

                            if (resToCheck[j] > 0) {

                                //pega a j-esima restricao
                                vint_set_clear(&res);
                                int n2 = lp_row(mip, j, idx, coef);
                                //vint_set_add(&res, idx, n2);

                                //se o clique domina a restricao
                                if (check_dominance(vetor, idx, n2) == 1) {
                                    if(resToCheck[k] == 2){
                                        string nome = "clq00" + to_string(contadorNome);
                                        contadorNome++; rowsAdded++;
                                        n = lp_row(mip, j, idx, coef);
                                        lp_add_row(mip, n, idx, coef, nome.c_str(), 'G', 1);
                                    }
                                    rowsToRemove[rowsRemoved] = j;
                                    resToCheck[j] = 0;
                                    rowsRemoved++;
                                    cout << "        domina a restricao " << lp_row_name(mip, j, name) << endl;
                                    saida << "        domina a restricao " << lp_row_name(mip, j, name) << "\n";
                                    //cout << rowsRemoved << "   " << rowsAdded << endl;
                                }

                            }

                        }
                        saida << "\n";

                        if (rowsRemoved > 0) {
                            string nome = "clq00" + to_string(contadorNome);
                            contadorNome++;
                            rowsAdded++;
                            lp_add_row(mip, clq_set_clique_size(clqSet, i),
                                       const_cast<int *>(clq_set_clique_elements(clqSet, i)), coef, nome.c_str(), 'L',
                                       1);
                        }

                    }

                        // tem complemento
                    else {
//                        string nome = "clq00" + to_string(contadorNome);
//                        contadorNome++;
//                        rowsAdded++;
//
//                        cout << endl << endl << nome << " COMPLEMENTO!!!" << endl << endl; cout.flush();
//
//                        int* elem = const_cast<int*>(clq_set_clique_elements(clqSet, i));
//
//                        for(int j=0; j<clq_set_clique_size(clqSet, i); j++) {
//                            if (elem[j] >= nVar) {
//                                elem[j] = elem[j]-nVar;
//                                coef[j] = -1.0;
//                            }
//                        }
//
//
//                        lp_add_row(mip, clq_set_clique_size(clqSet, i), elem, coef, nome.c_str(), 'L', 1);



                    }
                }
                endRestr=clock();
                timeRestr += (endRestr-startRest)/double(CLOCKS_PER_SEC)*1000;

            }

            cout.flush();

        }

    }


    if (rowsRemoved > 0) {
        lp_remove_rows(mip, rowsRemoved, rowsToRemove);
    }


    int endFor=clock();



    saida << "\n\nRestricoes removidas:   " << rowsRemoved
          << "\nRestricoes adicionadas: " << rowsAdded << "\n"
          << "\nTempo execucao: " << (endFor-startFor)/double(CLOCKS_PER_SEC)*1000 << "ms"
          << " (Total: " << (endFor-start)/double(CLOCKS_PER_SEC)*1000 << "ms)\n";

    cout << "\n\nRestricoes removidas:   " << rowsRemoved
         << "\nRestricoes adicionadas: " << rowsAdded << endl << endl
         << "\nTempo execucao: " << (endFor-startFor)/double(CLOCKS_PER_SEC)*1000 << "ms "
         << " (Total: " << (endFor-start)/double(CLOCKS_PER_SEC)*1000 << "ms)" << endl; cout.flush();



    fstream tempo("0000_tempos.txt", fstream::app);
    tempo //<< file << " "
            //<< "leitura: " << (endLeitura-startLeitura)/double(CLOCKS_PER_SEC)*1000 << " "
            << "grafo: " << (endGraph-startGraph)/double(CLOCKS_PER_SEC)*1000 << " "
            << "preprocessamento: " << (endPreprocess-startPreprocess)/double(CLOCKS_PER_SEC)*1000 << " "
            << "clique: " << timeClique << " "
            << "restricoes: " << timeRestr << " "
            << "for: " << (endFor-startFor)/double(CLOCKS_PER_SEC)*1000 << " "
            ;//<< "\n";


    fstream resultado("0000_resultados.txt", fstream::app);
    resultado << file << " "
              << rowsRemoved << " "
              << rowsAdded << " "
              << "\n";



    resultado.close();
    tempo.close();
    saida.close();

    cgraph_free(&cg);
    vint_set_clean( &clq );
    vint_set_clean( &clique );
    delete[] idx;
    delete[] vetor;
    delete[](coef);
    delete[](rowsToRemove);
    delete[] name;

}


int check_dominance(int* var, int* res, int sizeRes) {

    int d = 1;

    for(int i=0; i<sizeRes; i++) {
        d *= var[res[i]];
    }

    return d;

}