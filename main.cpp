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
#ifdef _OPENMP
#include <omp.h>
#endif

#define LIMITE_GREEDY 1000
//#define PARALELO

void clq_preproccess(LinearProgram *mip, string file = "resultado.lp");

int check_dominance(int* var, int* res, int sizeRes);

int check_clique(LinearProgram* mip, int i, int *idx, double *coef, int n);


int main(int argc, const char **argv) {
    string arqSaida = argv[1];
    while( arqSaida.back() != '.') {
        arqSaida.pop_back();
    }
    arqSaida.pop_back();
    arqSaida += "_new";

#ifdef DEBUG
    int startLeitura=clock();
#endif

    LinearProgram *lp = lp_create();
    lp_read(lp, argv[1]);

#ifdef DEBUG
    int endLeitura=clock();


    fstream tempo("0000_tempos.txt", fstream::app);
    tempo << "\n" << argv[1] << " "
          << "leitura: " << (endLeitura-startLeitura)/double(CLOCKS_PER_SEC)*1000 << " ";
    tempo.close();
#endif

    clq_preproccess(lp, argv[1]);

#ifdef DEBUG
    int end=clock();

    tempo.open("0000_tempos.txt", fstream::app);
    tempo << "total: " << (end-startLeitura)/double(CLOCKS_PER_SEC)*1000 << " ";
    tempo.close();
#endif

    lp_write_lp(lp, arqSaida.c_str());


    lp_free(&lp);

    return 0;
}


void clq_preproccess(LinearProgram *mip, string file) {
    int n=0, rowsAdded=0, contadorNome=0;
#ifdef DEBUG
    int startRest, endRestr, startClique, endClique, start=clock();
    double timeClique=0.0, timeRestr=0.0;
    char* nomeRes = new char[512];
#endif
    char* name = new char[512];

#ifdef DEBUG
    string arqPrint = file;
    while( arqPrint.back() != '.') {
        arqPrint.pop_back();
    }
    arqPrint.pop_back();
    arqPrint += "_print.txt";
    string arqRes = file;
    while( arqRes.back() != '.') {
        arqRes.pop_back();
    }
    arqRes.pop_back();
    arqRes += "_restricoes.txt";
#endif

#ifdef DEBUG
    int startGraph=clock();
#endif

    CGraph *cg = build_cgraph_lp(mip);
    const int nVertices = cgraph_size(cg);

#ifdef DEBUG
    int endGraph=clock();

    int startPreprocess=clock();
#endif

    int costs[nVertices];
    fill(costs, costs + nVertices, 1); /* preenchendo o custo de cada variavel com valor 1 */

    int nRes = lp_rows(mip);
    int nVar = lp_cols(mip);

#ifdef PARALELO
    int** rowsToRemoveT = new int*[omp_get_num_threads()];
    int** idxT = new int*[omp_get_num_threads()];
    double** coefT = new double*[omp_get_num_threads()];
    int** vetorT = new int*[omp_get_num_threads()];
    int* rowsRemovedT = new int[omp_get_num_threads()];
    for(int i=0; i<omp_get_num_threads(); i++) {
        vetorT[i] = new int[nVar];
        rowsToRemoveT[i] = new int[nVar];
        idxT[i] = new int[nVar];
        coefT[i] = new double[nVar];
        fill(coefT[i], coefT[i]+nVar, 1.0);
        rowsRemovedT[i] = 0;
    }
#else
    int* rowsToRemove = new int[nVar];
    int* idx = new int[nVar];
    double* coef = new double[nVar];
    int* vetor = new int[nVar];
    int rowsRemoved = 0;
#endif

    int resToCheck[nRes];
    for(int i=0; i<nRes; i++) {
#ifdef PARALELO
        n = lp_row(mip, i, idxT[0], coefT[0]);
        resToCheck[i] = check_clique(mip, i, idxT[0], coefT[0], n);
    }
    fill(coefT[0], coefT[0]+nVar, 1.0);
#else
        n = lp_row(mip, i, idx, coef);
        resToCheck[i] = check_clique(mip, i, idx, coef, n);
    }
    fill(coef, coef+nVar, 1.0);
#endif

#ifdef DEBUG
    int endPreprocess=clock();
#endif


#ifdef DEBUG
    ofstream saida(arqPrint);

    ofstream tempoRes(arqRes);
    tempoRes << file;

    fstream resMax("0000_restricoes.txt", fstream::app);
    resMax << "\n" << file << " ";
    int indiceMaior, nVarMaior, nExpansoes;
    char *resMaior = new char[512];
    double tempoMaior = 0.0;
#endif


    IntSet clq, clique, res;
    vint_set_init( &clq );
    vint_set_init( &clique );
    vint_set_init( &res );

#ifdef DEBUG
    int startFor=clock();
#endif

#ifdef PARALELO
#pragma omp parallel for
#endif

    for(int k=0; k<nRes; k++) {

#ifdef PARALELO
        int *vetor = vetorT[omp_get_thread_num()];
        int *idx = idxT[omp_get_thread_num()];
        double *coef = coefT[omp_get_thread_num()];
        int *rowsToRemove = rowsToRemoveT[omp_get_thread_num()];
        int *rowsRemoved = rowsToRemoveT[omp_get_thread_num()];
#endif

        if(resToCheck[k] > 0) {

            n = lp_row(mip, k, idx, coef);


#ifdef DEBUG
            tempoRes << "\n" << k << " " << lp_row_name(mip, k, nomeRes) << " " << n << " ";

            startClique=clock();
#endif

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
            int status;
            if(n > LIMITE_GREEDY) {
                status = clqe_extend(clqe, cg, &clique, cliqueW, CLQEM_PRIORITY_GREEDY);
            }
            else {
                status = clqe_extend(clqe, cg, &clique, cliqueW, CLQEM_EXACT);
            }



#ifdef DEBUG
            endClique=clock();
            timeClique += (endClique-startClique)/double(CLOCKS_PER_SEC)*1000;

            tempoRes << (endClique-startClique)/double(CLOCKS_PER_SEC)*1000 << " " << status << " ";

            if (tempoMaior < (endClique-startClique)/double(CLOCKS_PER_SEC)*1000) {
                indiceMaior = k;
                nVarMaior = n;
                nExpansoes = status;
                resMaior = nomeRes;
                tempoMaior = (endClique-startClique)/double(CLOCKS_PER_SEC)*1000;
            }
#endif



            if (status > 0) {

#ifdef DEBUG
                //printf("\n%s - lifting module found %d new cliques:\n", lp_row_name(mip, k, name), status);
                saida << "Restricao " << k << ", " << status << " novos cliques\n";
                //printf("   Original (%d): ", lp_row(mip, k, idx, coef));
                saida << "    Original" << " (" << lp_row(mip, k, idx, coef) << "): ";
                for (int j = 0; j < n; j++) {
                    saida << idx[j] << " ";
                    //printf("%d ", idx[j]);
                }
                saida << "\n";
                //printf("\n");
#endif


                /* novos cliques sao obtidos por meio da consulta seguinte */
                const CliqueSet *clqSet = clqe_get_cliques(clqe);

#ifdef DEBUG
                startRest=clock();
#endif
                for (int i=0; i<clq_set_number_of_cliques(clqSet); i++) {



                    // nao tem complemento
                    if(clq_set_clique_elements(clqSet, i)[clq_set_clique_size(clqSet, i)-1] < nVar) {


                        fill(vetor, vetor+nVar, 0);
#ifdef DEBUG
                        //printf("   %d (%d): ", i, clq_set_clique_size(clqSet, i));
                        saida << "    " << i << " (" << clq_set_clique_size(clqSet, i) << "): ";

                        tempoRes << clq_set_clique_size(clqSet, i)-n << " ";
#endif

                        for (int j=0; j<clq_set_clique_size(clqSet, i); j++) {
#ifdef DEBUG
                            saida << clq_set_clique_elements(clqSet, i)[j] << " ";
                            //printf("%d ", clq_set_clique_elements(clqSet, i)[j]);
#endif

                            //monta vetor de 0/1
                            vetor[clq_set_clique_elements(clqSet, i)[j]] = 1;
                        }

                        saida << "\n";
                        //printf("\n");

                        for (int j=0; j<nRes; j++) {

                            if (resToCheck[j] > 0) {

                                //pega a j-esima restricao
                                int n2 = lp_row(mip, j, idx, coef);

                                //se o clique domina a restricao
                                if (check_dominance(vetor, idx, n2) == 1) {


                                    rowsToRemove[rowsRemoved] = j;
                                    if (j!=k) {
                                        resToCheck[j] = 0;
                                    }
                                    rowsRemoved++;

#ifdef DEBUG
                                    //cout << "        domina a restricao " << lp_row_name(mip, j, name) << endl;
                                    saida << "        domina a restricao " << lp_row_name(mip, j, name) << "\n";
                                    //cout << rowsRemoved << "   " << rowsAdded << endl;
#endif
                                }

                            }

                        }
#ifdef DEBUG
                        saida << "\n";
#endif

                        if (rowsRemoved > 0) {
                            string nome = "clq00" + to_string(contadorNome);
                            contadorNome++;
                            rowsAdded++;
                            lp_add_row(mip, clq_set_clique_size(clqSet, i),
                                       const_cast<int *>(clq_set_clique_elements(clqSet, i)), coef, nome.c_str(), 'L',
                                       1);
                            if (resToCheck[k] == 2) {
                                nome = "clq00" + to_string(contadorNome);
                                contadorNome++;
                                rowsAdded++;
                                lp_add_row(mip, lp_row(mip, k, idx, coef), idx, coef, nome.c_str(), 'G', 1);
                            }
                            resToCheck[k] = 0;
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

#ifdef PARALELO
    int rowsRemoved = 0;
    for(int i=0; i<omp_get_num_threads(); i++) {
        lp_remove_rows(mip, rowsRemovedT[i], rowsToRemoveT[i]);
        rowsRemoved += rowsRemovedT[i];
    }
#else
    if(rowsRemoved > 0) {
        lp_remove_rows(mip, rowsRemoved, rowsToRemove);
    }
#endif

#ifdef DEBUG
    int endFor=clock();
#endif



#ifdef DEBUG



    resMax << indiceMaior << " " << resMaior << " " << nVarMaior << " " << tempoMaior << " " << nExpansoes;

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
    resMax.close();
#endif

    cgraph_free(&cg);
    vint_set_clean( &clq );
    vint_set_clean( &clique );
    delete[] name;
#ifdef PARALELO
    for(int i=0; i<omp_get_num_threads(); i++) {
        delete[] vetorT[i];
        delete[] idxT[i];
        delete[] coefT[i];
        delete[] rowsToRemoveT[i];
    }
    delete[] vetorT;
    delete[] idxT;
    delete[] coefT;
    delete[] rowsToRemoveT;
    delete[] rowsRemovedT;
#else
    delete[] vetor;
    delete[] idx;
    delete[] coef;
    delete[] rowsToRemove;
#endif

}


int check_dominance(int* var, int* res, int sizeRes) {

    int d = 1;

    for(int i=0; i<sizeRes; i++) {
        d *= var[res[i]];
    }

    return d;

}

int check_clique(LinearProgram* mip, int i, int *idx, double *coef, int n) {

    int min = 0;

    for(int i=1; i<n; i++) {

       if(coef[i] < coef[min]) {
           min = i;
       }

    }

    if(coef[min]*2 > lp_rhs(mip, i)) {

        if(lp_sense(mip, i) == 'L') {
            return 1;
        }

        else if(lp_sense(mip, i) == 'E') {
            return 2;
        }

    }

    else {
        return 0;
    }

}