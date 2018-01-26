#include <cstdio>
#include <iostream>
#include <cstdlib>
#include <algorithm>
#include <fstream>
#include <ctime>
#include <vector>

using namespace std;
extern "C" {
#include "lp.h"
#include "vint_set.h"
#include "clique_separation.h"
#include "clique_extender.h"
#include "build_cgraph.h"
#include "clique.h"
}

#define LIMITE_GREEDY 100000
#define LIMITE_ROWSADDED 500


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
    int n=0, rowsAdded=0, rowsAddedTotal=0, rowsRemoved=0, rowsRemovedTotal=0, contadorNome=0;
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

    int* rowsToRemove = new int[nRes];
    int* idx = new int[nVar];
    int* idxCompl = new int[nVar];
    double* coef = new double[nVar];
    double* coefCompl = new double[nVar];
    int* vetor = new int[nVar];

    vector<int> rowsAdded_start;
    vector<int> rowsAdded_idx;
    vector<double> rowsAdded_coef;
    vector<char> rowsAdded_sense;
    vector<double> rowsAdded_rhs;
    vector<const char*> rowsAdded_names;



    int resToCheck[nRes];
    for(int i=0; i<nRes; i++) {
        n = lp_row(mip, i, idx, coef);
        resToCheck[i] = check_clique(mip, i, idx, coef, n);
    }
    fill(coef, coef+nVar, 1.0);

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

    for(int k=0; k<nRes; k++) {


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
            int status=0;
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

                saida << "Restricao " << k << " (" << lp_row_name(mip, k, nomeRes) << "), " << status << " novos cliques\n";
                saida << "    Original" << " (" << lp_row(mip, k, idx, coef) << "): ";
                for (int j = 0; j < n; j++) {
                    saida << idx[j] << " (" << lp_col_name(mip, idx[j], nomeRes) << ") ";
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

                fill(vetor, vetor+nVar, 0);

                for (int i=0; i<clq_set_number_of_cliques(clqSet); i++) {


                    // nao tem complemento
                    if(clq_set_clique_elements(clqSet, i)[clq_set_clique_size(clqSet, i)-1] < nVar) {


#ifdef DEBUG
                        //printf("   %d (%d): ", i, clq_set_clique_size(clqSet, i));
                        saida << "    " << i << " (" << clq_set_clique_size(clqSet, i) << "): ";

                        tempoRes << clq_set_clique_size(clqSet, i)-n << " ";
#endif

                        for (int j=0; j<clq_set_clique_size(clqSet, i); j++) {
#ifdef DEBUG
                            saida << clq_set_clique_elements(clqSet, i)[j] << " (" << clq_set_clique_elements(clqSet, i)[j]-nVar
                                  << " " << lp_col_name(mip, clq_set_clique_elements(clqSet, i)[j]-nVar, nomeRes) << ") ";
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
                                    rowsRemoved++; rowsRemovedTotal++;

#ifdef DEBUG
                                    //cout << "        domina a restricao " << lp_row_name(mip, j, name) << endl;
                                    saida << "        domina a restricao " << lp_row_name(mip, j, name) << "\n";
#endif
                                }

                            }

                        }
#ifdef DEBUG
                        saida << "\n";
#endif

                        if (rowsRemoved > 0) {
                            rowsAdded_names.push_back(("clq00" + to_string(contadorNome)).c_str());
                            contadorNome++;
                            rowsAdded++; rowsAddedTotal++;
                            rowsAdded_start.push_back(rowsAdded_coef.size());
                            rowsAdded_sense.push_back('L');
                            rowsAdded_rhs.push_back(1);
                            for(int j=0; j<clq_set_clique_size(clqSet, i); j++) {
                                rowsAdded_idx.push_back(clq_set_clique_elements(clqSet, i)[j]);
                                rowsAdded_coef.push_back(1);
                            }

                            if (resToCheck[k] == 2) {
                                rowsAdded_names.push_back(("clq00" + to_string(contadorNome)).c_str());
                                contadorNome++;
                                rowsAdded++; rowsAddedTotal++;
                                rowsAdded_start.push_back(rowsAdded_coef.size());
                                rowsAdded_sense.push_back('G');
                                rowsAdded_rhs.push_back(1);
                                for(int j=0; j<clq_set_clique_size(clqSet, i); j++) {
                                    rowsAdded_idx.push_back(clq_set_clique_elements(clqSet, i)[j]);
                                    rowsAdded_coef.push_back(1);
                                }
                            }
                            resToCheck[k] = 0;
                        }

                    }

                        // tem complemento
                    else {

                        fill(vetor, vetor+nVar, 0);
                        fill(coefCompl, coefCompl+nVar, 1.0);
#ifdef DEBUG
                        //printf("   %d (%d): ", i, clq_set_clique_size(clqSet, i));
                        saida << "    " << i << " (" << clq_set_clique_size(clqSet, i) << "): ";

                        tempoRes << clq_set_clique_size(clqSet, i)-n << " ";
#endif

                        for (int j=0; j<clq_set_clique_size(clqSet, i); j++) {



                            if(clq_set_clique_elements(clqSet, i)[j] >= nVar && vetor[clq_set_clique_elements(clqSet, i)[j]-nVar] == 0){
#ifdef DEBUG
                                saida << clq_set_clique_elements(clqSet, i)[j] << " (" << clq_set_clique_elements(clqSet, i)[j]-nVar
                                      << " " << lp_col_name(mip, clq_set_clique_elements(clqSet, i)[j]-nVar, nomeRes) << ") ";
#endif
                                idxCompl[j] = clq_set_clique_elements(clqSet, i)[j]-nVar;
                                coefCompl[j] = -1.0;
                            }
                            else if (vetor[clq_set_clique_elements(clqSet, i)[j]] == 0){
#ifdef DEBUG
                                saida << clq_set_clique_elements(clqSet, i)[j] << " ("
                                      << lp_col_name(mip, clq_set_clique_elements(clqSet, i)[j], nomeRes) << ") ";
#endif
                                //monta vetor de 0/1
                                vetor[clq_set_clique_elements(clqSet, i)[j]] = 1;
                                idxCompl[j] = clq_set_clique_elements(clqSet, i)[j];
                            }

                        }

#ifdef DEBUG
                        saida << "\n";
                        //printf("\n");
#endif

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
                                    rowsRemoved++; rowsRemovedTotal++;

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
                            resToCheck[k] = 0;
                            rowsAdded_names.push_back(("clq00" + to_string(contadorNome)).c_str());
                            contadorNome++;
                            rowsAdded++; rowsAddedTotal++;
                            rowsAdded_start.push_back(rowsAdded_coef.size());
                            rowsAdded_sense.push_back('L');
                            rowsAdded_rhs.push_back(0);
                            for(int j=0; j<clq_set_clique_size(clqSet, i); j++) {
                                rowsAdded_idx.push_back(idxCompl[j]);
                                rowsAdded_coef.push_back(coefCompl[j]);
                            }
                        }


                    }

                }
                endRestr=clock();
                timeRestr += (endRestr-startRest)/double(CLOCKS_PER_SEC)*1000;

            } //if (status > 0)

            cout.flush(); saida.flush();

        } //if (resToCheck[k] > 0)
        
        if(rowsAdded_idx.size() > LIMITE_ROWSADDED) {
            //cout << rowsAdded_idx.size() << endl; cout.flush();
            rowsAdded_start.push_back(rowsAdded_coef.size());
            lp_add_rows(mip, rowsAdded, &rowsAdded_start[0], &rowsAdded_idx[0], &rowsAdded_coef[0], &rowsAdded_sense[0],
                        &rowsAdded_rhs[0], &rowsAdded_names[0]);
            rowsAdded = 0;
            rowsAdded_start.clear();
            rowsAdded_idx.clear();
            rowsAdded_coef.clear();
            rowsAdded_sense.clear();
            rowsAdded_rhs.clear();
            rowsAdded_names.clear();
        }

    } //for (k<nRes)

    if(rowsRemoved > 0) {
        lp_remove_rows(mip, rowsRemoved, rowsToRemove);
    }
    if(rowsAdded > 0) {
        rowsAdded_start.push_back(rowsAdded_coef.size());
        lp_add_rows(mip, rowsAdded, &rowsAdded_start[0], &rowsAdded_idx[0], &rowsAdded_coef[0], &rowsAdded_sense[0],
                    &rowsAdded_rhs[0], &rowsAdded_names[0]);
    }


#ifdef DEBUG
    int endFor=clock();
#endif



#ifdef DEBUG



    resMax << indiceMaior << " " << resMaior << " " << nVarMaior << " " << tempoMaior << " " << nExpansoes;

    saida << "\n\nRestricoes removidas:   " << rowsRemovedTotal
          << "\nRestricoes adicionadas: " << rowsAddedTotal << "\n"
          << "\nTempo execucao: " << (endFor-startFor)/double(CLOCKS_PER_SEC)*1000 << "ms"
          << " (Total: " << (endFor-start)/double(CLOCKS_PER_SEC)*1000 << "ms)\n";

    cout << "\n\nRestricoes removidas:   " << rowsRemovedTotal
         << "\nRestricoes adicionadas: " << rowsAddedTotal << endl << endl
         << "\nTempo execucao: " << (endFor-startFor)/double(CLOCKS_PER_SEC)*1000 << "ms "
         << " (Total: " << (endFor-start)/double(CLOCKS_PER_SEC)*1000 << "ms)" << endl; cout.flush();



    fstream tempo("0000_tempos.txt", fstream::app);
    tempo   << "grafo: " << (endGraph-startGraph)/double(CLOCKS_PER_SEC)*1000 << " "
            << "preprocessamento: " << (endPreprocess-startPreprocess)/double(CLOCKS_PER_SEC)*1000 << " "
            << "clique: " << timeClique << " "
            << "restricoes: " << timeRestr << " "
            << "for: " << (endFor-startFor)/double(CLOCKS_PER_SEC)*1000 << " ";


    fstream resultado("0000_resultados.txt", fstream::app);
    resultado << file << " "
              << rowsRemovedTotal << " "
              << rowsAddedTotal << " "
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
    delete[] vetor;
    delete[] idx;
    delete[] coef;
    delete[] rowsToRemove;
    delete[] idxCompl;
    delete[] coefCompl;

}


int check_dominance(int* var, int* res, int sizeRes) {

    int d = 1;

    for(int i=0; i<sizeRes; i++) {
        d *= var[res[i]];
    }

    return d;

}

int check_clique(LinearProgram* mip, int i, int *idx, double *coef, int n) {

    if(n==1) {
        return 0;
    }

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

    return 0;

}
