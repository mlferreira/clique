cmake_minimum_required(VERSION 3.8)
project(clique)

set(CMAKE_CXX_STANDARD 11)

set(CMAKE_CXX_FLAGS "-Wall -g -std=c++11 -pg -g3 -D_GLIBCXX_DEBUG -DDEBUG -DCBC -I/usr/include/coin -fopenmp ")

set(CMAKE_EXE_LINKER_FLAGS "-L/opt/gurobi750/linux64//lib/../lib -L/usr/lib/gcc/x86_64-pc-linux-gnu/7.2.0 -L/usr/lib/gcc/x86_64-pc-linux-gnu/7.2.0/../../../../lib -L/lib/../lib -L/usr/lib/../lib -L/opt/ibm/ILOG/CPLEX_Studio127//cplex/lib/x86-64_linux/static_pic -L/opt/gurobi750/linux64//lib -L/usr/lib/gcc/x86_64-pc-linux-gnu/7.2.0/../../.. -L/opt/gurobi750/linux64//lib/../lib -L/usr/lib/gcc/x86_64-pc-linux-gnu/7.2.0 -L/usr/lib/gcc/x86_64-pc-linux-gnu/7.2.0/../../../../lib -L/lib/../lib -L/usr/lib/../lib -L/opt/ibm/ILOG/CPLEX_Studio127//cplex/lib/x86-64_linux/static_pic -L/opt/gurobi750/linux64//lib -L/usr/lib/gcc/x86_64-pc-linux-gnu/7.2.0/../../..")

set(SOURCE_FILES
        BKGraph.cpp
        BKGraph.hpp
        BKVertex.cpp
        BKVertex.hpp
        bron_kerbosch.cpp
        bron_kerbosch.h
        build_cgraph.cpp
        build_cgraph.h
        cgraph.c
        cgraph.h
        clique.c
        clique.h
        clique_extender.c
        clique_extender.h
        lp.cpp
        lp.h
        macros.h
        main.cpp
        memory.c
        memory.h
        node_heap.c
        node_heap.h
        strutils.c
        strutils.h
        vectormgm.c
        vectormgm.h
        vint_set.c
        vint_set.h)

add_executable(clique ${SOURCE_FILES})

target_link_libraries(clique m -lgomp -lCbcSolver -lCbc -lpthread -lrt -lreadline -lncurses -lCgl -lOsiClp -lClpSolver -lClp -lreadline -lncurses -lOsi -lCoinUtils -lreadline -lncurses -lbz2 -lz -lm -lcoinlapack -lz -lpthread -lrt -lgfortran -lm -lquadmath -lcoinblas -lz -lpthread -lrt -lgfortran -lm -lquadmath)
