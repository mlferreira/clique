CXXFLAGS = -std=c++11 -ofast -DCBC -fopenmp 

LIBS = -L/opt/gurobi750/linux64//lib/../lib -L/usr/lib/gcc/x86_64-pc-linux-gnu/7.2.0 -L/usr/lib/gcc/x86_64-pc-linux-gnu/7.2.0/../../../../lib -L/lib/../lib -L/usr/lib/../lib -L/opt/ibm/ILOG/CPLEX_Studio127//cplex/lib/x86-64_linux/static_pic -L/opt/gurobi750/linux64//lib -L/usr/lib/gcc/x86_64-pc-linux-gnu/7.2.0/../../.. -L/opt/gurobi750/linux64//lib/../lib -L/usr/lib/gcc/x86_64-pc-linux-gnu/7.2.0 -L/usr/lib/gcc/x86_64-pc-linux-gnu/7.2.0/../../../../lib -L/lib/../lib -L/usr/lib/../lib -L/opt/ibm/ILOG/CPLEX_Studio127//cplex/lib/x86-64_linux/static_pic -L/opt/gurobi750/linux64//lib -L/usr/lib/gcc/x86_64-pc-linux-gnu/7.2.0/../../.. -lgomp -lCbcSolver -lCbc -lpthread -lrt -lreadline -lncurses -lCgl -lOsiClp -lClpSolver -lClp -lreadline -lncurses -lOsi -lCoinUtils -lreadline -lncurses -lbz2 -lz -lm -lcoinlapack -lz -lpthread -lrt -lgfortran -lm -lquadmath -lcoinblas -lz -lpthread -lrt -lgfortran -lm -lquadmath

all: BKGraph.o BKVertex.o bron_kerbosch.o build_cgraph.o cgraph.o clique.o clique_extender.o lp.o memory.o node_heap.o strutils.o vectormgm.o vint_set.o  main.cpp
	g++ -o main main.cpp BKGraph.o BKVertex.o bron_kerbosch.o build_cgraph.o cgraph.o clique.o clique_extender.o lp.o memory.o node_heap.o strutils.o vectormgm.o vint_set.o $(LIBS) $(CXXFLAGS)
	
BKGraph.o: BKGraph.cpp BKGraph.hpp
	g++ -o BKGraph.o -c BKGraph.cpp

BKVertex.o: BKVertex.cpp BKVertex.hpp
	g++ -o BKVertex.o -c BKVertex.cpp

bron_kerbosch.o: bron_kerbosch.cpp bron_kerbosch.h
	g++ -o bron_kerbosch.o -c bron_kerbosch.cpp

build_cgraph.o: build_cgraph.cpp build_cgraph.h
	g++ -o build_cgraph.o -c build_cgraph.cpp

cgraph.o: cgraph.c cgraph.h
	g++ -o cgraph.o -c cgraph.c

clique.o: clique.c clique.h
	g++ -o clique.o -c clique.c

clique_extender.o: clique_extender.c clique_extender.h
	g++ -o clique_extender.o -c clique_extender.c

lp.o: lp.cpp lp.h
	g++ -o lp.o -c lp.c

memory.o: memory.c memory.h
	g++ -o memory.o -c memory.c

node_heap.o: node_heap.c node_heap.h
	g++ -o node_heap.o -c node_heap.c

strutils.o: strutils.c strutils.h
	g++ -o strutils.o -c strutils.c

vectormgm.o: vectormgm.c vectormgm.h
	g++ -o vectormgm.o -c vectormgm.c

vint_set.o: vint_set.c vint_set.h
	g++ -o vint_set.o -c vint_set.c


clean:
	rm -f *.o main
	
clique:
	./cmake-build-debug/clique cmake-build-debug/instancias/air03.mps
	./cmake-build-debug/clique cmake-build-debug/instancias/air04.mps
	./cmake-build-debug/clique cmake-build-debug/instancias/air05.mps
	./cmake-build-debug/clique cmake-build-debug/instancias/clq1.lp
	./cmake-build-debug/clique cmake-build-debug/instancias/eil33-2.mps
	./cmake-build-debug/clique cmake-build-debug/instancias/eilA76.mps
	./cmake-build-debug/clique cmake-build-debug/instancias/eilA101.2.mps
	./cmake-build-debug/clique cmake-build-debug/instancias/eilB76.mps
	./cmake-build-debug/clique cmake-build-debug/instancias/eilB101.2.mps
	./cmake-build-debug/clique cmake-build-debug/instancias/eilB101.mps
	./cmake-build-debug/clique cmake-build-debug/instancias/eilC76.2.mps
	./cmake-build-debug/clique cmake-build-debug/instancias/eilC76.mps
	./cmake-build-debug/clique cmake-build-debug/instancias/eilD76.2.mps
	./cmake-build-debug/clique cmake-build-debug/instancias/eilD76.mps
	./cmake-build-debug/clique cmake-build-debug/instancias/stdc0025.mps
	./cmake-build-debug/clique cmake-build-debug/instancias/stdc1220.mps
	./cmake-build-debug/clique cmake-build-debug/instancias/stdc3282.mps
	./cmake-build-debug/clique cmake-build-debug/instancias/stdc6262p.mps
	./cmake-build-debug/clique cmake-build-debug/instancias/stdcd5113.mps
	./cmake-build-debug/clique cmake-build-debug/instancias/trd445c.mps
	./cmake-build-debug/clique cmake-build-debug/instancias/trdnc.mps
	./cmake-build-debug/clique cmake-build-debug/instancias/trdnc18.mps
	./cmake-build-debug/clique cmake-build-debug/instancias/trdta0010.mps
	./cmake-build-debug/clique cmake-build-debug/instancias/z26.lp
	./cmake-build-debug/clique cmake-build-debug/instancias/long_hidden01.mps
	./cmake-build-debug/clique cmake-build-debug/instancias/sprint01_j.mps
	./cmake-build-debug/clique cmake-build-debug/instancias/sprint_hidden01_j.mps
	./cmake-build-debug/clique cmake-build-debug/instancias/sprint_hint01_j.mps
	./cmake-build-debug/clique cmake-build-debug/instancias/trdtatl9220.mps
