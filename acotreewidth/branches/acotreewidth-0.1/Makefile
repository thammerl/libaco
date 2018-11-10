all: acotreewidth.o ants.o decomp.o graph.o ls.o util.o
	g++ -o acotreewidth acotreewidth.o ants.o decomp.o graph.o ls.o util.o

acotreewidth.o: acotreewidth.cpp decomp.h ants.h util.h
	g++ -c acotreewidth.cpp -I ./include

ants.o: ants.cpp ants.h graph.h util.h
	g++ -c ants.cpp

graph.o: graph.cpp graph.h
	g++ -c graph.cpp

decomp.o: decomp.cpp ants.h decomp.h graph.h util.h
	g++ -c decomp.cpp

ls.o: ls.cpp ls.h
	g++ -c ls.cpp

util.o: util.cpp util.h
	g++ -c util.cpp

clean: .
	rm *.o acotreewidth
