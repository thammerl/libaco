all: build lib bin bin/acotsp

build:
	mkdir -p build

lib:
	mkdir -p lib

bin:
	mkdir -p bin

bin/acotsp: libaco liblocalsearch build/acotsp.o build/tsp.o
	cp ../libaco/build/libaco.a lib/
	cp ../liblocalsearch/build/liblocalsearch.a lib/
	g++ build/acotsp.o build/tsp.o -L./lib/ -laco -o bin/acotsp

build/acotsp.o: src/acotsp.cpp include/acotsp/tsp.h include/libaco/ants.h include/libaco/graph.h include/libaco/util.h
	g++ -Wall -ansi -I include/ -c src/acotsp.cpp -o build/acotsp.o

build/tsp.o: src/tsp.cpp include/libaco/ants.h include/libaco/graph.h include/libaco/util.h include/liblocalsearch/localsearch.h
	g++ -Wall -ansi -I include/ -c src/tsp.cpp -o build/tsp.o

include/libaco/%.h: ../libaco/include/libaco/%.h
	cp ../libaco/include/libaco/$*.h include/libaco/

include/liblocalsearch/%.h: ../liblocalsearch/include/liblocalsearch/%.h
	cp ../liblocalsearch/include/liblocalsearch/$*.h include/liblocalsearch/

libaco:
	make -C ../libaco -f ../libaco/Makefile

liblocalsearch:
	make -C ../liblocalsearch -f ../liblocalsearch/Makefile

clean: .
	rm include/libaco/*.h
	rm -Rf build lib bin
