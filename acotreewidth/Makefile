all: build lib bin bin/acotreewidth

build:
	mkdir -p build

lib:
	mkdir -p lib

bin:
	mkdir -p bin

bin/acotreewidth: libaco liblocalsearch build/acotreewidth.o build/decomp.o
	cp ../libaco/build/libaco.a lib/
	cp ../liblocalsearch/build/liblocalsearch.a lib/
	g++ build/acotreewidth.o build/decomp.o -L./lib/ -laco -llocalsearch -o bin/acotreewidth -O3

build/acotreewidth.o: src/acotreewidth.cpp include/acotreewidth/decomp.h include/libaco/ants.h include/libaco/graph.h include/libaco/util.h include/liblocalsearch/localsearch.h
	g++ -Wall -ansi -I include/ -c src/acotreewidth.cpp -o build/acotreewidth.o -O3

build/decomp.o: src/decomp.cpp include/acotreewidth/decomp.h include/libaco/ants.h include/libaco/graph.h include/libaco/util.h include/liblocalsearch/localsearch.h
	g++ -Wall -ansi -I include/ -c src/decomp.cpp -o build/decomp.o -O3

include/libaco/%.h: ../libaco/include/libaco/%.h
	cp ../libaco/include/libaco/$*.h include/libaco/

include/liblocalsearch/%.h: ../liblocalsearch/include/liblocalsearch/%.h
	cp ../liblocalsearch/include/liblocalsearch/$*.h include/liblocalsearch/

libaco: .
	make -C ../libaco -f ../libaco/Makefile

liblocalsearch: .
	make -C ../liblocalsearch -f ../liblocalsearch/Makefile

clean: .
	rm include/libaco/*.h
	rm -Rf build lib bin
