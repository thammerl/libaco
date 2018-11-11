all: build build/liblocalsearch.a

build:
	mkdir -p build

build/liblocalsearch.a: build/localsearch.o
	ar rcs build/liblocalsearch.a build/*.o

build/localsearch.o: src/localsearch.cpp include/liblocalsearch/localsearch.h
	g++ -Wall -ansi -I ./include -c src/localsearch.cpp -o build/localsearch.o -O3

doc: ./src/ ./include/liblocalsearch/
	doxygen

clean: .
	rm -Rf build
