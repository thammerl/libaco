all: build lib bin bin/acotemplate

build:
	mkdir -p build

lib:
	mkdir -p lib

bin:
	mkdir -p bin

bin/acotemplate: libaco liblocalsearch build/acotemplate.o build/template.o
	cp ../libaco/build/libaco.a lib/
	cp ../liblocalsearch/build/liblocalsearch.a lib/
	g++ build/acotemplate.o build/template.o -L./lib/ -laco -o bin/acotemplate

build/acotemplate.o: src/acotemplate.cpp include/acotemplate/template.h include/libaco/ants.h include/libaco/graph.h include/libaco/util.h
	g++ -ansi -I include/ -c src/acotemplate.cpp -o build/acotemplate.o

build/template.o: src/template.cpp include/libaco/ants.h include/libaco/graph.h include/libaco/util.h include/liblocalsearch/localsearch.h
	g++ -ansi -I include/ -c src/template.cpp -o build/template.o

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
