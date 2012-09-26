HEADER := inc/VectorT.hh inc/VectorT.h inc/VectorT_inc.h inc/VectorT_inc.hh inc/config.h inc/Timer.hh

all: apmesh cpmesh

apmesh: apmesh.cpp
	g++ -GL -GLU apmesh.cpp -o apmesh

cpmesh: cpmesh.cpp
	g++ -GL -GLU cpmesh.cpp -o cpmesh

clean:
	-rm apmesh
	-rm cpmesh
