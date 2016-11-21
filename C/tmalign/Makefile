all: TMalign

TMalign: TMalign.cpp basic_fun.h  Kabsch.h  NW.h TMalign.h

	g++ -static -O3 -ffast-math -lm -o TMalign TMalign.cpp 

so: TMalign.cpp basic_fun.h  Kabsch.h  NW.h TMalign.h
	g++ -c TMalign.cpp -o TMalign.o -g -fPIC -ffast-math -O3 -lm
	g++ -fPIC --shared TMalign.o -o tmalign.so 

clean:
	rm -f TMalign
	rm -f TMalign.o
	rm -f tmalign.so
