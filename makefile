OBJ=slab.o push.o
LIB=
CPP=g++
OPT=

slab:	slab.o push.o
	$(CPP) $(OPT) -o slab $(OBJ) $(LIB)

slab.o:		slab.cpp slab.hpp MultiArrays.hpp
	$(CPP) $(OPT) -c slab.cpp

push.o:		push.cpp slab.hpp MultiArrays.hpp
	$(CPP) $(OPT) -c push.cpp

clean:
	rm *.csv *.o slab
