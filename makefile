OBJ=slab.o push.o poisson.o
LIB=-lfftw3f -lm
CPP=g++
OPT=
DEB_OPT = -g -O0 -fno-omit-frame-pointer -Wall -Wextra -fsanitize=address
DEB_LIB=-fsanitize=address

slab:	slab.o push.o poisson.o
	$(CPP) $(OPT) -o slab $(OBJ) $(LIB)

slab.o:		slab.cpp slab.hpp MultiArrays.hpp
	$(CPP) $(OPT) -c slab.cpp

push.o:		push.cpp slab.hpp MultiArrays.hpp
	$(CPP) $(OPT) -c push.cpp

poisson.o: poisson.cpp slab.hpp MultiArrays.hpp
	$(CPP) $(OPT) -c poisson.cpp

clean:
	rm *.csv *.o *.dat slab
