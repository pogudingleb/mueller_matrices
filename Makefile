EIGEN_PATH = ./eigen
CFLAGS = -std=c++11 -O2 -I $(EIGEN_PATH)

all: matrices

parsing.o:
	g++ parsing.cpp $(CFLAGS) -c -o parsing.o

matrices: parsing.o
	g++ main.cpp parsing.o $(CFLAGS) -o matrices

clean:
	rm matrices
	rm *.o
