all: matrices

parsing.o:
	g++ parsing.cpp -std=c++11 -c -o parsing.o

matrices: parsing.o
	g++ main.cpp parsing.o -std=c++11 -o matrices

clean:
	rm matrices
	rm *.o
