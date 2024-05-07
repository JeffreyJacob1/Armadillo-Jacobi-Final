# Compiler
CPP = g++

# Compiler flags
CPPFLAGS = -O3 -std=c++17

main:
	g++ -std=c++17 -I./include -o main src/main.cpp

Jacobi-Test:
	g++ -std=c++17 -O3 -o Jacobi-Test Jacobi-Test.cpp -larmadillo

all: Jacobi-Test main

