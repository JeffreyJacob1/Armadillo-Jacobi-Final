# Armadillo-Jacobi-Final

This repository contains the Final Project code for Chem 279. The final project is a robust replacment for the Armadillo Linear Algebra library, and the key accomplishemt is the creation of a custom Jacobi Diagonalizer. The infrastructure from Armadillo is made from scratch, custom Templated matrix and vector classes, replacing arma::mat and arma::vec classes, which are used in all the custom funtions. In addition to the Jacobi Diagonalizer many more Armaadillo funtions are replaced to buiuld the supporting code, arma::accu, arma mat mul, and is Symmetric are examples. 

The matrix.h file in the Include directory contains all of the funtions and can be used in a main to solve the eigen values and eigen vectors of any symmetric matrix. main.cpp in the src folder is a example of this. Jacobi-Test.cpp is a comprehensive test of armadillo eig_sym vs our Jacobi Diagonalizer, and also validates the correctness.

A pdf of the powerpoint presentation with graphs is incuded

make all to generates the two execs main, Jacobi-Test
