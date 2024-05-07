// g++ -std=c++17 -I./include -o main src/main.cpp


#include "matrix.h"

#include <vector>
#include <cmath>
#include <stdexcept>
#include <limits>
#include <algorithm>
#include <numeric>
#include <iostream>












using namespace std;
int main() {
    // Define a simple symmetric matrix
    Matrix<double> symMatrix(3, 3);
    symMatrix(0, 0) = 2;
    symMatrix(0, 1) = -1;
    symMatrix(0, 2) = 0;
    symMatrix(1, 0) = -1;
    symMatrix(1, 1) = 2;
    symMatrix(1, 2) = -1;
    symMatrix(2, 0) = 0;
    symMatrix(2, 1) = -1;
    symMatrix(2, 2) = 2;

    // Prepare vectors for eigenvalues and matrices for eigenvectors
    Vector<double> eigenvalues(3);  // Initialize with default size and values
    Matrix<double> eigenvectors(3, 3);

    try {
        // Compute eigenvalues and eigenvectors
        eigen_sym(eigenvalues, eigenvectors, symMatrix);

        // Output the results
        cout << "Eigenvalues:\n";
        for (size_t i = 0; i < eigenvalues.size(); i++) {
            cout << eigenvalues(i) << "\n";
        }

        cout << "\nEigenvectors:\n";
        for (size_t i = 0; i < eigenvectors.nRows(); i++) {
            for (size_t j = 0; j < eigenvectors.nCols(); j++) {
                cout << eigenvectors(i, j) << " ";
            }
            cout << "\n";
        }
    } catch (const std::exception& e) {
        cerr << "Error: " << e.what() << endl;
    }

    return 0;
}