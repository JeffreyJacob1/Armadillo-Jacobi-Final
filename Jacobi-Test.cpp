#include <vector>
#include <cmath>
#include <stdexcept>
#include <limits>
#include <algorithm>
#include <numeric>
#include <iostream>

using namespace std;

template<typename T>
class Matrix {
protected:
    vector<T> data;
    size_t rows, cols;

public:
    Matrix() : rows(0), cols(0) {}
    Matrix(size_t nRows, size_t nCols, T initVal = T()) : rows(nRows), cols(nCols), data(nRows * nCols, initVal) {}

    T& operator()(size_t i, size_t j) {
        return data[i * cols + j];
    }

    const T& operator()(size_t i, size_t j) const {
        return data[i * cols + j];
    }

    size_t nRows() const { return rows; }
    size_t nCols() const { return cols; }

    void resize(size_t newRows, size_t newCols, T initVal = T()) {
        vector<T> newData(newRows * newCols, initVal);
        size_t minRows = min(rows, newRows);
        size_t minCols = min(cols, newCols);
        for (size_t i = 0; i < minRows; ++i)
            for (size_t j = 0; j < minCols; ++j)
                newData[i * newCols + j] = data[i * cols + j];
        data.swap(newData);
        rows = newRows;
        cols = newCols;
    }

    void reset_size(size_t newRows, size_t newCols, T initVal = T()) {
        resize(newRows, newCols, initVal);
    }

    // Summation of all elements in the matrix
    T sum() const {
        return accumulate(data.begin(), data.end(), T(0));
    }

    // Norm calculation (only makes sense for vector-like matrices)
    T norm() const {
        if (cols != 1 && rows != 1)
            throw logic_error("Norm calculation is only defined for vector-like matrices.");
        T sum = 0;
        for (const auto& value : data) {
            sum += value * value;
        }
        return sqrt(sum);
    }
};

template<typename T>
class Vector : public Matrix<T> {
public:
    Vector() : Matrix<T>(0, 1) {}
    explicit Vector(size_t nRows, T initVal = T()) : Matrix<T>(nRows, 1, initVal) {}

    T& operator()(size_t i) {
        return Matrix<T>::data[i];
    }

    const T& operator()(size_t i) const {
        return Matrix<T>::data[i];
    }

    void resize(size_t newSize, T initVal = T()) {
        Matrix<T>::resize(newSize, 1, initVal);
    }

    void reset_size(size_t newSize, T initVal = T()) {
        resize(newSize, initVal);
    }

    size_t size() const {
        return Matrix<T>::rows;
    }

    // Dot product between two vectors
    T dot(const Vector<T>& other) const {
        if (this->size() != other.size())
            throw invalid_argument("Vectors must be the same size to compute dot product.");
        T result = 0;
        for (size_t i = 0; i < this->size(); i++) {
            result += (*this)(i) * other(i);
        }
        return result;
    }

    // Distance to another vector
    T distanceTo(const Vector<T>& other) const {
        Vector<T> diff = *this - other;
        return diff.norm();
    }

    // Overload subtraction operator to facilitate distance calculation
    Vector<T> operator-(const Vector<T>& other) const {
        if (this->size() != other.size())
            throw invalid_argument("Vectors must be the same size to compute subtraction.");
        Vector<T> result(this->size());
        for (size_t i = 0; i < this->size(); i++) {
            result(i) = (*this)(i) - other(i);
        }
        return result;
    }
};




template<typename T>
Matrix<T> operator*(const Matrix<T>& lhs, const Matrix<T>& rhs) {
    if (lhs.nCols() != rhs.nRows()) {
        throw std::invalid_argument("Matrix dimensions must agree for multiplication.");
    }

    Matrix<T> result(lhs.nRows(), rhs.nCols(), 0);
    for (size_t i = 0; i < lhs.nRows(); ++i) {
        for (size_t k = 0; k < lhs.nCols(); ++k) { // common dimension
            for (size_t j = 0; j < rhs.nCols(); ++j) {
                result(i, j) += lhs(i, k) * rhs(k, j);
            }
        }
    }
    return result;
}

template<typename T>
Matrix<T> operator+(const Matrix<T>& lhs, const Matrix<T>& rhs) {
    if (lhs.nRows() != rhs.nRows() || lhs.nCols() != rhs.nCols()) {
        throw std::invalid_argument("Matrix dimensions must agree for addition.");
    }
    Matrix<T> result(lhs.nRows(), lhs.nCols());
    for (size_t i = 0; i < lhs.nRows(); ++i) {
        for (size_t j = 0; j < lhs.nCols(); ++j) {
            result(i, j) = lhs(i, j) + rhs(i, j);
        }
    }
    return result;
}

template<typename T>
Matrix<T> operator-(const Matrix<T>& lhs, const Matrix<T>& rhs) {
    if (lhs.nRows() != rhs.nRows() || lhs.nCols() != rhs.nCols()) {
        throw std::invalid_argument("Matrix dimensions must agree for subtraction.");
    }
    Matrix<T> result(lhs.nRows(), lhs.nCols());
    for (size_t i = 0; i < lhs.nRows(); ++i) {
        for (size_t j = 0; j < lhs.nCols(); ++j) {
            result(i, j) = lhs(i, j) - rhs(i, j);
        }
    }
    return result;
}

template<typename T>
T sumVector(const Vector<T>& v) {
    T sum = 0;
    for (size_t i = 0; i < v.size(); i++) {
        sum += v(i);
    }
    return sum;
}


template<typename T>
bool isSymmetric(const Matrix<T>& matrix) {
    if (matrix.nRows() != matrix.nCols()) return false;
    size_t n = matrix.nRows();
    for (size_t i = 0; i < n; ++i)
        for (size_t j = 0; j < i; ++j)
            if (matrix(i, j) != matrix(j, i))
                return false;
    return true;
}

template<typename T>
Matrix<double> doubleMatrix(const Matrix<T>& matrix) {
    Matrix<double> result(matrix.nRows(), matrix.nCols());
    for (size_t i = 0; i < matrix.nRows(); ++i)
        for (size_t j = 0; j < matrix.nCols(); ++j)
            result(i, j) = static_cast<double>(matrix(i, j));
    return result;
}



void eigensort(Vector<double>& eigenvalues, Matrix<double>* eigenvectors) {
    vector<size_t> idx(eigenvalues.nRows());
    iota(idx.begin(), idx.end(), 0);

    sort(idx.begin(), idx.end(), [&eigenvalues](size_t i, size_t j) { return eigenvalues(i) < eigenvalues(j); });

    Vector<double> sortedEigenvalues(eigenvalues.nRows());
    Matrix<double> sortedEigenvectors(eigenvalues.nRows(), eigenvalues.nRows());

    for (size_t i = 0; i < idx.size(); ++i) {
        sortedEigenvalues(i) = eigenvalues(idx[i]);
        for (size_t j = 0; j < eigenvalues.nRows(); ++j)
            sortedEigenvectors(j, i) = (*eigenvectors)(j, idx[i]);
    }

    eigenvalues = sortedEigenvalues;
    *eigenvectors = sortedEigenvectors;
}

/**
 * @brief Rotation function needed for Jacobi method.
 * @param matrixA The matrix to rotate.
 * @param sin The sine of the rotation angle.
 * @param tau The tangent of half the rotation angle.
 * @param i The row index of the first element to rotate.
 * @param j The column index of the first element to rotate.
 * @param k The row index of the second element to rotate.
 * @param l The column index of the second element to rotate.
*/
void rotate(Matrix<double>& matrixA, double sin, double tau, int i, int j, int k, int l) {
    // Store the current values of the elements to be rotated
    double g = matrixA(i, j);
    double h = matrixA(k, l);

    // Update the elements using the rotation formula
    matrixA(i, j) = g - sin * (h + g * tau);
    matrixA(k, l) = h + sin * (g - h * tau);
}

/**
 * @brief Jacobi rotation function.
 * @param matrixA The matrix to rotate.
 * @param i The row index of the first element to rotate.
 * @param j The column index of the first element to rotate.
 * @param eigenvalues The vector to store the eigenvalues.
 * @param eigenvectors The matrix to store the eigenvectors.
 * @param z An auxiliary vector used in the rotation process.
 * @param thresh The threshold value for determining if a rotation is needed.
*/
void jacobiRotate(Matrix<double>& matrixA, int i, int j, Vector<double>& eigenvalues, \
                  Matrix<double>& eigenvectors, Vector<double>& z, double thresh) {
    // Calculate the off-diagonal element and the difference between eigenvalues
    double offdiag = 100.0 * abs(matrixA(i, j));
    double h = eigenvalues(j) - eigenvalues(i);

    // Variables for the rotation angle and its trigonometric values
    double theta, t, cos, sin, tau;
    double epsilon = numeric_limits<double>::epsilon();

    // Check if the off-diagonal element is small enough to be considered zero
    if (offdiag <= epsilon * abs(eigenvalues(i)) && offdiag <= epsilon * abs(eigenvalues(j))) {
        matrixA(i, j) = 0.0;
    }
    // Check if the off-diagonal element is larger than the threshold
    else if (abs(matrixA(i, j)) > thresh) {
        h = eigenvalues(j) - eigenvalues(i);

        // Calculate the tangent of the rotation angle
        if (offdiag <= epsilon * abs(h)) {
            t = matrixA(i, j) / h;
        }
        else {
            theta = 0.5 * h / matrixA(i, j);
            t = 1.0 / (abs(theta) + sqrt(1.0 + theta * theta));
            if (theta < 0.0) {
                t = -t;
            }
        }

        // Calculate the cosine and sine of the rotation angle
        cos = 1.0 / sqrt(1.0 + t * t);
        sin = t * cos;
        tau = sin / (1.0 + cos);
        h = t * matrixA(i, j);

        // Update the eigenvalues and the auxiliary vector
        z(i) -= h;
        z(j) += h;
        eigenvalues(i) -= h;
        eigenvalues(j) += h;

        // Set the off-diagonal element to zero after rotation
        matrixA(i, j) = 0.0;

        // Apply the rotation to the remaining elements in the rows and columns
        for (int k = 0; k < i; k++) {
            rotate(matrixA, sin, tau, k, i, k, j);
        }
        for (int k = i + 1; k < j; k++) {
            rotate(matrixA, sin, tau, i, k, k, j);
        }
        for (int k = j + 1; k < matrixA.nRows(); k++) {
            rotate(matrixA, sin, tau, i, k, j, k);
        }

        // Update the eigenvectors matrix with the rotation
        for (int k = 0; k < matrixA.nRows(); k++) {
            rotate(eigenvectors, sin, tau, k, i, k, j);
        }
    }
}

/**
 * @brief Computes the eigenvalues and eigenvectors of a symmetric matrix using the Jacobi method.
 * @tparam T The data type of the input matrix.
 * @param eigenvalues The vector to store the computed eigenvalues.
 * @param eigenvectors The matrix to store the computed eigenvectors.
 * @param matrix The input symmetric matrix.
*/
template <typename T>
void eigen_sym(Vector<double>& eigenvalues, Matrix<double>& eigenvectors, const Matrix<T>& matrix) {
    // Check if the input matrix is symmetric
    if (!isSymmetric(matrix)) {
        throw invalid_argument("Matrix must be symmetric.");
    }

    // Convert the input matrix to a double-precision matrix
    Matrix<double> matrixA = doubleMatrix(matrix);
    const int size = matrix.nRows();

    // Initialize the eigenvalues and eigenvectors
    eigenvalues.reset_size(size);
    eigenvectors.reset_size(size, size);
    int rotations = 0;
    const double epsilon = numeric_limits<double>::epsilon();
    double threshold;

    // Initialize the eigenvectors matrix to the identity matrix
    for (int i = 0; i < size; i++) {
        eigenvectors(i, i) = 1.0;
    }

    // Initialize the auxiliary vectors
    Vector<double> b = Vector<double>(size);
    Vector<double> z = Vector<double>(size);

    // Copy the diagonal elements of the matrix to the eigenvalues vector
    for (int i = 0; i < size; i++) {
        eigenvalues(i) = matrixA(i, i);
        b(i) = eigenvalues(i);
        z(i) = 0.0;
    }

    // Perform iterations until convergence or maximum iterations reached
    for (int iter = 1; iter <= 50; iter++) {
        // Calculate the sum of the absolute values of the off-diagonal elements
        double sum = 0.0;
        for (int i = 0; i < size - 1; i++) {
            for (int j = i + 1; j < size; j++) {
                sum += abs(matrixA(i, j));
            }
        }

        // If the sum is zero, the matrix is already diagonal
        if (sum == 0.0) {
            eigensort(eigenvalues, &eigenvectors);
            return;
        }

        // Set the threshold value based on the iteration number
        if (iter < 4) {
            threshold = 0.2 * sum / (size * size);
        }
        else {
            threshold = 0.0;
        }

        // Perform rotations for each pair of off-diagonal elements
        for (int i = 0; i < size - 1; i++) {
            for (int j = i + 1; j < size; j++) {
                double offdiag = 100.0 * abs(matrixA(i, j));

                // Check if the off-diagonal element is small enough to be considered zero
                if (iter > 4 && offdiag <= epsilon * abs(eigenvalues(i)) && \
                    offdiag <= epsilon * abs(eigenvalues(j))) {
                    matrixA(i, j) = 0.0;
                }
                // Check if the off-diagonal element is larger than the threshold
                else if (abs(matrixA(i, j)) > threshold) {
                    // Perform the Jacobi rotation
                    jacobiRotate(matrixA, i, j, eigenvalues, eigenvectors, z, threshold);
                    rotations++;
                }
            }
        }

        // Update the eigenvalues and reset the auxiliary vector
        for (int i = 0; i < size; i++) {
            b(i) += z(i);
            eigenvalues(i) = b(i);
            z(i) = 0.0;
        }
    }

    // If the maximum iterations are reached without convergence, throw an exception
    throw runtime_error("Too many iterations, Matrix did not converge.");
}

#include <armadillo>
#include <iostream>
#include <vector>
#include <chrono>


using namespace arma;

int main() {
    // Define sizes of matrices to test
    vector<size_t> sizes = {2, 3, 5, 10, 20,50,100};

    // Loop over each size
    for (size_t n : sizes) {
        cout << "Testing matrix size: " << n << "x" << n << endl;

        // Generate symmetric matrix
        Matrix<double> symMatrix(n, n);
        mat armaSymMatrix(n, n);

        // Fill matrix with random values and make it symmetric
        for (size_t i = 0; i < n; i++) {
            for (size_t j = i; j < n; j++) {
                double value = rand() % 10;  // Random values between 0 and 9
                symMatrix(i, j) = symMatrix(j, i) = value;
                armaSymMatrix(i, j) = armaSymMatrix(j, i) = value;
            }
        }

        // Prepare vectors for eigenvalues and matrices for eigenvectors
        Vector<double> eigenvalues(n);
        Matrix<double> eigenvectors(n, n);
        vec armaEigenvalues;
        mat armaEigenvectors;

        // Measure and run custom implementation
        auto start = chrono::high_resolution_clock::now();
        eigen_sym(eigenvalues, eigenvectors, symMatrix);
        auto end = chrono::high_resolution_clock::now();
        chrono::duration<double> elapsed = end - start;
        cout << "Custom implementation took " << elapsed.count() << " seconds." << endl;

        // Measure and run Armadillo's eig_sym
        start = chrono::high_resolution_clock::now();
        eig_sym(armaEigenvalues, armaEigenvectors, armaSymMatrix);
        end = chrono::high_resolution_clock::now();
        elapsed = end - start;
        cout << "Armadillo eig_sym took " << elapsed.count() << " seconds." << endl;

        // Display results from both implementations
        cout << "Custom eigenvalues:\n";
        for (size_t i = 0; i < eigenvalues.size(); i++) {
            cout << eigenvalues(i) << " ";
        }
        cout << "\nArmadillo eigenvalues:\n";
        armaEigenvalues.print();

        cout << "\n---\n";
    }

    return 0;
}

// g++ -std=c++17 -O3 -o test2 jacobi-test2.cpp -larmadillo