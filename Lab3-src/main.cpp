#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <cmath>

#define PRECISION 50
#define SET_STREAM_PRECISION(stream) \
    (stream).setf(ios::fixed); \
    (stream) << setprecision(PRECISION)

using namespace std;

typedef vector<double> Vector;
typedef vector<Vector>  Matrix;

void Print(const Matrix& A);
void Print(const Vector& x);
void Print(const string& filename, const Matrix& A);
Vector VectorSub(const Vector& a, const Vector& b);
Vector MatrixByVector(const Matrix& A, const Vector& x);
double Norm(const Vector& x);
double NormInf(const Vector& x);
double NormInf(const Matrix& A);
Matrix Zeros(const size_t size);
void SetJacobiIterationsParameters(const Matrix& A, const Vector& b, Matrix& C, Vector& g);
double GetMuStopConditionNumber(const Matrix& C);
Vector SeidelIterations(const Matrix& C, const Vector& g, const Vector& x0, const double eps, int& n);
Vector SolveSLAUWithSeidelJacobiMethod(const Matrix& A, const Vector& b, const Vector& x0, const double eps, int& n);

// 7 points
void FixedStartVectorExperiment(const string& filename, const Matrix& A, const Vector& actualSolution, const Vector& x0,
    const double eps0, const size_t steps);

int main() {
    Matrix A = {
        { 1122, 57, 315, 36, 17, 38, 92, 43, 87, 53 },
        { 4, 4331, -14, 31, 53, 323, 45, 32, 16, 75 },
        { 8, 3, 5641, 15, 67, 51, 79, 31, 34, 64 },
        { 18, 1, 239, 4676, 83, 31, 13, 43, 57, 13 },
        { 14, 43, 12, 1, 6653, 57, 86, 43, 24, 10 },
        { 11, 21, 35, 41, 50, 3462, 71, 83, 92, 14 },
        { 91, 83, 74, 61, 57, 49, 3635, 22, 14, 90 },
        { 46, 273, 12, 14, 53, 67, 53, 3513, 61, 87 },
        { 14, 53, 468, 54, 67, 97, 34, 62, 5612, 10 },
        { 65, 13, 53, 16, 774, 80, 45, 89, 314, 1535 }
    };
    for (auto& row : A)
        for (auto& elem : row)
            elem /= 1e2;
    Vector actualSolution = {
        7, 9, 3, 2, 1, 5, 8, 4, 6, 10
    };
    Vector b = MatrixByVector(A, actualSolution);
    Vector x0 = {
        10, 2, 4, 3, 5, 6, 7, 8, 9, 1
    };
    int n = 0;
    double eps = 0.1;
    cout << "A = " << endl;
    Print(A);
    cout << "b = " << endl;
    Print(b);
    Vector x = SolveSLAUWithSeidelJacobiMethod(A, b, x0, eps, n);
    cout << "x = " << endl;
    Print(x);
    cout << "got x with " << n << " iterations for required accuracy of " << eps << " from x0 =" << endl;
    Print(x0);
    cout << "Actual solution x* = " << endl;
    Print(actualSolution);
    Vector error = VectorSub(x, actualSolution);
    Vector residual = VectorSub(MatrixByVector(A, x), b);
    cout << "Infinite norm of x-x* is " << NormInf(error) << endl;
    cout << "(2nd norm of x-x* is " << Norm(error) << ")" << endl;
    cout << "Infinite norm of Ax-b is " << NormInf(residual) << endl;
    cout << "(2nd norm of Ax-b is " << Norm(residual) << ")" << endl;

    Print("../csv/matrix.csv", A);

    FixedStartVectorExperiment("../csv/accuracy_dependencies.csv", A, actualSolution, x0, 0.1, 10);

    return 0;
}

void Print(const Matrix& A) {
    for (auto& row : A) {
        for (auto& elem : row)
            cout << "\t" << elem;
        cout << endl;
    }
}

void Print(const Vector& x) {
    for (auto& elem : x)
        cout << "\t" << elem;
    cout << endl;
}

void Print(const string& filename, const Matrix& A) {
    ofstream file(filename);
    if (!file.is_open()) {
        cerr << "Could not open " << filename << " file!" << endl;
        return;
    }
    SET_STREAM_PRECISION(file);
    for (auto iter_row : A) {
        for (auto iter = iter_row.begin(); iter < iter_row.end(); iter++) {
            file << *iter;
            if (iter + 1 < iter_row.end())
                file << ";";
        }
        file << endl;
    }
    file.close();
}

Vector VectorSub(const Vector& a, const Vector& b) {
    const size_t size = a.size();
    Vector res(size);
    for (size_t i = 0; i < size; i++)
        res[i] = a[i] - b[i];
    return res;
}

Vector MatrixByVector(const Matrix& A, const Vector& x) {
    size_t size = x.size();
    Vector res(size);
    for (size_t i = 0; i < size; i++) {
        res[i] = 0;
        for (size_t j = 0; j < size; j++)
            res[i] += A[i][j] * x[j];
    }
    return res;
}

double Norm(const Vector& x) {
    double normSqr = 0;
    for (auto& elem : x)
        normSqr += elem * elem;
    return sqrt(normSqr);
}

double NormInf(const Vector& x) {
    double max = 0;
    for (auto& elem : x) {
        double elem_abs = abs(elem);
        if (elem_abs > max)
            max = elem_abs;
    }
    return max;
}

double NormInf(const Matrix& A) {
    double maxRowSum = 0;
    for (auto& row : A) {
        double curRowSum = 0;
        for (auto& elem : row)
            curRowSum += abs(elem);
        if (curRowSum > maxRowSum)
            maxRowSum = curRowSum;
    }
    return maxRowSum;
}

Matrix Zeros(const size_t size) {
    Matrix res = Matrix(size);
    for (auto& row : res)
        for (size_t i = 0; i < size; i++)
            row.push_back(0);
    return res;
}

void SetJacobiIterationsParameters(const Matrix& A, const Vector& b, Matrix& C, Vector& g) {
    const size_t size = A.size();
    for (size_t i = 0; i < size; i++)
        for (size_t j = 0; j < size; j++)
            C[i][j] = i == j ? 0 : -A[i][j] / A[i][i];
    for (size_t i = 0; i < size; i++)
        g[i] = b[i] / A[i][i];
}

double GetMuStopConditionNumber(const Matrix& C) {
    const size_t size = C.size();
    double mu = 0;
    for (size_t i = 0; i < size; i++) {
        double alpha = 0, beta = 0;
        for (size_t j = 0; j < i; j++)
            alpha += abs(C[i][j]);
        for (size_t j = i; j < size; j++)
            beta += abs(C[i][j]);
        double curMu = beta / (1 - alpha);
        if (curMu > mu)
            mu = curMu;
    }
    return mu;
}

Vector SeidelIterations(const Matrix& C, const Vector& g, const Vector& x0, const double eps, int& n) {
    const size_t size = x0.size();
    Vector xn = x0, xn_prev = x0;
    double mu = GetMuStopConditionNumber(C);
    n = 0;
    do {
        xn_prev = xn;
        for (size_t i = 0; i < size; i++) {
            xn[i] = 0;
            for (size_t j = 0; j < i; j++)
                xn[i] += C[i][j] * xn[j];
            for (size_t j = i; j < size; j++)
                xn[i] += C[i][j] * xn_prev[j];
            xn[i] += g[i];
        }
        n++;
    } while (NormInf(VectorSub(xn, xn_prev)) >= (eps * (1 - mu) / mu));
    return xn;
}

Vector SolveSLAUWithSeidelJacobiMethod(const Matrix& A, const Vector& b, const Vector& x0, const double eps, int& n) {
    const size_t size = A.size();
    Matrix C = Zeros(size);
    Vector g(size);
    SetJacobiIterationsParameters(A, b, C, g);
    double normC = NormInf(C);
    if (normC >= 1)
        cout << "Infinite norm of C = " << normC << " > 1, method will diverge!" << endl;
    return SeidelIterations(C, g, x0, eps, n);
}

void FixedStartVectorExperiment(const string& filename, const Matrix& A, const Vector& actualSolution, const Vector& x0,
    const double eps0, const size_t steps) {
    ofstream fileRes(filename);
    if (!fileRes.is_open()) {
        cout << "Could not open " << filename << " file!" << endl;
        return;
    }
    SET_STREAM_PRECISION(fileRes);
    fileRes << "eps; error-inf-norm; error-2nd-norm; residual-inf-norm; residual-2nd-norm; iters" << endl;
    double eps = eps0;
    Vector b = MatrixByVector(A, actualSolution);
    for (size_t i = 0; i < steps; i++) {
        int n = 0;
        Vector x = SolveSLAUWithSeidelJacobiMethod(A, b, x0, eps, n);
        Vector error = VectorSub(x, actualSolution);
        Vector residual = VectorSub(MatrixByVector(A, x), b);
        fileRes << eps << "; " << NormInf(error) << "; " << Norm(error) << "; " << NormInf(residual) << "; " << Norm(residual) << "; " << n << endl;
        eps /= 10;
    }
    fileRes.close();
}