#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>

#define PRECISION 30
//#define DEBUG_OUT

using namespace std;

typedef vector<vector<double>> matrix_t;
typedef vector<double> vector_t;

matrix_t zeros(const size_t size);
vector_t vector_add(const vector_t& a, const vector_t& b);
vector_t vector_div(const vector_t& a, const vector_t& b);
vector_t matrix_by_vector(const matrix_t& A, const vector_t& x);
void LU_factorization(const matrix_t& A, matrix_t& L, matrix_t& U);
vector_t forward_substituion(const matrix_t& A, const vector_t& b);
vector_t backward_substituion(const matrix_t& A, const vector_t& b);
vector_t solve_SLAU(const matrix_t& A, const vector_t& b);
void matrix_print(const matrix_t& A);
void vector_print(const vector_t& x);
void matrix_print_csv(const string filename, const matrix_t& A);
void vector_print_csv(const string filename, const vector_t& x);

// For hilbelrt matrix experiment in matlab
matrix_t hilbert(const size_t size);
vector_t hilbert_vector(const size_t size);
void hilbert_expetiment(const string filename_hilbert_solutions, const string filename_hilbert_vectors, const size_t n0, const size_t steps);

// For matrix condition number experiment
matrix_t eye(const size_t size);
matrix_t eye_plus_num(const size_t size, const double num);
void eye_cond_num_experiment(const string filename_eye_cond_num_solutions,
    const string filename_eye_cond_vector, const string filename_eye_cond_num_matrices,
    const vector_t& vector, const size_t size,
    const double num0, const size_t steps);

int main(int* argc, const char** argv) {
    matrix_t A = { {2, 7, 5, 36, 17, 38, 92, 43, 87, 53},
                   {4, 1, -14, 31, 53, 23, 45, 32, 16, 75},
                   {8, 3, 1, 15, 67, 51, 79, 31, 34, 64},
                   {18, 1, 39, 76, 83, 31, 13, 43, 57, 83},
                   {14, 43, 12, 1, 53, 57, 86, 43, 24, 10},
                   {11, 21, 35, 41, 50, 62, 71, 83, 92, 14},
                   {91, 83, 74, 61, 57, 49, 35, 22, 14, 90},
                   {46, 73, 12, 14, 53, 67, 753, 15, 61, 87} ,
                   {14, 53, 68, 54, 67, 97, 34, 15, 52, 10},
                   {65, 13, 53, 16, 74, 80, 45, 89, 14, 15} };
    vector_t b = { 1, 2, 4, 7, 3, 14, 61, 21, 18, 31 };

    vector_t x = solve_SLAU(A, b);
    cout << "Solution is:" << endl;
    vector_print(x);

    cout << "A*x is:" << endl;
    vector_print(matrix_by_vector(A, x));

    //export data for matlab analysis
    matrix_print_csv("../csv/matrix.csv", A);
    vector_print_csv("../csv/vector.csv", b);
    vector_print_csv("../csv/solution.csv", x);

    //experiment on condition number of matrix
    vector_t delta_b = { 0.1, 0.2, 0.05, 0.06, 0.03, 0.21, 0.13, 0.02, 0.1, 0.13 };
    vector_t delta_solution = vector_div(x, solve_SLAU(A, vector_add(b, delta_b)));

    vector_print_csv("../csv/delta_vector.csv", delta_b);
    vector_print_csv("../csv/delta_solution.csv", delta_solution);

    //experiment on hilbert matrix
    hilbert_expetiment("../csv/hilbert_solutions.csv", "../csv/hilbert_vectors.csv", 1, 25);

    //expetiment on eye condition number
    eye_cond_num_experiment("../csv/eye_cond_num_solutions.csv",
        "../csv/eye_cond_num_vector.csv", "../csv/eye_cond_num_matrices.csv", b, 10, 10, 10);

    return 0;
}

matrix_t zeros(const size_t size) {
    matrix_t matrix(size);
    for (auto iter = matrix.begin(); iter < matrix.end(); iter++) {
        for (size_t i = 0; i < matrix.size(); i++)
            iter->push_back(0);
    }
    return matrix;
}

vector_t vector_add(const vector_t& a, const vector_t& b) {
    const int N = a.size();
    vector_t res(N);
    for (int i = 0; i < N; i++)
        res[i] = a[i] + b[i];
    return res;
}

vector_t vector_div(const vector_t& a, const vector_t& b) {
    const int N = a.size();
    vector_t res(N);
    for (int i = 0; i < N; i++)
        res[i] = a[i] - b[i];
    return res;
}

vector_t matrix_by_vector(const matrix_t& A, const vector_t& x) {
    size_t size = x.size();
    vector_t vector(size);
    for (size_t i = 0; i < size; i++)
        for (size_t j = 0; j < size; j++)
            vector[i] += A[i][j] * x[j];
    return vector;
}

void LU_factorization(const matrix_t& A, matrix_t& L, matrix_t& U) {
    const int N = A.size();
    for (int m = 0; m < N; m++) {
        for (int j = m; j < N; j++) {
            U[m][j] = A[m][j];
            for (int k = 0; k < m; k++)
                U[m][j] -= L[m][k] * U[k][j];
        }
        for (int i = m; i < N; i++) {
            L[i][m] = A[i][m];
            for (int k = 0; k < m; k++)
                L[i][m] -= L[i][k] * U[k][m];
            L[i][m] /= U[m][m];
        }
    }
}

vector_t forward_substituion(const matrix_t& A, const vector_t& b) {
    const int N = A.size();
    vector_t x(N);
    for (int i = 0; i < N; i++) {
        x[i] = b[i];
        for (int j = 0; j < i; j++)
            x[i] -= A[i][j] * x[j];
        x[i] /= A[i][i];
    }
    return x;
}

vector_t backward_substituion(const matrix_t& A, const vector_t& b) {
    const int N = A.size();
    vector_t x(N);
    for (int i = N - 1; i >= 0; i--) {
        x[i] = b[i];
        for (int j = N - 1; j > i; j--)
            x[i] -= A[i][j] * x[j];
        x[i] /= A[i][i];
    }
    return x;
}

vector_t solve_SLAU(const matrix_t& A, const vector_t& b) {
#ifdef DEBUG_OUT
    cout << "Solving a SLAU with matrix A:" << endl;
    matrix_print(A);
    cout << "and vector b:" << endl;
    vector_print(b);
#endif

    const int N = A.size();
    matrix_t L = zeros(N);
    matrix_t U = zeros(N);
    LU_factorization(A, L, U);

#ifdef DEBUG_OUT
    cout << "matrix L:" << endl;
    matrix_print(L);
    cout << "matrix U:" << endl;
    matrix_print(U);
#endif

    vector_t y = forward_substituion(L, b);
#ifdef DEBUG_OUT
    cout << "vector y: Ly = b:" << endl;
    vector_print(y);
#endif

    vector_t x = backward_substituion(U, y);
#ifdef DEBUG_OUT
    cout << "vector x: Ux = y:" << endl;
    vector_print(x);
#endif

    return x;
}

void matrix_print(const matrix_t& A) {
    for (auto iter_row : A) {
        for (auto iter : iter_row)
            cout << '\t' << iter;
        cout << endl;
    }
}

void vector_print(const vector_t& x) {
    for (auto iter = x.begin(); iter < x.end(); iter++)
        cout << '\t' << *iter << endl;
}

void matrix_print_csv(const string filename, const matrix_t& A) {
    ofstream file(filename);
    if (!file.is_open()) {
        cerr << "File " << filename << " could not be opened!" << endl;
        return;
    }
    file.setf(ios::fixed);
    for (auto iter_row : A) {
        for (auto iter = iter_row.begin(); iter < iter_row.end(); iter++) {
            file << setprecision(PRECISION) << *iter;
            if (iter + 1 < iter_row.end())
                file << ";";
        }
        file << endl;
    }
    file.close();
}

void vector_print_csv(const string filename, const vector_t& x) {
    ofstream file(filename);
    if (!file.is_open()) {
        cerr << "File " << filename << " could not be opened!" << endl;
        return;
    }
    file.setf(ios::fixed);
    for (auto iter = x.begin(); iter < x.end(); iter++) {
        file << setprecision(PRECISION) << *iter << endl;
    }
    file.close();
}

matrix_t hilbert(const size_t size) {
    matrix_t matrix = zeros(size);
    for (size_t i = 0; i < size; i++)
        for (size_t j = 0; j < size; j++)
            matrix[i][j] = 1 / (i + j + 1.0);
    return matrix;
}

vector_t hilbert_vector(const size_t size) {
    vector_t vector(size);
    for (size_t i = 0; i < size; i++)
        // just for example, could be any other formula
        vector[i] = 1.0 - (1.0 / (i + 2.0));
    return vector;
}

void hilbert_expetiment(const string filename_hilbert_solutions, const string filename_hilbert_vectors, const size_t n0, const size_t steps) {
    ofstream file_sol(filename_hilbert_solutions);
    if (!file_sol.is_open()) {
        cerr << "File " << filename_hilbert_solutions << " could not be opened!" << endl;
        return;
    }
    ofstream file_vec(filename_hilbert_vectors);
    if (!file_sol.is_open()) {
        cerr << "File " << filename_hilbert_vectors << " could not be opened!" << endl;
        file_sol.close();
        return;
    }
    file_sol.setf(ios::fixed);
    file_vec.setf(ios::fixed);
    for (size_t i = 0; i < steps; i++) {
        size_t n = n0 + i;
        matrix_t A = hilbert(n);
        vector_t b = hilbert_vector(n);
        vector_t solution = solve_SLAU(A, b);
        file_sol << n << ";";
        for (auto iter = solution.begin(); iter < solution.end(); iter++) {
            file_sol << setprecision(PRECISION) << *iter;
            if (iter + 1 < solution.end())
                file_sol << ";";
        }
        for (size_t j = 0; j < n0 + steps - n; j++) {
            file_sol << 0;
            if (j < n0 + steps - 1 - n)
                file_sol << ";";
        }
        file_sol << endl;
        file_vec << n << ";";
        for (auto iter = b.begin(); iter < b.end(); iter++) {
            file_vec << setprecision(PRECISION) << *iter;
            if (iter + 1 < b.end())
                file_vec << ";";
        }
        for (size_t j = 0; j < n0 + steps - n; j++) {
            file_vec << 0;
            if (j < n0 + steps - 1 - n)
                file_vec << ";";
        }
        file_vec << endl;
    }
    file_sol.close();
    file_vec.close();
}

matrix_t eye(const size_t size) {
    matrix_t eye = zeros(size);
    for (size_t i = 0; i < size; i++)
        eye[i][i] = 1;
    return eye;
}

matrix_t eye_plus_num(const size_t size, const double num) {
    matrix_t res = eye(size);
    for (size_t i = 0; i < size; i++)
        for (size_t j = 0; j < size; j++)
            res[i][j] += num;
    for (size_t i = 0; i < size; i++)
        for (size_t j = 0; j < i; j++)
            res[i][j] -= 1.0;
    return res;
}

void eye_cond_num_experiment(const string filename_eye_cond_num_solutions,
    const string filename_eye_cond_vector, const string filename_eye_cond_num_matrices,
    const vector_t& vector, const size_t size, 
    const double num0, const size_t steps) {
    vector_print_csv(filename_eye_cond_vector, vector);
    ofstream file_eye_cond_num_solutions(filename_eye_cond_num_solutions);
    if (!file_eye_cond_num_solutions.is_open()) {
        cerr << "File " << filename_eye_cond_num_solutions << " could not be opened!" << endl;
        return;
    }
    ofstream file_eye_cond_num_matrices(filename_eye_cond_num_matrices);
    if (!file_eye_cond_num_matrices.is_open()) {
        cerr << "File " << filename_eye_cond_num_matrices << " could not be opened!" << endl;
        return;
    }
    file_eye_cond_num_solutions.setf(ios::fixed);
    file_eye_cond_num_matrices.setf(ios::fixed);
    double num = num0;
    for (size_t i = 0; i < steps; i++) {
        matrix_t matrix = eye_plus_num(size, num);
        vector_t x = solve_SLAU(matrix, vector);
        for (size_t j = 0; j < size - 1; j++)
            file_eye_cond_num_solutions << setprecision(PRECISION) << x[j] << ";";
        file_eye_cond_num_solutions << setprecision(PRECISION) << x[size - 1] << endl;
        for (size_t j = 0; j < size; j++) {
            for (size_t k = 0; k < size - 1; k++)
                file_eye_cond_num_matrices << setprecision(PRECISION) << matrix[j][k] << ";";
            file_eye_cond_num_matrices << setprecision(PRECISION) << matrix[j][size - 1] << endl;
        }
        num *= 10;
    }
}