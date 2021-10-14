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
matrix_t matrix_by_matrix(const matrix_t& A, const matrix_t& B);
matrix_t generate_orthogonal_matrix(const vector_t& w);
double norm(const vector_t& x);
void LU_factorization(const matrix_t& A, matrix_t& L, matrix_t& U);
vector_t forward_substituion(const matrix_t& A, const vector_t& b);
vector_t backward_substituion(const matrix_t& A, const vector_t& b);
vector_t solve_SLAU(const matrix_t& A, const vector_t& b);
void matrix_print(const matrix_t& A);
void vector_print(const vector_t& x);
void matrix_print_csv(const string filename, const matrix_t& A);
void vector_print_csv(const string filename, const vector_t& x);

// For matrix condition number experiment
matrix_t eye(const size_t size);
matrix_t eye_cond_num(const size_t size, const vector_t& w1, const vector_t& w2, const double num);
void eye_cond_experiment(const string filename_eye_cond_solutions,
    const string filename_eye_cond_vectors,
    const string filename_eye_cond_matrices,
    const string filename_eye_cond_actual_solution,
    const vector_t& actual_solution, const vector_t& w1, const vector_t& w2,
    const double num0, const size_t steps);

void good_and_bad_matrices_experiment(const string filename_good_and_bad_results,
    const string filename_good_matrix,
    const string filename_bad_matrix,
    const vector_t& max_delta_vector,
    const vector_t& vector, const size_t steps, const size_t size,
    const matrix_t& good_matrix, const matrix_t& bad_matrix);

int main(int* argc, const char** argv) {
    matrix_t A = { {2, 7, 5, 36, 17, 38, 92, 43, 87, 53},
                   {4, 1, -14, 31, 53, 23, 45, 32, 16, 75},
                   {8, 3, 1, 15, 67, 51, 79, 31, 34, 64},
                   {18, 1, 39, 76, 83, 31, 13, 43, 57, 83},
                   {14, 43, 12, 1, 53, 57, 86, 43, 24, 10},
                   {11, 21, 35, 41, 50, 62, 71, 83, 92, 14},
                   {91, 83, 74, 61, 57, 49, 35, 22, 14, 90},
                   {46, 73, 12, 14, 53, 67, 753, 15, 61, 87},
                   {14, 53, 68, 54, 67, 97, 34, 15, 52, 10},
                   {65, 13, 53, 16, 74, 80, 45, 89, 14, 15} };
    vector_t actual_x = { 7, 9, 3, 2, 1, 5, 8, 4, 6, 10 };
    vector_t b = matrix_by_vector(A, actual_x);

    vector_t x = solve_SLAU(A, b);
    cout << "Solution is:" << endl;
    vector_print(x);

    matrix_print_csv("../csv/matrix.csv", A);
    vector_print_csv("../csv/vector.csv", b);
    vector_print_csv("../csv/actual_solution.csv", actual_x);
    vector_print_csv("../csv/solution.csv", x);

    cout << "Fact error is: " << norm(vector_div(x, actual_x));

    //experiment on condition number of matrix inequality
    vector_t delta_b = { 0.1, 0.2, 0.05, 0.06, 0.03, 0.21, 0.13, 0.02, 0.1, 0.13 };
    vector_t delta_solution = vector_div(x, solve_SLAU(A, vector_add(b, delta_b)));

    vector_print_csv("../csv/delta_vector.csv", delta_b);
    vector_print_csv("../csv/delta_solution.csv", delta_solution);
    
    //experiment on eye condition number
    vector_t w1 = { 56, 9, 3, 2, 1, 5, 8, 4, 6, 10 };
    vector_t w2 = { 3, 2, 3, 9, 32, 5, 8, 75, 4, 1 };
    eye_cond_experiment("../csv/eye_cond_solutions.csv",
        "../csv/eye_cond_vectors.csv",
        "../csv/eye_cond_matrices.csv",
        "../csv/eye_cond_actual_solution.csv",
        actual_x, w1, w2, 10, 10);
  
    //experiment on 10 points
    vector_t good_and_bad_exp_vec = { 75, 91, 32, 25, 16, 51, 85, 48, 60, 102 };
    vector_t delta_vec = { 0.135, 0.3, 0.531, 0.535, 0.6536, 0.131, 0.35, 0.3563, 0.35, 0.3342 };
    vector_t w11 = w1;
    vector_t w12 = w2;
    vector_t w21 = { 27, 23, 163, 363, 35, 351, 312, 13, 35, 31 };
    vector_t w22 = { 14, 98, 53, 378, 56, 57, 32, 103, 32, 57 };
    good_and_bad_matrices_experiment("../csv/good_and_bad_res.csv", 
        "../csv/good_matrix.csv", "../csv/bad_matrix.csv", delta_vec, good_and_bad_exp_vec, 30, 10, eye_cond_num(10, w11, w12, 10), eye_cond_num(10, w21, w22, 1000));

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

matrix_t matrix_by_matrix(const matrix_t& A, const matrix_t& B) {
    const size_t size = A.size();
    matrix_t res = zeros(size);
    for (size_t i = 0; i < size; i++)
        for (size_t j = 0; j < size; j++) {
            double temp = 0;
            for (size_t k = 0; k < size; k++)
                temp += A[i][k] * B[k][j];
            res[i][j] = temp;
        }
    return res;
}

matrix_t generate_orthogonal_matrix(const vector_t& w) {
    const size_t size = w.size();
    matrix_t res = eye(size);
    double norm_w = norm(w);
    for (size_t i = 0; i < size; i++)
        for (size_t j = 0; j < size; j++)
            res[i][j] -= 2 * w[i] * w[j] / (norm_w * norm_w);
    return res;
}

double norm(const vector_t& x) {
    double norm_sqr = 0;
    for (auto elem : x)
        norm_sqr += elem * elem;
    return sqrt(norm_sqr);
}

void LU_factorization(const matrix_t& A, matrix_t& L, matrix_t& U) {
    const int N = A.size();
    for (int m = 0; m < N; m++) {
        for (int j = m; j < N; j++) {
            U[m][j] = A[m][j];
            for (int k = 0; k < m; k++)
                U[m][j] -= L[m][k] * U[k][j];
        }
        for (int i = m + 1; i < N; i++) {
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
    matrix_t L = eye(N);
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

matrix_t eye(const size_t size) {
    matrix_t eye = zeros(size);
    for (size_t i = 0; i < size; i++)
        eye[i][i] = 1;

    return eye;
}

matrix_t eye_cond_num(const size_t size, const vector_t& w1, const vector_t& w2, const double num) {
    matrix_t res = eye(size);
    for (size_t i = 0; i < size; i++)
        res[i][i] = num - i * (double)(num - 1) / (size - 1);
    // generating orthogonal matrices
    matrix_t Q1 = generate_orthogonal_matrix(w1);
    matrix_t Q2 = generate_orthogonal_matrix(w2);
    // multiplying Q1*res
    matrix_t Q1_by_res = matrix_by_matrix(Q1, res);
    // multiplying Q1_by_res * Q2
    matrix_t Q1_by_res_by_Q2 = matrix_by_matrix(Q1_by_res, Q2);
    return Q1_by_res_by_Q2;
}

void eye_cond_experiment(const string filename_eye_cond_solutions,
    const string filename_eye_cond_vectors,
    const string filename_eye_cond_matrices,
    const string filename_eye_cond_actual_solution,
    const vector_t& actual_solution, const vector_t& w1, const vector_t& w2,
    const double num0, const size_t steps) {
    const size_t size = actual_solution.size();
    vector_print_csv(filename_eye_cond_actual_solution, actual_solution);
    ofstream file_eye_cond_solutions(filename_eye_cond_solutions);
    if (!file_eye_cond_solutions.is_open()) {
        cerr << "File " << filename_eye_cond_solutions << " could not be opened!" << endl;
        return;
    }
    ofstream file_eye_cond_matrices(filename_eye_cond_matrices);
    if (!file_eye_cond_matrices.is_open()) {
        cerr << "File " << filename_eye_cond_matrices << " could not be opened!" << endl;
        return;
    }
    ofstream file_eye_cond_vectors(filename_eye_cond_vectors);
    if (!file_eye_cond_vectors.is_open()) {
        cerr << "File " << filename_eye_cond_vectors << " could not be opened!" << endl;
        return;
    }
    file_eye_cond_solutions.setf(ios::fixed);
    file_eye_cond_matrices.setf(ios::fixed);
    file_eye_cond_vectors.setf(ios::fixed);
    double num = num0;
    for (size_t i = 0; i < steps; i++) {
        matrix_t matrix = eye_cond_num(size, w1, w2, num);
        vector_t vector = matrix_by_vector(matrix, actual_solution);
        vector_t x = solve_SLAU(matrix, vector);
        for (size_t j = 0; j < size - 1; j++)
            file_eye_cond_solutions << setprecision(PRECISION) << x[j] << ";";
        file_eye_cond_solutions << setprecision(PRECISION) << x[size - 1] << endl;
        for (size_t j = 0; j < size - 1; j++)
            file_eye_cond_vectors << setprecision(PRECISION) << vector[j] << ";";
        file_eye_cond_vectors << setprecision(PRECISION) << vector[size - 1] << endl;
        for (size_t j = 0; j < size; j++) {
            for (size_t k = 0; k < size - 1; k++)
                file_eye_cond_matrices << setprecision(PRECISION) << matrix[j][k] << ";";
            file_eye_cond_matrices << setprecision(PRECISION) << matrix[j][size - 1] << endl;
        }
        num *= 10;
    }
}

void good_and_bad_matrices_experiment(const string filename_good_and_bad_results,
    const string filename_good_matrix,
    const string filename_bad_matrix,
    const vector_t& max_delta_vector,
    const vector_t& vector, const size_t steps, const size_t size,
    const matrix_t& good_matrix, const matrix_t& bad_matrix) {
    matrix_print_csv(filename_good_matrix, good_matrix);
    matrix_print_csv(filename_bad_matrix, bad_matrix);
    ofstream file_res(filename_good_and_bad_results);
    if (!file_res.is_open()) {
        cerr << "File " << filename_good_and_bad_results << " could not be opened!" << endl;
        return;
    }
    file_res.setf(ios::fixed);
    for (size_t i = 0; i < steps; i++) {
        vector_t delta_vector(size);
        for (size_t k = 0; k < size; k++)
            delta_vector[k] = max_delta_vector[k] * (i + 1.0) / steps;
        file_res << setprecision(PRECISION) << norm(delta_vector) / norm(vector) << ";";
        // solving for good matrix
        vector_t solution = solve_SLAU(good_matrix, vector);
        vector_t x = solve_SLAU(good_matrix, vector_add(vector, delta_vector));
        vector_t delta_x = vector_div(solution, x);
        file_res << setprecision(PRECISION) << norm(delta_x) / norm(solution) << ";";
        // solving for bad matrix
        solution = solve_SLAU(bad_matrix, vector);
        x = solve_SLAU(bad_matrix, vector_add(vector, delta_vector));
        delta_x = vector_div(solution, x);
        file_res << setprecision(PRECISION) << norm(delta_x) / norm(solution);
        file_res << endl;
    }
}