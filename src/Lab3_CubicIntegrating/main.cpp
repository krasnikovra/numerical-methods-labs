#define _USE_MATH_DEFINES

#include <iostream>
#include <fstream>
#include <vector>
#include <iomanip>
#include <math.h>

#define ROOT "../../"
#define PRECISION 50
#define SET_STREAM_PRECISION(stream) \
    (stream).setf(ios::fixed); \
    (stream) << setprecision(PRECISION)

using namespace std;

using MathFunc = double (*)(double);

double Pow(double a, int b);
double CubicNewtonCotesIntegral(MathFunc f, const double a, const double b, const size_t intervalsCount);
double EvaluateIntegralWithRungesAccuracy(MathFunc f, const double a, const double b, const double eps, size_t* q = nullptr);

void WriteErrorOnEps(const string& filename, MathFunc f, MathFunc F, const double a, const double b, const double eps0, const size_t steps);
void WriteQOnEps(const string& filename, MathFunc f, MathFunc F, const double a, const double b, const double eps0, const size_t steps);
void WriteErrorOnH(const string& filename, MathFunc f, MathFunc F, const double a, const double b, const size_t steps);
void WriteConstApprox(const string& filename, MathFunc f, MathFunc F, const double a, const double b, const size_t steps);

int main() {
    const auto f = [](double x) { return Pow(x, 5) - 5.2 * Pow(x, 3) + 5.5 * Pow(x, 2) - 7 * x - 3.5; };
    const auto F = [](double x) { return Pow(x, 6) / 6 - 5.2 * Pow(x, 4) / 4 + 5.5 * Pow(x, 3) / 3 - 7 * Pow(x, 2) / 2 - 3.5 * x; };
    const double a = -3;
    const double b = 0.7;
    const double eps = 0.0001;
    const double eps0 = 0.1;
    const size_t epsSteps = 12;
    const size_t hSteps = 12;
    cout << EvaluateIntegralWithRungesAccuracy(f, a, b, eps) << endl;
    try {
        WriteErrorOnEps(ROOT"csv/error_on_eps.csv", f, F, a, b, eps0, epsSteps);
        WriteQOnEps(ROOT"csv/q_on_eps.csv", f, F, a, b, eps0, epsSteps);
        WriteErrorOnH(ROOT"csv/error_on_h.csv", f, F, a, b, hSteps);
        WriteConstApprox(ROOT"csv/const_approx.csv", f, F, a, b, hSteps);
    }
    catch (const string& err) {
        cout << "Error occured:" << endl << err << endl;
    }
    return 0;
}

double Pow(double a, int b) {
    double res = 1;
    for (int i = 0; i < b; i++)
        res *= a;
    return res;
}

double CubicNewtonCotesIntegral(MathFunc f, const double a, const double b, const size_t intervalsCount) {
    const size_t n = 3 * intervalsCount;
    double res = 0;
    const double h = (b - a) / n;
    res += f(a) + 3 * f(a + (n - 2) * h) + 3 * f(a + (n - 1) * h) + f(b);
    for (size_t i = 1; i <= n / 3 - 1; i++)
        res += 3 * f(a + (3 * i - 2) * h) + 3 * f(a + (3 * i - 1) * h) + 2 * f(a + (3 * i) * h);
    res *= 3.0 * h / 8.0;
    return res;
}

double EvaluateIntegralWithRungesAccuracy(MathFunc f, const double a, const double b, const double eps, size_t* q) {
    size_t intervalsCount = 1;
    size_t iters = 0;
    const int k = 4; // method's order
    const double coef = Pow(2, k) - 1;
    double iPrev = CubicNewtonCotesIntegral(f, a, b, intervalsCount);
    double i = iPrev;
    do {
        intervalsCount *= 2;
        ++iters;
        iPrev = i;
        i = CubicNewtonCotesIntegral(f, a, b, intervalsCount);
    } while (abs(i - iPrev) >= eps * coef);
    if (q != nullptr)
        *q = iters;
    return i;
}

void WriteErrorOnEps(const string& filename, MathFunc f, MathFunc F, const double a, const double b, const double eps0, const size_t steps) {
    ofstream file(filename);
    if (!file.is_open())
        throw string("File ") + filename + string(" could not be opened.");
    SET_STREAM_PRECISION(file);
    file << "eps;err" << endl;
    double eps = eps0;
    const double exactIntegral = F(b) - F(a);
    for (size_t i = 0; i < steps; i++) {
        file << eps << ";" << abs(exactIntegral - EvaluateIntegralWithRungesAccuracy(f, a, b, eps)) << endl;
        eps /= 10;
    }
    file.close();
}

void WriteQOnEps(const string& filename, MathFunc f, MathFunc F, const double a, const double b, const double eps0, const size_t steps) {
    ofstream file(filename);
    if (!file.is_open())
        throw string("File ") + filename + string(" could not be opened.");
    SET_STREAM_PRECISION(file);
    file << "eps;q" << endl;
    const double exactIntegral = F(b) - F(a);
    double eps = eps0;
    for (size_t i = 0; i < steps; i++) {
        size_t q = 0;
        EvaluateIntegralWithRungesAccuracy(f, a, b, eps, &q);
        file << eps << ";" << q << endl;
        eps /= 10;
    }
    file.close();
}

void WriteErrorOnH(const string& filename, MathFunc f, MathFunc F, const double a, const double b, const size_t steps) {
    ofstream file(filename);
    if (!file.is_open())
        throw string("File ") + filename + string(" could not be opened.");
    SET_STREAM_PRECISION(file);
    file << "h;err" << endl;
    const double exactIntegral = F(b) - F(a);
    size_t n = 1;
    for (size_t i = 0; i < steps; i++) {
        file << (b - a) / n << ";" << abs(exactIntegral - CubicNewtonCotesIntegral(f, a, b, n)) << endl;
        n *= 2;
    }
    file.close();
}

void WriteConstApprox(const string& filename, MathFunc f, MathFunc F, const double a, const double b, const size_t steps) {
    ofstream file(filename);
    if (!file.is_open())
        throw string("File ") + filename + string(" could not be opened.");
    SET_STREAM_PRECISION(file);
    file << "h;const" << endl;
    const double exactIntegral = F(b) - F(a);
    size_t n = 1;
    for (size_t i = 0; i < steps; i++) {
        file << (b - a) / n << ";" << abs(exactIntegral - CubicNewtonCotesIntegral(f, a, b, n)) / Pow((b - a) / n, 4) << endl;
        n *= 2;
    }
    file.close();
}