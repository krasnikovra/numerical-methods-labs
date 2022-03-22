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
#define T -0.9

using namespace std;
using MathFunc = double (*)(double);
using Grid = vector<double>;

const Grid LEGENDRE_GRID = Grid{ -sqrt(0.6), 0, sqrt(0.6) };
const Grid LEGENDRE_QUOT = Grid{ 5.0, 8.0, 5.0 };

double Pow(double a, size_t b);
Grid MakeUniformGrid(const double a, const double b, const size_t pointsCount);
Grid TranslateLegendreGrid(const double a, const double b);
double SimpleGaussIntegral(MathFunc f, const double a, const double b);
double GaussIntegral(MathFunc f, const double a, const double b, const size_t intervalsCount);
double EvaluateIntegralWithRungesAccuracy(MathFunc f, const double a, const double b, const double eps, size_t* q = nullptr, double* richardson = nullptr);

void WriteErrorOnEps(const string& filename, MathFunc f, const double accInt, const double a, const double b, const double eps0, const size_t steps);
void WriteQOnEps(const string& filename, MathFunc f, const double a, const double b, const double eps0, const size_t steps);
void WriteErrorOnH(const string& filename, MathFunc f, const double accInt, const double a, const double b, const size_t steps);

double f(double x) {
    return (Pow(x, 5) - 5.2 * Pow(x, 3) + 5.5 * Pow(x, 2) - 7 * x - 3.5) * cos(0.4 * x);
}

double F(double x) {
    return cos(0.4 * x) * (30471.875 + 68.75 * x - 2441.25 * Pow(x, 2) + 31.25 * Pow(x, 4)) +
        sin(0.4 * x) * (-180.625 + 12188.75 * x + 13.75 * Pow(x, 2) - 325.5 * Pow(x, 3) + 2.5 * Pow(x, 5));
}

double fN(double x) {
    return x < T ? f(x) : 2 * f(T) - f(x);
}

int main() {
    const double a = -3.0;
    const double b = 0.7;
    const double eps0 = 0.1;
    const size_t epsSteps = 11;
    const size_t hSteps = 11;

    const double accInt = F(b) - F(a);
    const double accIntN = F(T) - F(a) + 2 * b * f(T) - F(b) - 2 * T * f(T) + F(T);

    try {
        WriteErrorOnEps(ROOT"csv/error_on_eps.csv", f, accInt, a, b, eps0, epsSteps);
        WriteQOnEps(ROOT"csv/q_on_eps.csv", f, a, b, eps0, epsSteps);
        WriteErrorOnH(ROOT"csv/error_on_h.csv", f, accInt, a, b, hSteps);

        WriteErrorOnEps(ROOT"csv/error_on_epsn.csv", fN, accIntN, a, b, eps0, epsSteps);
        WriteQOnEps(ROOT"csv/q_on_epsn.csv", fN, a, b, eps0, epsSteps);
        WriteErrorOnH(ROOT"csv/error_on_hn.csv", fN, accIntN, a, b, hSteps);

        cout << EvaluateIntegralWithRungesAccuracy(fN, a, b, 1e-5) << endl;
        cout << accIntN << endl;
    }
    catch (const string& err) {
        cout << "Error occured:" << endl << err << endl;
    }
    return 0;
}

double Pow(double a, size_t b) {
    double res = 1;
    for (size_t i = 0; i < b; i++)
        res *= a;
    return res;
}

Grid MakeUniformGrid(const double a, const double b, const size_t pointsCount) {
    Grid res;
    for (size_t i = 0; i < pointsCount; i++)
        res.push_back(a + (b - a) * i / (double)(pointsCount - 1));
    return res;
}

Grid TranslateLegendreGrid(const double a, const double b) {
    Grid res;
    for (auto& t : LEGENDRE_GRID)
        res.push_back((b + a) / 2 + (b - a) / 2 * t);
    return res;
}

double SimpleGaussIntegral(MathFunc f, const double a, const double b) {
    double res = 0;
    Grid grid = TranslateLegendreGrid(a, b);
    for (size_t i = 0; i < grid.size(); i++)
        res += LEGENDRE_QUOT[i] * f(grid[i]);
    return res * (b - a) / 18;
}

double GaussIntegral(MathFunc f, const double a, const double b, const size_t intervalsCount) {
    double res = 0;
    Grid intervalPartition = MakeUniformGrid(a, b, intervalsCount + 1);
    for (size_t i = 0; i + 1 < intervalPartition.size(); i++)
        res += SimpleGaussIntegral(f, intervalPartition[i], intervalPartition[i + 1]);
    return res;
}

double EvaluateIntegralWithRungesAccuracy(MathFunc f, const double a, const double b, const double eps, size_t* q, double* richardson) {
    size_t intervalsCount = 1;
    size_t iters = 0;
    const double coef = 63; //2^6 - 1
    double iPrev = GaussIntegral(f, a, b, intervalsCount);
    double i = iPrev;
    do {
        intervalsCount *= 2;
        ++iters;
        iPrev = i;
        i = GaussIntegral(f, a, b, intervalsCount);
    } while (abs(i - iPrev) >= eps * coef);
    if (q != nullptr)
        *q = iters;
    if (richardson != nullptr)
        *richardson = ((coef + 1) * i - iPrev) / coef;
    return i;
}

void WriteErrorOnEps(const string& filename, MathFunc f, const double accInt, const double a, const double b, const double eps0, const size_t steps) {
    ofstream file(filename);
    if (!file.is_open())
        throw string("File ") + filename + string(" could not be opened.");
    SET_STREAM_PRECISION(file);
    file << "eps;err;err_richardson" << endl;
    double eps = eps0;
    for (size_t i = 0; i < steps; i++) {
        double rich = 0.0;
        file << eps << ";" << abs(accInt - EvaluateIntegralWithRungesAccuracy(f, a, b, eps, nullptr, &rich)) << ";";
        file << abs(accInt - rich) << endl;
        eps /= 10;
    }
    file.close();
}

void WriteQOnEps(const string& filename, MathFunc f, const double a, const double b, const double eps0, const size_t steps) {
    ofstream file(filename);
    if (!file.is_open())
        throw string("File ") + filename + string(" could not be opened.");
    SET_STREAM_PRECISION(file);
    file << "eps;q" << endl;
    double eps = eps0;
    for (size_t i = 0; i < steps; i++) {
        size_t q = 0;
        EvaluateIntegralWithRungesAccuracy(f, a, b, eps, &q);
        file << eps << ";" << q << endl;
        eps /= 10;
    }
    file.close();
}

void WriteErrorOnH(const string& filename, MathFunc f, const double accInt, const double a, const double b, const size_t steps) {
    ofstream file(filename);
    if (!file.is_open())
        throw string("File ") + filename + string(" could not be opened.");
    SET_STREAM_PRECISION(file);
    file << "h;err" << endl;
    size_t n = 1;
    for (size_t i = 0; i < steps; i++) {
        file << (b - a) / n << ";" << abs(accInt - GaussIntegral(f, a, b, n)) << endl;
        n *= 2;
    }
    file.close();
}