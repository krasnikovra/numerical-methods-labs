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
using Grid = vector<double>;

const Grid LEGENDRE_GRID = Grid{ -sqrt(0.6), 0, sqrt(0.6) };
const Grid LEGENDRE_QUOT = Grid{ 5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0 };

double Pow(double a, size_t b);
Grid MakeUniformGrid(const double a, const double b, const size_t pointsCount);
Grid TranslateLegendreGrid(const double a, const double b);
double SimpleGaussIntegral(MathFunc f, const double a, const double b);
double GaussIntegral(MathFunc f, const double a, const double b, const size_t intervalsCount);
double EvaluateIntegralWithRungesAccuracy(MathFunc f, const double a, const double b, const double eps, size_t* q = nullptr);

void WriteErrorOnEps(const string& filename, MathFunc f, MathFunc F, const double a, const double b, const double eps0, const size_t steps);
void WriteQOnEps(const string& filename, MathFunc f, MathFunc F, const double a, const double b, const double eps0, const size_t steps);

double f(double x) {
    return (Pow(x, 5) - 5.2 * Pow(x, 3) + 5.5 * Pow(x, 2) - 7 * x - 3.5) * cos(0.4 * x);
}

double F(double x) {
    return (243775 * cos((2 * x) / 5.0)) / 8.0 - (1445 * sin((2 * x) / 5.0)) / 8.0 + (275 * x * cos((2 * x) / 5.0)) / 4.0
        + (48755 * x * sin((2 * x) / 5.0)) / 4.0 - (9765 * Pow(x, 2) * cos((2 * x) / 5.0)) / 4.0 + (125 * Pow(x, 4) * cos((2 * x) / 5.0)) / 4.0
        + (55 * Pow(x, 2) * sin((2 * x) / 5.0)) / 4.0 - (651 * Pow(x, 3) * sin((2 * x) / 5.0)) / 2.0 + (5 * Pow(x, 5) * sin((2 * x) / 5.0)) / 2.0;
}

int main() {
    const double a = -3.0;
    const double b = 0.7;
    const double eps0 = 0.1;
    const size_t epsSteps = 11;
    try {
        WriteErrorOnEps(ROOT"csv/error_on_eps.csv", f, F, a, b, eps0, epsSteps);
        WriteQOnEps(ROOT"csv/q_on_eps.csv", f, F, a, b, eps0, epsSteps);
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
    return res * (b - a) / 2;
}

double GaussIntegral(MathFunc f, const double a, const double b, const size_t intervalsCount) {
    double res = 0;
    Grid intervalPartition = MakeUniformGrid(a, b, intervalsCount + 1);
    for (size_t i = 0; i + 1 < intervalPartition.size(); i++)
        res += SimpleGaussIntegral(f, intervalPartition[i], intervalPartition[i + 1]);
    return res;
}

double EvaluateIntegralWithRungesAccuracy(MathFunc f, const double a, const double b, const double eps, size_t* q) {
    size_t intervalsCount = 1;
    size_t iters = 0;
    const size_t coef = 31; //2^5 - 1
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
        size_t iters = 0;
        file << eps << ";" << abs(exactIntegral - EvaluateIntegralWithRungesAccuracy(f, a, b, eps, &iters)) << endl;
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