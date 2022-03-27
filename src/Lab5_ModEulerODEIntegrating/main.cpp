#include <functional>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <string>
#include <cassert>
#include <cmath>

#define ROOT "../../"
#define PRECISION 50
#define SET_STREAM_PRECISION(stream) \
    (stream).setf(ios::fixed); \
    (stream) << setprecision(PRECISION)

using namespace std;

struct Point {
    double x, y;
    Point(const double x_, const double y_) : x(x_), y(y_) {}
};

using Grid = vector<double>;
using GridFunc = vector<Point>;
using MathFunc = function<double(double)>;
using ODEFunc = function<double(double, double)>;

double Pow(const double a, const size_t b);
GridFunc SimpleModEulerODE(const ODEFunc& f, const double y0, const double a, const double b, const size_t n);
double SimpleModEulerODELastY(const ODEFunc& f, const double y0, const double a, const double b, const size_t n);
GridFunc ModEulerODERungesAcc(const ODEFunc& f, const double y0, const double a, const double b, const size_t n, const double eps, size_t* maxIters = nullptr);
double MaxError(const MathFunc& f, const GridFunc& gridF);

void WriteGridFunc(const string& filename, const GridFunc& gridFunc);
void WriteODESolutions(const string& filenameGeneric, const ODEFunc& f, const double y0, const double a, const double b, const vector<size_t>& ns);
void WriteErrorOnEps(const string& filename, const MathFunc& sol, const ODEFunc& f,
    const double y0, const double a, const double b, const size_t n, const double eps0, const size_t steps);
void WriteItersOnEps(const string& filename, const ODEFunc& f,
    const double y0, const double a, const double b, const size_t n, const double eps0, const size_t steps);
void WriteErrorOnDy(const string& filename, const MathFunc& sol, const ODEFunc& f, const double y0,
    const double a, const double b, const size_t n, const double eps, const double dyStep, const size_t dySteps);

void WriteErrorOnH(const string& filename, const MathFunc& sol, const ODEFunc& f,
    const double y0, const double a, const double b, const size_t n_, const size_t nSteps);

int main() {
    // не удовлетворяет условию Липшица по y в любом прямоугольнике [1,2]x[0,d] !
    ODEFunc ode = [](double x, double y) { return 4 * y / x + 2 * x * sqrt(y); };
    MathFunc y = [](double x) { return Pow(x, 4) * Pow(log(x) + 1, 2); };
    const double y0 = 1;
    const double a = 1;
    const double b = 2;
    const vector<size_t> ns = { 5, 10 };
    const size_t n = 10;
    const double eps0 = 0.1;
    const size_t epsSteps = 11;
    const double eps = 1e-5;
    const double dyStep = 1e-2;
    const size_t dySteps = 100;
    const size_t nInit = 1;
    const size_t nSteps = 15;

    GridFunc ans = SimpleModEulerODE(ode, y0, a, b, 10);
    for (auto& point : ans)
        cout << point.x << "\t" << point.y << "\t" << y(point.x) << endl;

    try {
        WriteODESolutions(ROOT"csv/ans", ode, y0, a, b, ns);
        WriteErrorOnEps(ROOT"csv/err_on_eps.csv", y, ode, y0, a, b, n, eps0, epsSteps);
        WriteItersOnEps(ROOT"csv/iters_on_eps.csv", ode, y0, a, b, n, eps0, epsSteps);
        WriteErrorOnDy(ROOT"csv/err_on_dy.csv", y, ode, y0, a, b, n, eps, dyStep, dySteps);
        WriteErrorOnH(ROOT"csv/err_on_h.csv", y, ode, y0, a, b, nInit, nSteps);
    }
    catch (const string& err) {
        cout << "Error occured:" << endl << err << endl;
    }

    return 0;
}

double Pow(const double a, const size_t b) {
    double res = 1;
    for (size_t i = 0; i < b; i++)
        res *= a;
    return res;
}

GridFunc SimpleModEulerODE(const ODEFunc& f, const double y0, const double a, const double b, const size_t n) {
    const double h = (b - a) / (double)n;
    double x = a;
    double y = y0;
    GridFunc res = { Point(a, y0) };
    for (size_t i = 0; i < n; i++) {
        y += h * f(x + h / 2.0, y + h / 2.0 * f(x, y));
        x += h;
        res.push_back(Point(x, y));
    }
    return res;
}

double SimpleModEulerODELastY(const ODEFunc& f, const double y0, const double a, const double b, const size_t n) {
    const double h = (b - a) / (double)n;
    double x = a;
    double y = y0;
    for (size_t i = 0; i < n; i++) {
        y += h * f(x + h / 2.0, y + h / 2.0 * f(x, y));
        x += h;
    }
    return y;
}

GridFunc ModEulerODERungesAcc(const ODEFunc& f, const double y0, const double a, const double b, const size_t n, const double eps, size_t* maxIters) {
    const double coef = 3.0;
    const double h = (b - a) / (double)n;
    double x = a;
    double y = y0;
    size_t mIters = 0;
    GridFunc res;
    for (size_t i = 0; i <= n; i++) {
        res.push_back(Point(x, y));
        size_t n = 1, iters = 0;
        double yPrev = SimpleModEulerODELastY(f, y, x, x + h, n);
        double yCur = yPrev;
        do {
            yPrev = yCur;
            n *= 2; ++iters;
            yCur = SimpleModEulerODELastY(f, y, x, x + h, n);
        } while (abs(yPrev - yCur) > coef * eps);
        y = yCur;
        x += h;
        if (iters > mIters)
            mIters = iters;
    }
    if (maxIters)
        *maxIters = mIters;
    return res;
}

double MaxError(const MathFunc& f, const GridFunc& gridF) {
    double res = 0;
    for (auto& p : gridF) {
        double err = abs(f(p.x) - p.y);
        if (err > res)
            res = err;
    }
    return res;
}

void WriteGridFunc(const string& filename, const GridFunc& gridFunc) {
    ofstream file(filename);
    if (!file.is_open())
        throw string("File ") + filename + string(" could not be opened.");
    SET_STREAM_PRECISION(file);
    file << "x;y" << endl;
    for (auto& point : gridFunc)
        file << point.x << ";" << point.y << endl;
    file.close();
}

void WriteODESolutions(const string& filenameGeneric,
                       const ODEFunc& f,
                       const double y0,
                       const double a,
                       const double b,
                       const vector<size_t>& ns) {
    for (size_t i = 0; i < ns.size(); i++)
        WriteGridFunc(filenameGeneric + to_string(i + 1) + ".csv", SimpleModEulerODE(f, y0, a, b, ns[i] - 1));
}

void WriteErrorOnEps(const string& filename,
                     const MathFunc& sol, 
                     const ODEFunc& f, 
                     const double y0, 
                     const double a, 
                     const double b, 
                     const size_t n, 
                     const double eps0, 
                     const size_t steps) {
    ofstream file(filename);
    if (!file.is_open())
        throw string("File ") + filename + string(" could not be opened.");
    SET_STREAM_PRECISION(file);
    file << "eps;err;err(1st)" << endl;
    double eps = eps0;
    for (size_t i = 0; i < steps; i++) {
        GridFunc ans = ModEulerODERungesAcc(f, y0, a, b, n, eps);
        file << eps << ";" << MaxError(sol, ans) << ";"
            << abs(sol(ans[1].x)-ans[1].y) << endl;
        eps /= 10;
    }
    file.close();
}

void WriteItersOnEps(const string& filename,
                     const ODEFunc& f,
                     const double y0,
                     const double a, 
                     const double b, 
                     const size_t n, 
                     const double eps0, 
                     const size_t steps) {
    ofstream file(filename);
    if (!file.is_open())
        throw string("File ") + filename + string(" could not be opened.");
    SET_STREAM_PRECISION(file);
    file << "eps;maxIters;" << endl;
    double eps = eps0;
    for (size_t i = 0; i < steps; i++) {
        size_t iters = 0;
        GridFunc ans = ModEulerODERungesAcc(f, y0, a, b, n, eps, &iters);
        file << eps << ";" << iters << endl;
        eps /= 10;
    }
    file.close();
}

void WriteErrorOnDy(const string& filename, 
                    const MathFunc& sol, 
                    const ODEFunc& f,
                    const double y0, 
                    const double a, 
                    const double b, 
                    const size_t n, 
                    const double eps, 
                    const double dyStep, 
                    const size_t dySteps) {
    ofstream file(filename);
    if (!file.is_open())
        throw string("File ") + filename + string(" could not be opened.");
    SET_STREAM_PRECISION(file);
    file << "dy;err;" << endl;
    double dy = 0;
    for (size_t i = 0; i < dySteps; i++) {
        file << dy << ";" << MaxError(sol, ModEulerODERungesAcc(f, y0 + dy, a, b, n, eps)) << endl;
        dy += dyStep;
    }
    file.close();
}

void WriteErrorOnH(const string& filename,
                   const MathFunc& sol,
                   const ODEFunc& f,
                   const double y0,
                   const double a,
                   const double b,
                   const size_t n_,
                   const size_t nSteps) {
    ofstream file(filename);
    if (!file.is_open())
        throw string("File ") + filename + string(" could not be opened.");
    SET_STREAM_PRECISION(file);
    file << "h;err;" << endl;
    size_t n = n_;
    for (size_t i = 0; i < nSteps; i++) {
        file << (b - a) / (double)n << ";" << MaxError(sol, SimpleModEulerODE(f, y0, a, b, n)) << endl;
        n *= 2;
    }
    file.close();
}
