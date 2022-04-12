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
using ODEMethod = function<GridFunc(const ODEFunc&, const double, const double, const double, const size_t)>;

double Pow(const double a, const size_t b);
double EulerNextY(const ODEFunc& f, const double x, const double y, const double h);
GridFunc AdamsODE(const ODEFunc& f, const double y0, const double a, const double b, const size_t n);
GridFunc EulerAdamsPredictorCorrectorODE(const ODEFunc& f, const double y0, const double a, const double b, const size_t n);
double MaxError(const MathFunc& f, const GridFunc& gridF);

void WriteGridFunc(const string& filename, const GridFunc& gridFunc);
void WriteODESolutions(const string& filenameGeneric, const ODEMethod& Method, const ODEFunc& f, const double y0, const double a, const double b, const vector<size_t>& ns);
void WriteErrorOnH(const string& filename, const ODEMethod& Method, const MathFunc& sol, const ODEFunc& f,
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
    const size_t nInit = 1;
    const size_t nSteps = 15;

    cout << "Adams" << endl;
    GridFunc ans = AdamsODE(ode, y0, a, b, 10);
    for (auto& point : ans)
        cout << point.x << "\t" << point.y << "\t" << y(point.x) << endl;

    cout << "Euler predictor - Adams corrector" << endl;
    ans = EulerAdamsPredictorCorrectorODE(ode, y0, a, b, 10);
    for (auto& point : ans)
        cout << point.x << "\t" << point.y << "\t" << y(point.x) << endl;

    try {
        WriteODESolutions(ROOT"csv/ans", AdamsODE, ode, y0, a, b, ns);
        WriteErrorOnH(ROOT"csv/err_on_h.csv", AdamsODE, y, ode, y0, a, b, nInit, nSteps);
        WriteODESolutions(ROOT"csv/anspc", EulerAdamsPredictorCorrectorODE, ode, y0, a, b, ns);
        WriteErrorOnH(ROOT"csv/err_on_hpc.csv", EulerAdamsPredictorCorrectorODE, y, ode, y0, a, b, nInit, nSteps);
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

double EulerNextY(const ODEFunc& f, const double x, const double y, const double h) {
    return y + h * f(x + h / 2.0, y + h / 2.0 * f(x, y));
}

GridFunc AdamsODE(const ODEFunc& f, const double y0, const double a, const double b, const size_t n) {
    const double h = (b - a) / (double)n;
    GridFunc res = { Point(a, y0) };
    double xPrev = a;
    double yPrev = y0;
    double xCur = a + h;
    double yCur = EulerNextY(f, xPrev, yPrev, h);
    double fPrev = f(xPrev, yPrev);
    res.push_back(Point(xCur, yCur));
    for (size_t i = 1; i < n; i++) {
        double fCur = f(xCur, yCur);
        double yNext = yCur  + h / 2.0 * (3 * fCur - fPrev);
        double xNext = xCur + h;
        yPrev = yCur;
        xPrev = xCur;
        fPrev = fCur;
        yCur = yNext;
        xCur = xNext;
        res.push_back(Point(xCur, yCur));
    }
    return res;
}

GridFunc EulerAdamsPredictorCorrectorODE(const ODEFunc& f, const double y0, const double a, const double b, const size_t n) {
    const double h = (b - a) / (double)n;
    GridFunc res = { Point(a, y0) };
    double xCur = a;
    double yCur = y0;
    for (size_t i = 0; i < n; i++) {
        double xNext = xCur + h;
        // mod euler predictor
        double yNextPredict = EulerNextY(f, xCur, yCur, h);
        // corrector
        double yNext = yCur + h / 2.0 * (f(xNext, yNextPredict) + f(xCur, yCur));
        yCur = yNext;
        xCur = xNext;
        res.push_back(Point(xCur, yCur));
    }
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
                       const ODEMethod& Method,
                       const ODEFunc& f,
                       const double y0,
                       const double a,
                       const double b,
                       const vector<size_t>& ns) {
    for (size_t i = 0; i < ns.size(); i++)
        WriteGridFunc(filenameGeneric + to_string(i + 1) + ".csv", Method(f, y0, a, b, ns[i] - 1));
}

void WriteErrorOnH(const string& filename,
                   const ODEMethod& Method,
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
        file << (b - a) / (double)n << ";" << MaxError(sol, Method(f, y0, a, b, n)) << endl;
        n *= 2;
    }
    file.close();
}