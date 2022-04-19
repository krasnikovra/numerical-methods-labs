#include <functional>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <algorithm>
#include <string>
#include <cassert>
#include <cmath>

#define ROOT "../../"
#define PRECISION 50
#define SET_STREAM_PRECISION(stream) \
    (stream).setf(ios::fixed); \
    (stream) << setprecision(PRECISION)

using namespace std;

using Vector = vector<double>;

struct Point {
    double x;
    Vector y;
    Point(const double x_, const Vector& y_) : x(x_), y(y_) {}
};

using GridFunc = vector<Point>;
using MathFunc = function<double(const double)>;
using ODEFunc = function<Vector(const double, const Vector&)>;

struct SecOrderODE {
    MathFunc p, q, r, f;
    SecOrderODE(const MathFunc& p_, const MathFunc& q_, const MathFunc& r_, const MathFunc& f_) :
        p(p_), q(q_), r(r_), f(f_) {}
};

double Pow(const double a, const size_t b);
GridFunc ModEulerODE(const ODEFunc& f, const Vector y0, const double a, const double b,
                     const size_t n);
ODEFunc MakeODEFuncOf2ndOrderODE(const SecOrderODE& ode);
MathFunc MakeZeroMathFunc();
GridFunc SolveBordaryValueProblem(const SecOrderODE& ode, const double a, const double b,
                                  const Vector& alpha, const Vector& beta, const double A,
                                  const double B, const size_t n);
double MaxError(const MathFunc& f, const GridFunc& gridFunc, const size_t derivativeOrder = 0);

void WriteGridFunc(const string& filename, const GridFunc& gridFunc, const size_t derivativeOrder = 0);
void WriteSolutions(const string& filenameGeneric, const SecOrderODE& ode, const double a, const double b,
                    const Vector& alpha, const Vector& beta, const double A, const double B,
                    const vector<size_t>& ns, const size_t derivativeOrder = 0);
void WriteErrorOnDA(const string& filename, const MathFunc& sol, const SecOrderODE& ode,
                    const double a, const double b, const Vector& alpha,
                    const Vector& beta, const double A, const double B, 
                    const size_t n, const double dA0, const size_t dASteps,
                    const size_t derivativeOrder = 0);
void WriteErrorOnDB(const string& filename, const MathFunc& sol, const SecOrderODE& ode,
                    const double a, const double b, const Vector& alpha,
                    const Vector& beta, const double A, const double B,
                    const size_t n, const double dB0, const size_t dBSteps,
                    const size_t derivativeOrder = 0);
void WriteErrorOnH(const string& filename, const MathFunc& sol, const SecOrderODE& ode,
                   const double a, const double b, const Vector& alpha,
                   const Vector& beta, const double A, const double B,
                   const size_t n0, const size_t nSteps, const size_t derivativeOrder = 0);

int main() {
    const MathFunc solEx = [](const double x) -> double { return exp(x); };
    const MathFunc solExDer = [](const double x) -> double { return exp(x); };

    const MathFunc p = [](const double x) -> double { return 1.0; };
    const MathFunc q = [](const double x) -> double { return 1.0 + Pow(sin(x), 2); };
    const MathFunc r = [](const double x) -> double { return Pow(cos(x), 2); };
    const MathFunc f = [](const double x) -> double { return 3.0 * exp(x); };
    const double a = 0.0;
    const double b = 1.0;
    const Vector alpha = { 1.0, 1.0 };
    const Vector beta = { 1.0, 0.0 };
    const double A = alpha[0] * solEx(a) + alpha[1] * solExDer(a);
    const double B = beta[0] * solEx(b) + beta[1] * solExDer(b);
    const size_t nTest = 10;

    //----------------------- test --------------------------
    const SecOrderODE ode(p, q, r, f);
    GridFunc sol = SolveBordaryValueProblem(ode, a, b, alpha, beta, A, B, nTest);

    cout << "x" << "\ty" << "\te^x" << endl;
    for (auto& point : sol)
        cout << point.x << "\t" << point.y[0] << "\t" << solEx(point.x) << endl;
    //--------------------------------------------------------

    //--------------------- research -------------------------
    const vector<size_t> ns = { 5, 10 };
    const size_t n0 = 1;
    const size_t nSteps = 15;
    const size_t n = 10;
    const double dA0 = A;
    const size_t dASteps = 10;
    const double dB0 = B;
    const size_t dBSteps = 10;

    try {
        WriteSolutions(ROOT"csv/ans", ode, a, b, alpha, beta, A, B, ns);
        WriteErrorOnH(ROOT"csv/err_on_h.csv", solEx, ode, a, b, alpha, beta, A, B, n0, nSteps);
        WriteErrorOnDA(ROOT"csv/err_on_da.csv", solEx, ode, a, b, alpha, beta, A, B, n, dA0, dASteps);
        WriteErrorOnDB(ROOT"csv/err_on_db.csv", solEx, ode, a, b, alpha, beta, A, B, n, dB0, dBSteps);
    }
    catch (const exception& e) {
        cout << "An error occured!" << endl;
        cout << e.what() << endl;
    }

    cout << "Success!" << endl;
    //--------------------------------------------------------

    return 0;
}

Vector operator*(const double a, const Vector& vec) {
    Vector res;
    for (auto& x : vec)
        res.push_back(x * a);
    return res;
}

Vector operator+(const Vector& a, const Vector& b) {
    if (a.size() != b.size())
        throw exception("Vectors adding should have same size");
    Vector res;
    for (size_t i = 0; i < a.size(); i++)
        res.push_back(a[i] + b[i]);
    return res;
}

GridFunc operator*(const double a, const GridFunc& f) {
    GridFunc res;
    for (auto& point : f)
        res.push_back(Point(point.x, a * point.y));
    return res;
}

GridFunc operator+(const GridFunc& f, const GridFunc& g) {
    GridFunc res;
    for (size_t i = 0; i < f.size(); i++) {
        if (f[i].x != g[i].x)
            throw exception("GridFuncs adding should have same x grid");
        res.push_back(Point(f[i].x, f[i].y + g[i].y));
    }
    return res;
}

double Pow(const double a, const size_t b) {
    double res = 1;
    for (size_t i = 0; i < b; i++)
        res *= a;
    return res;
}

GridFunc ModEulerODE(const ODEFunc& f, const Vector y0,
                     const double a, const double b, const size_t n) {
    const double h = (b - a) / (double)n;
    double x = a;
    Vector y = y0;
    GridFunc res = { Point(a, y)};
    for (size_t i = 0; i < n; i++) {
        y = y + h * f(x + h / 2.0, y + h / 2.0 * f(x, y));
        x = x + h;
        res.push_back(Point(x, y));
    }
    return res;
}

ODEFunc MakeODEFuncOf2ndOrderODE(const SecOrderODE& ode) {
    return [ode](const double x, const Vector& y) -> Vector { 
        return { y[1], ode.f(x) - ode.r(x) / ode.p(x) * y[0] - ode.q(x) / ode.p(x) * y[1] }; 
    };
}

MathFunc MakeZeroMathFunc() {
    return [](const double x) -> double { return 0.0; };
}

GridFunc SolveBordaryValueProblem(const SecOrderODE& ode,
                                  const double a,
                                  const double b,
                                  const Vector& alpha,
                                  const Vector& beta,
                                  const double A,
                                  const double B,
                                  const size_t n) {
    const Vector u0 = {
        alpha[0] / (Pow(alpha[0], 2) + Pow(alpha[1], 2)) * A,
        alpha[1] / (Pow(alpha[0], 2) + Pow(alpha[1], 2)) * A,
    };
    const Vector v0 = {
        alpha[1],
        -alpha[0]
    };
    GridFunc u = ModEulerODE(
        MakeODEFuncOf2ndOrderODE(ode),
        u0, a, b, n
    );
    GridFunc v = ModEulerODE(
        MakeODEFuncOf2ndOrderODE(SecOrderODE(ode.p, ode.q, ode.r, MakeZeroMathFunc())),
        u0, a, b, n
    );
    const double c = (B - beta[0] * u[u.size() - 1].y[0] - beta[1] * u[u.size() - 1].y[1]) /
        (beta[0] * v[v.size() - 1].y[0] + beta[1] * v[v.size() - 1].y[1]);
    return u + c * v;
}

double MaxError(const MathFunc& f,
                const GridFunc& gridFunc,
                const size_t derivativeOrder) {
    double res = 0;
    for (auto& p : gridFunc) {
        double err = abs(f(p.x) - p.y[derivativeOrder]);
        if (err > res)
            res = err;
    }
    return res;
}

void WriteGridFunc(const string& filename,
                   const GridFunc& gridFunc,
                   const size_t derivativeOrder) {
    ofstream file(filename);
    if (!file.is_open())
        throw runtime_error(
            string("File ") + filename + string(" could not be opened.")
        );
    SET_STREAM_PRECISION(file);
    file << "x;y" << endl;
    for (auto& point : gridFunc)
        file << point.x << ";" << point.y[derivativeOrder] << endl;
    file.close();
}


void WriteSolutions(const string& filenameGeneric,
                    const SecOrderODE& ode,
                    const double a,
                    const double b,
                    const Vector& alpha,
                    const Vector& beta,
                    const double A,
                    const double B,
                    const vector<size_t>& ns,
                    const size_t derivativeOrder) {
    for (size_t i = 0; i < ns.size(); i++)
        WriteGridFunc(
            filenameGeneric + to_string(i + 1) + ".csv",
            SolveBordaryValueProblem(ode, a, b, alpha, beta, A, B, ns[i]),
            derivativeOrder
        );
}

void WriteErrorOnDA(const string& filename,
                    const MathFunc& sol,
                    const SecOrderODE& ode,
                    const double a,
                    const double b,
                    const Vector& alpha,
                    const Vector& beta,
                    const double A,
                    const double B,
                    const size_t n,
                    const double dA0,
                    const size_t dASteps,
                    const size_t derivativeOrder) {
    ofstream file(filename);
    if (!file.is_open())
        throw string("File ") + filename + string(" could not be opened.");
    SET_STREAM_PRECISION(file);
    file << "dy;err;" << endl;
    double dA = dA0;
    for (size_t i = 0; i < dASteps; i++) {
        file << dA << ";" <<
            MaxError(
                sol,
                SolveBordaryValueProblem(ode, a, b, alpha, beta, A + dA, B, n),
                derivativeOrder
            ) << endl;
        dA /= 10;
    }
    file.close();
}

void WriteErrorOnDB(const string& filename,
                    const MathFunc& sol,
                    const SecOrderODE& ode,
                    const double a,
                    const double b,
                    const Vector& alpha,
                    const Vector& beta,
                    const double A,
                    const double B,
                    const size_t n,
                    const double dB0,
                    const size_t dBSteps,
                    const size_t derivativeOrder) {
    ofstream file(filename);
    if (!file.is_open())
        throw string("File ") + filename + string(" could not be opened.");
    SET_STREAM_PRECISION(file);
    file << "dy;err;" << endl;
    double dB = dB0;
    for (size_t i = 0; i < dBSteps; i++) {
        file << dB << ";" <<
            MaxError(
                sol,
                SolveBordaryValueProblem(ode, a, b, alpha, beta, A, B + dB, n),
                derivativeOrder
            ) << endl;
        dB /= 10;
    }
    file.close();
}

void WriteErrorOnH(const string& filename,
                   const MathFunc& sol,
                   const SecOrderODE& ode,
                   const double a,
                   const double b,
                   const Vector& alpha,
                   const Vector& beta,
                   const double A,
                   const double B,
                   const size_t n0,
                   const size_t nSteps,
                   const size_t derivativeOrder) {
    ofstream file(filename);
    if (!file.is_open())
        throw string("File ") + filename + string(" could not be opened.");
    SET_STREAM_PRECISION(file);
    file << "h;err;" << endl;
    size_t n = n0;
    for (size_t i = 0; i < nSteps; i++) {
        file << (b - a) / (double)n << ";" <<
            MaxError(
                sol, 
                SolveBordaryValueProblem(ode, a, b, alpha, beta, A, B, n), 
                derivativeOrder
            ) << endl;
        n *= 2;
    }
    file.close();
}
