#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <cmath>

#define ROOT "../../"
#define PRECISION 50
#define SET_STREAM_PRECISION(stream) \
    (stream).setf(ios::fixed); \
    (stream) << setprecision(PRECISION)

using namespace std;

struct TridiagonalMatrix {
	vector<double> b, c, d;
};

struct Data {
	vector<double> x, y;
};

using Grid = vector<double>;
using Vector = vector<double>;
using MathFunc = double (*)(double);

Grid MakeUniformGrid(double a, double b, size_t n, double alpha = 1.0);
Data MakeDataOutOfFunc(const Grid& grid, MathFunc f);
TridiagonalMatrix MakeTridiagonalMatrixOutOfData(const Data& data);
Vector GetRightPartVector(const Data& data);
Vector ThomasAlghorithm(const TridiagonalMatrix& A, const Vector& r);
Vector GetSplineMParams(const Data& data);
double EvaluateIthSplinePoly(double x, const Vector& M, const Data& data, const size_t i);
double EvaluateSpline(double x, const Vector& M, const Data& data);

void WriteSplinesValues(const char* filenameVal, const char* filenameGrid, double a, double b, size_t pointsDensity, const vector<size_t>& ns, MathFunc f, double alpha = 1.0);
void WriteMaxMidpointErrorOnN(const char* filenameErr, double a, double b, size_t n1, size_t n2, MathFunc f, double alpha = 1.0);
void WriteMaxMidpointErrorOnAlpha(const char* filenameErr, double a, double b, size_t n, double alpha1, double alpha2, double step, MathFunc f);

int main(void) {
	MathFunc f = [](double x) { return 1 / tan(x) - x; };
	const double a = -3;
	const double b = -1;
	const size_t plotPointsCount = 1000;
	const vector<size_t> n{ 4,5,6 };
	const size_t n1 = 4;
	const size_t n2 = 100;
	const size_t nForAlphaExp = 5;

	WriteSplinesValues(ROOT"csv/values.csv", ROOT"csv/grids.csv", a, b, plotPointsCount, n, f);
	WriteMaxMidpointErrorOnN(ROOT"csv/error_on_n.csv", a, b, n1, n2, f);

	WriteMaxMidpointErrorOnAlpha(ROOT"csv/error_on_alpha.csv", a, b, nForAlphaExp, 0.01, 3, 1e-3, f);

	const double alpha = 2;

	WriteSplinesValues(ROOT"csv/valuesparam.csv", ROOT"csv/gridsparam.csv", a, b, plotPointsCount, n, f, alpha);
	WriteMaxMidpointErrorOnN(ROOT"csv/error_on_nparam.csv", a, b, n1, n2, f, alpha);

	return 0;
}

Grid MakeUniformGrid(double a, double b, size_t n, double alpha) {
	Grid res;
	for (size_t i = 0; i < n; i++)
		res.push_back(a + (b - a) * pow(i / (double)(n - 1), alpha));
	return res;
}

Data MakeDataOutOfFunc(const Grid& grid, MathFunc f) {
	Data res;
	for (auto& x : grid) {
		res.x.push_back(x);
		res.y.push_back(f(x));
	}
	return res;
}

TridiagonalMatrix MakeTridiagonalMatrixOutOfData(const Data& data) {
	TridiagonalMatrix res;
	const size_t n = data.x.size() - 1;
	double h1 = data.x[1] - data.x[0];
	double h2 = data.x[2] - data.x[1];
	res.b.push_back(0);
	res.c.push_back(2 + h1 / (h1 + h2) + (h1 * h1) / ((h1 + h2) * h2)); // from not-a-knot spline condition
	res.d.push_back(h2 / (h1 + h2) - (h1 * h1) / ((h1 + h2) * h2)); // from not-a-knot spline condition
	for (size_t i = 2; i < n - 1; i++) {
		double hi = data.x[i] - data.x[i - 1];
		double hip1 = data.x[i + 1] - data.x[i];
		res.b.push_back(hi / (hi + hip1));
		res.c.push_back(2);
		res.d.push_back(hip1 / (hi + hip1));
	}
	double hnm1 = data.x[n - 1] - data.x[n - 2];
	double hn = data.x[n] - data.x[n - 1];
	res.b.push_back(hnm1 / (hn + hnm1) - (hn * hn) / ((hn + hnm1) * hnm1));
	res.c.push_back(2 + hn / (hnm1 + hn) + (hn * hn) / ((hn + hnm1) * hnm1));
	res.d.push_back(0);
	return res;
}

Vector GetRightPartVector(const Data& data) {
	Vector res;
	const size_t n = data.x.size() - 1;
	for (size_t i = 1; i < n; i++) {
		double hi = data.x[i] - data.x[i - 1];
		double hip1 = data.x[i + 1] - data.x[i];
		double ci = 6 / (hi + hip1) * ((data.y[i + 1] - data.y[i]) / hip1 - (data.y[i] - data.y[i - 1]) / hi);
		res.push_back(ci);
	}
	return res;
}

Vector ThomasAlghorithm(const TridiagonalMatrix& A, const Vector& r) {
	const size_t n = r.size() - 1;
	Vector deltas, lambdas;
	deltas.push_back(-A.d[0] / A.c[0]);
	lambdas.push_back(r[0] / A.c[0]);
	for (size_t i = 1; i < n; i++) {
		deltas.push_back(-A.d[i] / (A.b[i] * deltas[i - 1] + A.c[i]));
		lambdas.push_back((r[i] - A.b[i] * lambdas[i - 1] ) / (A.b[i] * deltas[i - 1] + A.c[i]));
	}
	deltas.push_back(0);
	lambdas.push_back((r[n] - A.b[n] * lambdas[n - 1]) / (A.b[n] * deltas[n - 1] + A.c[n]));

	Vector res;
	res.resize(n + 1);
	res[n] = lambdas[n];
	for (size_t i = n - 1; i >= 0; i--) {
		res[i] = deltas[i] * res[i + 1] + lambdas[i];
		if (i == 0)
			break;
	}
	return res;
}

Vector GetSplineMParams(const Data& data) {
	TridiagonalMatrix A = MakeTridiagonalMatrixOutOfData(data);
	Vector r = GetRightPartVector(data);
	Vector secDerParams = ThomasAlghorithm(A, r);
	const size_t n = data.x.size() - 1;
	double h1 = data.x[1] - data.x[0];
	double h2 = data.x[2] - data.x[1];
	double M0 = (1 + h1 / h2) * secDerParams[0] - h1 / h2 * secDerParams[1];
	double hnm1 = data.x[n - 1] - data.x[n - 2];
	double hn = data.x[n] - data.x[n - 1];
	double Mn = (1 + hn / hnm1) * secDerParams[secDerParams.size() - 1] - hn / hnm1 * secDerParams[secDerParams.size() - 2];
	secDerParams.insert(secDerParams.begin(), M0);
	secDerParams.insert(secDerParams.end(), Mn);
	return secDerParams;
}

double EvaluateIthSplinePoly(double x, const Vector& M, const Data& data, const size_t i) {
	return M[i - 1] * (data.x[i] - x) * (data.x[i] - x) * (data.x[i] - x) / (6 * (data.x[i] - data.x[i - 1])) +
		M[i] * (x - data.x[i - 1]) * (x - data.x[i - 1]) * (x - data.x[i - 1]) / (6 * (data.x[i] - data.x[i - 1])) +
		((data.y[i] - data.y[i - 1]) / (data.x[i] - data.x[i - 1]) - (data.x[i] - data.x[i - 1]) / 6 * (M[i] - M[i - 1])) * (x - data.x[i - 1]) +
		data.y[i - 1] - M[i - 1] * (data.x[i] - data.x[i - 1]) * (data.x[i] - data.x[i - 1]) / 6;
}

double EvaluateSpline(double x, const Vector& M, const Data& data) {
	for (size_t i = 1; i < data.x.size(); i++) {
		if (x <= data.x[i])
			return EvaluateIthSplinePoly(x, M, data, i);
	}
	return nan("nan");
}

double EvaluateMaxMidpointsError(const Data& data, const Vector& M, MathFunc f) {
	double res = 0;
	for (size_t i = 0; i + 1 < data.x.size(); i++) {
		double x = (data.x[i] + data.x[i + 1]) / 2.0;
		double err = abs(EvaluateSpline(x, M, data) - f(x));
		if (err > res)
			res = err;
	}
	return res;
}

void WriteSplinesValues(const char* filenameVal, const char* filenameGrid, double a, double b, size_t pointsDensity, const vector<size_t>& ns, MathFunc f, double alpha) {
	ofstream fileVal(filenameVal);
	if (!fileVal.is_open())
		throw "Error while opening the file";
	ofstream fileGrid(filenameGrid);
	if (!fileGrid.is_open())
		throw "Error while opening the file";
	SET_STREAM_PRECISION(fileVal);
	SET_STREAM_PRECISION(fileGrid);
	Grid xGrid = MakeUniformGrid(a, b, pointsDensity);
	for (auto& n : ns) {
		fileVal << n << ";";
		fileGrid << n << ";";
	}
	fileVal << endl;
	fileGrid << endl;
	for (auto& x : xGrid)
		fileVal << x << ";";
	fileVal << endl;
	for (auto& n : ns) {
		Grid grid = MakeUniformGrid(a, b, n, alpha);
		for (auto& x : grid)
			fileGrid << x << ";";
		fileGrid << endl;
		Data data = MakeDataOutOfFunc(grid, f);
		Vector M = GetSplineMParams(data);
		for (auto& x : xGrid)
			fileVal << EvaluateSpline(x, M, data) << ";";
		fileVal << endl;
	}
	fileVal.close();
	fileGrid.close();
}

void WriteMaxMidpointErrorOnN(const char* filenameErr, double a, double b, size_t n1, size_t n2, MathFunc f, double alpha) {
	ofstream fileErr(filenameErr);
	if (!fileErr.is_open())
		throw "Error while opening the file";
	SET_STREAM_PRECISION(fileErr);

	for (size_t i = n1; i <= n2; i++)
		fileErr << i << ";";
	fileErr << endl;
	for (size_t i = n1; i <= n2; i++) {
		Grid grid = MakeUniformGrid(a, b, i, alpha);
		Data data = MakeDataOutOfFunc(grid, f);
		Vector M = GetSplineMParams(data);
		fileErr << EvaluateMaxMidpointsError(data, M, f) << ";";
	}
	fileErr << endl;
	fileErr.close();
}

void WriteMaxMidpointErrorOnAlpha(const char* filenameErr, double a, double b, size_t n, double alpha1, double alpha2, double step, MathFunc f) {
	ofstream fileErr(filenameErr);
	if (!fileErr.is_open())
		throw "Error while opening the file";
	SET_STREAM_PRECISION(fileErr);

	for (double alpha = alpha1; alpha <= alpha2; alpha += step)
		fileErr << alpha << ";";
	fileErr << endl;
	for (double alpha = alpha1; alpha <= alpha2; alpha += step) {
		Grid grid = MakeUniformGrid(a, b, n, alpha);
		Data data = MakeDataOutOfFunc(grid, f);
		Vector M = GetSplineMParams(data);
		fileErr << EvaluateMaxMidpointsError(data, M, f) << ";";
	}
	fileErr << endl;
	fileErr.close();
}