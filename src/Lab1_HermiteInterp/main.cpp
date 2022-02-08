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

struct DataCell {
	double x, y, der;
	DataCell(double _x, double _y, double _der) : x(_x), y(_y), der(_der) {};
};

using Data = vector<DataCell>;
using Grid = vector<double>;
using MathFunc = double (*)(double);

Grid MakeUniformGrid(double a, double b, size_t n, double alpha = 1.0);
Data MakeDataOutOfFunc(const Grid& grid, MathFunc f, MathFunc fDer);
double EvaluateHermitePolynomial(double x, const Data& data);
double EvaluateMaxMidpointsError(const Data& data, MathFunc f);
double EvaluateRootPolynomial(double x, const Grid& grid);
size_t Factorial(size_t n); 
double EvaluateTheoreticError(double x, double etha, MathFunc fNPlusOneDer, Grid grid);

void WritePolynomialsValues(const char* filenameVal, const char* filenameGrid, double a, double b, size_t pointsDensity, const vector<size_t>& ns, MathFunc f, MathFunc fDer, double alpha = 1.0);
void WriteMaxMidpointErrorOnN(const char* filenameErr, double a, double b, size_t n1, size_t n2, MathFunc f, MathFunc fDer, double alpha = 1.0);
void WriteMaxMidpointErrorOnAlpha(const char* filenameErr, double a, double b, size_t n, double alpha1, double alpha2, double step, MathFunc f, MathFunc fDer);
void WriteTheoreticError(const char* filenameErr, double etha, size_t pointsDensity, const Grid& grid, MathFunc fNthDer);

int main() {
	/*
	Data testData = {
		{0,1,1},
		{1,3,2},
		{2,7,5},
	};
	*/
	const double a = -3.0;
	const double b = -1.0;
	const size_t n = 3;
	const size_t n1 = 2;
	const size_t n2 = 38;
	const size_t pointsDensity = 1000;
	const MathFunc f = [](double x) { return 1 / tan(x) - x; };
	const MathFunc fDer = [](double x) { return -1 / (sin(x) * sin(x)) - 1; };
	const MathFunc f6thDer = [](double x) { return 4 * (56 * cos(2 * x) + cos(4 * x) + 123) / (tan(x) * pow(sin(x), 6)); };
	
	WritePolynomialsValues(ROOT"csv/values.csv", ROOT"csv/grids.csv", a, b, pointsDensity, { 3, 4, 5 }, f, fDer);
	WriteMaxMidpointErrorOnN(ROOT"csv/error_on_n.csv", a, b, n1, n2, f, fDer);
	WriteTheoreticError(ROOT"csv/error_theoretic.csv", (a + b)/ 5, pointsDensity, MakeUniformGrid(a, b, 3), f6thDer);
	
	const double alpha = 1.23;

	WritePolynomialsValues(ROOT"csv/valuesparam.csv", ROOT"csv/gridsparam.csv", a, b, pointsDensity, { 3, 4, 5 }, f, fDer, alpha);
	WriteMaxMidpointErrorOnN(ROOT"csv/error_on_nparam.csv", a, b, n1, 35, f, fDer, alpha);
	WriteMaxMidpointErrorOnAlpha(ROOT"csv/error_on_alpha.csv", a, b, 3, 0.1, 3, 0.01, f, fDer);

	return 0;
}

Grid MakeUniformGrid(double a, double b, size_t n, double alpha) {
	Grid res;
	for (size_t i = 0; i < n; i++)
		res.push_back(a + (b - a) * pow(i / (double)(n - 1), alpha));
	return res;
}

Data MakeDataOutOfFunc(const Grid& grid, MathFunc f, MathFunc fDer) {
	Data res;
	for (auto& x : grid)
		res.push_back(DataCell(x, f(x), fDer(x)));
	return res;
}

double EvaluateHermitePolynomial(double x, const Data& data) {
	double res = 0;
	for (size_t j = 0; j < data.size(); j++) {
		double prod = 1, sum = 0;
		for (size_t k = 0; k < data.size(); k++) {
			if (k == j)
				continue;
			prod *= (x - data[k].x) / (data[j].x - data[k].x);
			sum += (x - data[j].x) / (data[j].x - data[k].x);
		}
		prod = prod * prod;
		res += ((x - data[j].x) * data[j].der + (1 - 2 * sum) * data[j].y) * prod;
	}
	return res;
}

double EvaluateMaxMidpointsError(const Data& data, MathFunc f) {
	double res = 0;
	for (size_t i = 0; i + 1 < data.size(); i++) {
		double x = (data[i].x + data[i + 1].x) / 2.0;
		double err = abs(EvaluateHermitePolynomial(x, data) - f(x));
		if (err > res)
			res = err;
	}
	return res;
}

double EvaluateRootPolynomial(double x, const Grid& grid) {
	double res = 1;
	for (auto& xi : grid)
		res *= (x - xi);
	return res;
}

size_t Factorial(size_t n) {
	size_t res = 1;
	for (size_t i = 1; i <= n; i++)
		res *= i;
	return res;
}

double EvaluateTheoreticError(double x, double etha, MathFunc f2NthDer, Grid grid) {
	return abs(f2NthDer(etha) * pow(EvaluateRootPolynomial(x, grid), 2) / (double)Factorial(2 * grid.size()));
}

void WritePolynomialsValues(const char* filenameVal, const char* filenameGrid, double a, double b, size_t pointsDensity, const vector<size_t>& ns, MathFunc f, MathFunc fDer, double alpha) {
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
		Data data = MakeDataOutOfFunc(grid, f, fDer);
		for (auto& x : xGrid)
			fileVal << EvaluateHermitePolynomial(x, data) << ";";
		fileVal << endl;
	}
	fileVal.close();
	fileGrid.close();
}

void WriteMaxMidpointErrorOnN(const char* filenameErr, double a, double b, size_t n1, size_t n2, MathFunc f, MathFunc fDer, double alpha) {
	ofstream fileErr(filenameErr);
	if (!fileErr.is_open())
		throw "Error while opening the file";
	SET_STREAM_PRECISION(fileErr);

	for (size_t i = n1; i <= n2; i++)
		fileErr << i << ";";
	fileErr << endl;
	for (size_t i = n1; i <= n2; i++) {
		Grid grid = MakeUniformGrid(a, b, i, alpha);
		Data data = MakeDataOutOfFunc(grid, f, fDer);
		fileErr << EvaluateMaxMidpointsError(data, f) << ";";
	}
	fileErr << endl;
	fileErr.close();
}

void WriteMaxMidpointErrorOnAlpha(const char* filenameErr, double a, double b, size_t n, double alpha1, double alpha2, double step, MathFunc f, MathFunc fDer) {
	ofstream fileErr(filenameErr);
	if (!fileErr.is_open())
		throw "Error while opening the file";
	SET_STREAM_PRECISION(fileErr);

	for (double alpha = alpha1; alpha <= alpha2; alpha += step)
		fileErr << alpha << ";";
	fileErr << endl;
	for (double alpha = alpha1; alpha <= alpha2; alpha += step) {
		Grid grid = MakeUniformGrid(a, b, n, alpha);
		Data data = MakeDataOutOfFunc(grid, f, fDer);
		fileErr << EvaluateMaxMidpointsError(data, f) << ";";
	}
	fileErr << endl;
	fileErr.close();
}

void WriteTheoreticError(const char* filenameErr, double etha, size_t pointsDensity, const Grid& grid, MathFunc fNthDer) {
	ofstream fileErr(filenameErr);
	if (!fileErr.is_open())
		throw "Error while opening the file";
	SET_STREAM_PRECISION(fileErr);
	
	fileErr << grid.size() << endl;
	for (auto& x : MakeUniformGrid(grid[0], grid[grid.size() - 1], pointsDensity))
		fileErr << EvaluateTheoreticError(x, etha, fNthDer, grid) << ";";
	fileErr << endl;
	fileErr.close();
}