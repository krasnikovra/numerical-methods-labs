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

Grid MakeUniformGrid(double a, double b, size_t n);
Data MakeDataOutOfFunc(Grid grid, MathFunc f, MathFunc fDer);
double EvaluateHermitePolynomial(double x, Data data);
double EvaluateMaxMidpointsError(Data data, MathFunc f);
double EvaluateRootPolynomial(double x, Grid grid);
size_t Factorial(size_t n);
double EvaluateTheoreticError(double x, MathFunc fNPlusOneDer, Grid grid);

void WritePolynomialsValues(const char* filenameVal, const char* filenameGrid, double a, double b, size_t pointsDensity, vector<size_t> ns, MathFunc f, MathFunc fDer);
void WriteMaxMidpointErrorOnN(const char* filenameErr, double a, double b, size_t n1, size_t n2, MathFunc f, MathFunc fDer);
void WriteTheoreticError(const char* filenameErr, size_t pointsDensity, Grid grid, MathFunc fNthDer);

int main() {
	/*
	Data testData = {
		{0,1,1},
		{1,3,2},
		{2,7,5},
	};
	vector<double> testCoef = EvaluateHermitePolynomialCoefficients(testData);
	*/
	const double a = -3.0;
	const double b = -1.0;
	const size_t n = 3;
	const size_t n1 = 2;
	const size_t n2 = 38;
	const size_t pointsDensity = 1000;
	const MathFunc f = [](double x) { return 1 / tan(x) - x; };
	const MathFunc fDer = [](double x) { return -1 / (sin(x) * sin(x)) - 1; };
	const MathFunc f3thDer = [](double x) { return -2 / (sin(x) * sin(x)) * (2 / (tan(x) * tan(x)) + 1 / (sin(x) * sin(x))); };
	
	WritePolynomialsValues(ROOT"csv/values.csv", ROOT"csv/grids.csv", a, b, pointsDensity, { 3, 4, 5 }, f, fDer);
	WriteMaxMidpointErrorOnN(ROOT"csv/error_on_n.csv", a, b, n1, n2, f, fDer);
	WriteTheoreticError(ROOT"csv/error_theoretic.csv", pointsDensity, MakeUniformGrid(a, b, 3), f3thDer);

	return 0;
}

Grid MakeUniformGrid(double a, double b, size_t n) {
	Grid res;
	for (size_t i = 0; i < n; i++)
		res.push_back(a + (b - a) * i / (double)(n - 1));
	return res;
}

Data MakeDataOutOfFunc(Grid grid, MathFunc f, MathFunc fDer) {
	Data res;
	for (auto& x : grid)
		res.push_back(DataCell(x, f(x), fDer(x)));
	return res;
}

double EvaluateHermitePolynomial(double x, Data data) {
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

double EvaluateMaxMidpointsError(Data data, MathFunc f) {
	double res = 0;
	for (size_t i = 0; i + 1 < data.size(); i++) {
		double x = (data[i].x + data[i + 1].x) / 2.0;
		double err = abs(EvaluateHermitePolynomial(x, data) - f(x));
		if (err > res)
			res = err;
	}
	return res;
}

double EvaluateRootPolynomial(double x, Grid grid) {
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

double EvaluateTheoreticError(double x, MathFunc fNthDer, Grid grid) {
	return abs(fNthDer(x) * EvaluateRootPolynomial(x, grid) / (double)Factorial(grid.size()));
}

void WritePolynomialsValues(const char* filenameVal, const char* filenameGrid, double a, double b, size_t pointsDensity, vector<size_t> ns, MathFunc f, MathFunc fDer) {
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
		Grid grid = MakeUniformGrid(a, b, n);
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

void WriteMaxMidpointErrorOnN(const char* filenameErr, double a, double b, size_t n1, size_t n2, MathFunc f, MathFunc fDer) {
	ofstream fileErr(filenameErr);
	if (!fileErr.is_open())
		throw "Error while opening the file";
	SET_STREAM_PRECISION(fileErr);

	for (size_t i = n1; i <= n2; i++)
		fileErr << i << ";";
	fileErr << endl;
	for (size_t i = n1; i <= n2; i++) {
		Grid grid = MakeUniformGrid(a, b, i);
		Data data = MakeDataOutOfFunc(grid, f, fDer);
		fileErr << EvaluateMaxMidpointsError(data, f) << ";";
	}
	fileErr << endl;
	fileErr.close();
}

void WriteTheoreticError(const char* filenameErr, size_t pointsDensity, Grid grid, MathFunc fNthDer) {
	ofstream fileErr(filenameErr);
	if (!fileErr.is_open())
		throw "Error while opening the file";
	SET_STREAM_PRECISION(fileErr);
	
	fileErr << grid.size() << endl;
	for (auto& x : MakeUniformGrid(grid[0], grid[grid.size() - 1], pointsDensity))
		fileErr << EvaluateTheoreticError(x, fNthDer, grid) << ";";
	fileErr << endl;
	fileErr.close();
}