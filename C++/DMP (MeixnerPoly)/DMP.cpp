#include <iostream>
#include <vector>
#include <symengine/symbol.h>
#include <symengine/matrix.h>
#include <symengine/lambda_double.h>
#include <symengine/basic.h>
#include <symengine/symbol.h>
#include <symengine/add.h>
#include <symengine/mul.h>
#include <gsl/gsl_sf_hyperg.h>  header
#include <gsl/gsl_sf_gamma.h>  
#include <cmath>
#include <chrono>
#include <Eigen/Dense> 


using namespace std;
using namespace SymEngine;

// Function to compute the matrix
vector<vector<double>> DMP_(int N, int Ord, double beta, double c) {
	int cor = 1;
	vector<vector<double>> MN(N, vector<double>(Ord, 0.0));

	// Define symbolic variables
	RCP<const Symbol> nn = symbol("nn");
	RCP<const Symbol> xx = symbol("xx");
	RCP<const Symbol> bb = symbol("bb");
	RCP<const Symbol> cc = symbol("cc");

	// Define the symbolic expression
	RCP<const Basic> r = mul({ rf(bb, nn), hyper({neg(nn), neg(xx)}, {bb}, sub(integer(1), div(integer(1), cc))),
							  sqrt(div(mul({pow(cc, xx), gamma(add(bb, xx)), pow(cc, nn), pow(sub(integer(1), cc), bb)}),
									   mul({factorial(xx), gamma(bb), factorial(nn), rf(bb, nn)}))) });

	// Define key indices
	int n1 = N / 2;
	int x1 = N / 2;
	int n0 = N / 2 - 1;
	int x0 = N / 2 - 1;

	// Initialize key elements
	map_basic_basic subs = { {nn, integer(x1)}, {xx, integer(n1)}, {bb, real_double(beta)}, {cc, real_double(c)} };
	MN[n1 + cor - 1][x1 + cor - 1] = eval_double(*r->subs(subs));

	subs = { {nn, integer(n1)}, {xx, integer(x0)}, {bb, real_double(beta)}, {cc, real_double(c)} };
	MN[n1 + cor - 1][x0 + cor - 1] = eval_double(*r->subs(subs));

	MN[n0 + cor - 1][x1 + cor - 1] = MN[n1 + cor - 1][x0 + cor - 1];

	subs = { {nn, integer(n0)}, {xx, integer(x0)}, {bb, real_double(beta)}, {cc, real_double(c)} };
	MN[n0 + cor - 1][x0 + cor - 1] = eval_double(*r->subs(subs));

	// Vectorized computations for "Initials" section
	for (int x = x0; x <= x1; ++x) {
		for (int n = n1 + 1; n < N; ++n) {
			double A = c / (c - 1);
			double B = (x - x * c - n + 1 - n * c + c - beta * c) / (1 - c);
			double C = (n - 1) * (n - 2 + beta) / (1 - c);
			double D = sqrt(c / (n * (beta + n - 1)));
			double E = sqrt(c * c / (n * (n - 1) * (beta + n - 2) * (beta + n - 1)));
			MN[n + cor - 1][x + cor - 1] = (B * D) / A * MN[n - 1 + cor - 1][x + cor - 1] + (C * E) / A * MN[n - 2 + cor - 1][x + cor - 1];
			MN[x + cor - 1][n + cor - 1] = MN[n + cor - 1][x + cor - 1];
		}
	}

	// Optimize A1 using vectorized computations
	for (int n = n0 - 1; n < N; ++n) {
		for (int x = x1; x < n; ++x) {
			double c1 = (n * (c - 1) + x + (x + beta) * c) * sqrt(1 / (c * (x + 1) * (x + beta)));
			double c2 = -sqrt((x * (x + beta - 1)) / ((x + 1) * (x + beta)));
			MN[n + cor - 1][x + 1 + cor - 1] = c1 * MN[n + cor - 1][x + cor - 1] + c2 * MN[n + cor - 1][x - 1 + cor - 1];
			MN[x + 1 + cor - 1][n + cor - 1] = MN[n + cor - 1][x + 1 + cor - 1];
		}
	}

	// Optimize A2
	for (int n = n1 - 1; n < N; ++n) {
		int x = x0;
		while (x >= 1) {
			double c1 = (n * (c - 1) + x + (x + beta) * c) * sqrt(1 / (c * (x + 1) * (x + beta)));
			double c2 = -sqrt((x * (x + beta - 1)) / ((x + 1) * (x + beta)));
			MN[n + cor - 1][x - 1 + cor - 1] = -c1 / c2 * MN[n + cor - 1][x + cor - 1] + 1 / c2 * MN[n + cor - 1][x + 1 + cor - 1];
			MN[x - 1 + cor - 1][n + cor - 1] = MN[n + cor - 1][x - 1 + cor - 1];
			if (abs(MN[n + cor - 1][x - 1 + cor - 1]) < 1e-7 && abs(MN[n + cor - 1][x + cor - 1]) < 1e-5) {
				break;
			}
			x--;
		}
	}

	// Optimize A3
	for (int x = n1 - 2; x >= 0; --x) {
		int n = x + 2;
		while (n >= 2) {
			double A = c / (c - 1);
			double B = (x - x * c - n + 1 - n * c + c - beta * c) / (1 - c);
			double C = (n - 1) * (n - 2 + beta) / (1 - c);
			double D = sqrt(c / (n * (beta + n - 1)));
			double E = sqrt(c * c / (n * (n - 1) * (beta + n - 2) * (beta + n - 1)));
			MN[n - 2 + cor - 1][x + cor - 1] = A / (C * E) * MN[n + cor - 1][x + cor - 1] - (B * D) / (C * E) * MN[n - 1 + cor - 1][x + cor - 1];
			MN[x + cor - 1][n - 2 + cor - 1] = MN[n - 2 + cor - 1][x + cor - 1];
			if (abs(MN[n - 2 + cor - 1][x + cor - 1]) < 1e-6 && abs(MN[n - 1 + cor - 1][x + cor - 1]) < 1e-4) {
				break;
			}
			n--;
		}
	}

	return MN;
}



int main() {
	std::cout << " Starting" << std::endl;
	int NN = 1000; // Example size of the signal
	double beta = NN / 4;
	double c = 0.1;
	int Ord = NN; // Example order of polynomials

	// Start timing
	auto start = std::chrono::high_resolution_clock::now();

	auto CH = DMP_(NN, Ord, beta, c);

	// Stop timing
	auto end = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> elapsed = end - start;

	std::cout << "Computation time: " << elapsed.count() << " seconds" << std::endl;

	return 0;
}