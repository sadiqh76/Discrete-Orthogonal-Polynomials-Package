#include <iostream>
#include <vector>
#include <cmath>
#include <chrono> // For timing
#include <fstream> // For saving the generated matrix

std::vector<std::vector<double>> DCP_(int N, int Order, int Parameter_p) {
	// Preallocate NxN matrix
	std::vector<std::vector<double>> R(N, std::vector<double>(N, 0.0));

	// Initial values
	double aa = static_cast<double>(Parameter_p);
	Parameter_p = int(aa);
	int cor = 0;
	double init = std::sqrt(std::exp(-aa + (aa - 1) * std::log(aa) - std::lgamma(aa)));

	R[Parameter_p - 1 + cor][0 + cor] = init;
	R[Parameter_p + 0 + cor][0 + cor] = init;

	// Compute R(0, x) using R(0, 0)
	for (int n = Parameter_p - 1; n >= 1; --n) {
		R[n - 1 + cor][0 + cor] = R[n + cor][0 + cor] / std::sqrt(aa / double(n));
	}

	for (int n = Parameter_p + 1; n <= N - 1; ++n) {
		R[n + cor][0 + cor] = R[n - 1 + cor][0 + cor] * std::sqrt(aa / double(n));
	}

	// Compute R(1, x) using R(0, x)
	for (int n = 1; n <= N - 1; ++n) {
		R[n + cor][1 + cor] = R[n + cor][0 + cor] * (aa - double(n)) / std::sqrt(aa);
	}

	// Compute R(x, n) 
	for (int x = 1; x <= N - 2; ++x) {
		for (int n = x + 1; n <= N - 1; ++n) {
			double c1 = (aa - double(n) + double(x)) * std::sqrt(1.0 / (aa * (double(x) + 1)));
			double c2 = -std::sqrt(double(x) / (double(x) + 1));
			R[n + cor][x + 1 + cor] = c1 * R[n + cor][x + cor] + c2 * R[n + cor][x - 1 + cor];
		}
	}

	// Using symmetry relation
	for (int n = 1; n <= N - 1; ++n) {
		for (int x = 0; x <= n; ++x) {
			R[x + cor][n + cor] = R[n + cor][x + cor];
		}
	}

	// Trim to Order × N in case of Order less than N
	if (Order < N) {
		R.resize(Order);  // Trim rows to Order
		return R;
	}

	return R;
}


// Function to save a matrix to a CSV file
void saveMatrixToCSV(const std::vector<std::vector<double>>& matrix, const std::string& filename) {
	std::ofstream file(filename);
	if (file.is_open()) {
		for (const auto& row : matrix) {
			for (size_t j = 0; j < row.size(); ++j) {
				file << row[j];
				if (j < row.size() - 1) {
					file << ","; 
				}
			}
			file << "\n"; 
		}
		file.close();
		std::cout << "Matrix saved to " << filename << std::endl;
	}
	else {
		std::cerr << "Unable to open file " << filename << std::endl;
	}
}


int main() {
	std::cout << " Starting" << std::endl;
	int Poly_Size = 1500; // Example size of the signal
	int Ord = 1500; 		 // Example Orderer of polynomials
	double p = Poly_Size/6; // Example parameter

	// Start timing
	auto start = std::chrono::high_resolution_clock::now();

	auto CH = DCP_(Poly_Size , Ord, p);
	
	std::cout << CH.size()<< " x ";
	std::cout << CH[0].size()<< std::endl;;
	
	// Stop timing
	auto end = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> elapsed = end - start;

	std::cout << "Computation time: " << elapsed.count() << " seconds" << std::endl;
	// Saving the generated polynomial
	saveMatrixToCSV(CH,"a.csv");


	return 0;
}