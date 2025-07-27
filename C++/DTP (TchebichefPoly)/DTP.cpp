#include <iostream>
#include <vector>
#include <cmath>
#include <chrono>

using namespace std;

std::vector<std::vector<double>> DTP_(int N, int Ord) {
    // Initialize matrix with Ord rows and N cols
    std::vector<std::vector<double>> R(Ord, std::vector<double>(N, 0.0));

    R[0][0] = 1.0 / std::sqrt(double(N));

    for (int n = 1; n < Ord; ++n)
        R[n][0] = -std::sqrt((double(N) - n) / (double(N) + n)) *
                  std::sqrt((2.0 * n + 1) / (2.0 * n - 1)) * R[n - 1][0];

    for (int n = 0; n < Ord; ++n)
        R[n][1] = (1.0 + n * (n + 1.0) / (1.0 - double(N))) * R[n][0];

    int Ord1 = std::min(static_cast<int>(std::floor(N / 4.0)), Ord);

    for (int x = 2; x < N / 2; ++x) {
        for (int n = 0; n < Ord1; ++n) {
            double b1 = -n * (n + 1.0) - (2.0 * x - 1.0) * (x - N - 1.0) - x;
            double lamda1 = b1 / (x * (N - x));
            double lamda2 = ((x - 1.0) * (x - N - 1.0)) / (x * (N - x));
            R[n][x] = lamda1 * R[n][x - 1] + lamda2 * R[n][x - 2];
        }
    }

    for (int n = Ord1; n < Ord; ++n) {
        int xx = static_cast<int>(std::floor(0.5 * N - std::sqrt(std::pow(N * 0.5, 2) - std::pow(n / 2.0, 2)))) % Ord;
        for (int x = xx; x < N / 2; ++x) {
            double a1 = (2.0 / n) * std::sqrt((4 * std::pow(n, 2) - 1) / (std::pow(N, 2) - std::pow(n, 2)));
            double a2 = ((1 - N) / double(n)) * std::sqrt((4 * std::pow(n, 2) - 1) / (std::pow(N, 2) - std::pow(n, 2)));
            double a3 = ((1 - n) / double(n)) * std::sqrt((2.0 * n + 1) / (2.0 * n - 3)) *
                        std::sqrt((std::pow(N, 2) - std::pow(n - 1, 2)) / (std::pow(N, 2) - std::pow(n, 2)));
            R[n][x] = a1 * x * R[n - 1][x] + a2 * R[n - 1][x] + a3 * R[n - 2][x];
        }
    }

    for (int n = Ord - 1; n >= Ord1; --n) {
        int xx = static_cast<int>(std::floor(0.5 * N - std::sqrt(std::pow(N * 0.5, 2) - std::pow(n / 2.0, 2)))) % Ord;
        for (int x = xx + 1; x >= 2; --x) {
            double b1 = -n * (n + 1.0) - (2.0 * x - 1.0) * (x - N - 1.0) - x;
            double lamda1 = b1 / (x * (N - x));
            double lamda2 = ((x - 1.0) * (x - N - 1.0)) / (x * (N - x));
            R[n][x - 2] = (R[n][x] - lamda1 * R[n][x - 1]) / lamda2;
            if (std::abs(R[n][x - 2]) > std::abs(R[n][x - 1]))
                break;
        }
    }

    for (int x = N / 2; x < N; ++x)
        for (int n = 0; n < Ord; ++n)
            R[n][x] = R[n][N - 1 - x] / std::pow(-1.0, n);

    return R;
}


int main() {
    int Poly_size = 2000;
    int Ord = 2000;

	// Start timing
    auto start = std::chrono::high_resolution_clock::now();
    auto T = DTP_(Poly_size, Ord);
    auto end = std::chrono::high_resolution_clock::now();

	// End timing
    std::chrono::duration<double> duration = end - start;
    std::cout << "Computation time: " << duration.count() << " seconds\n";

	std::cout << "Matrix size: " << T.size() << " rows × "
	    << (T.empty() ? 0 : T[0].size()) << " columns" << std::endl;
		
    return 0;
}
