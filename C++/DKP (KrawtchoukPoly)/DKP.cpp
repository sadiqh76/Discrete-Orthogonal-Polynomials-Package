#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>  
#include <chrono>

using namespace std;

std::vector<std::vector<double>> DKP_(int N, double Parameter_p, int Ord) {
    std::vector<std::vector<double>> R(N, std::vector<double>(N, 0.0));

    int p_signal = 0;
    if (Parameter_p > 0.5) {
        Parameter_p = 1 - Parameter_p;
        p_signal = 1;
    }

    R[0][0] = pow(1.0 - Parameter_p, (double(N) - 1.0) / 2.0);

    for (int x = 1; x < N; ++x) {
        R[0][x] = (sqrt(Parameter_p * (double(N) - double(x))) /
                  sqrt(double(x) * (1.0 - Parameter_p))) * R[0][x - 1];
    }

    for (int x = 1; x < N - 1; ++x) {
        R[1][x] = ((-double(x) + Parameter_p * (double(N) - 1.0)) / (Parameter_p * (double(N) - 1))) *
                   sqrt((double(N) - 1.0) * Parameter_p / (1.0 - Parameter_p)) * R[0][x];
    }

    for (int n = 2; n <= static_cast<int>(ceil(N / 2.0)); ++n) {
        for (int x = n; x < N - n; ++x) {
            double an = ((double(N) - 2 * double(n) + 1) * Parameter_p + double(n) - x - 1) /
                        sqrt(Parameter_p * double(n) * (1 - Parameter_p) * (double(N) - double(n)));
            double bn = sqrt((double(n) - 1) * (double(N) - double(n) + 1) /
                            (double(n) * (double(N) - double(n))));
            R[n][x] = an * R[n - 1][x] - bn * R[n - 2][x];
        }
    }

    for (int x = 0; x < static_cast<int>(floor(N / 2.0)); ++x) {
        for (int n = x + 1; n < N - x; ++n) {
            R[n][x] = R[x][n];
        }
    }

    for (int n = 1; n <= N; ++n) {
        for (int x = N - n + 2; x <= N; ++x) {
            if (n - 1 < N && x - 1 < N && N - n < N && N - x < N) {
                R[n - 1][x - 1] = pow(-1, N + n + x - 1) * R[N - n][N - x];
            }
        }
    }

    if (p_signal == 1) {
        for (auto& row : R) {
            std::reverse(row.begin(), row.end());
        }
        for (int i = 0; i < N; ++i) {
            double sign = pow(-1, i + 1);
            for (int j = 0; j < N; ++j) {
                R[i][j] *= sign;
            }
        }
    }

    // Trim to Ord × N in case Ord < N
    if (Ord < N) {
        R.resize(Ord);
    }

    return R;
}


int main() {
    int N = 2500;
    double p = 0.35;
    int Ord = 2500;
	
	// Start timing
    auto start = chrono::high_resolution_clock::now();

    vector<vector<double>> R = DKP_(N, p, Ord);
	
	// End timing
    auto end = chrono::high_resolution_clock::now();
    chrono::duration<double> duration = end - start;
    
	cout << "Computation time: " << duration.count() << " seconds" << endl;


    cout << "Matrix size: " << R.size() << " rows × "
         << (R.empty() ? 0 : R[0].size()) << " columns" << endl;

    return 0;
}
