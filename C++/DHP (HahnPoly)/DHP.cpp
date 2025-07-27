#include <iostream>
#include <chrono>
#include <vector>
#include <cmath>

using namespace std;

// Function to compute the Discrete Hahn Polynomials
std::vector<std::vector<double>> DHP_(int N, double alpha, double beta, int Order) {
    int rows = std::min(N, N);
    std::vector<std::vector<double>> R(rows, std::vector<double>(N, 0.0));

    int M1 = static_cast<int>(floor(0.4 * N));
    int M = std::min(M1, rows);

    double Nd = (double)N;

    R[0][0] = exp((lgamma(alpha + beta + 2) + lgamma(Nd + alpha) -
                  (lgamma(alpha + 1) + lgamma(Nd + alpha + beta + 1))) / 2);

    for (int n = 0; n < M - 1; ++n) {
        double nd = (double)n;
        R[n + 1][0] = -sqrt((2 * nd + 3 + alpha + beta) * (nd + beta + 1) * (nd + alpha + beta + 1) * (-nd - 1 + Nd) /
                           ((nd + 1 + alpha + beta + Nd) * (2 * nd + alpha + beta + 1) * (nd + 1 + alpha) * (nd + 1))) * R[n][0];
    }

    for (int n = 0; n < M; ++n) {
		double nd = (double)n;
        R[n][1] = (((Nd - nd - 1) * beta - nd * nd - (alpha + 1) * nd + Nd - 1) * sqrt((beta + 1) * (Nd - 1) / (alpha + Nd - 1)) / ((beta + 1) * (Nd - 1))) * R[n][0];
    }

    for (int x = 2; x < N / 2; ++x) {
        double xd = (double)x;
        double B = sqrt((Nd - xd) * (beta + xd) * (Nd + alpha - xd) * xd);
        double AA = (-2 * xd * xd + (2 * Nd + alpha - beta + 2) * xd + (beta - 1) * Nd - alpha - 1) / B;
        double CC = -sqrt((beta + xd - 1) * (Nd - xd + 1) * (xd - 1) * (Nd + alpha - xd + 1)) / B;
        for (int n = 0; n < M; ++n) {
            double nd = (double)n;
            double BB = -nd * (alpha + beta + nd + 1) / B;
            R[n][x] = (AA + BB) * R[n][x - 1] + CC * R[n][x - 2];
        }
    }

    for (int n = M; n < rows; ++n) {
        double nd = (double)n;
        double DA = sqrt((alpha + beta + 2 * nd + 1) * (alpha + beta + 2 * nd - 1) * pow(alpha + beta + 2 * nd, 2) /
                         (nd * (alpha + beta + nd) * (Nd - nd) * (alpha + nd) * (beta + nd) * (alpha + beta + Nd + nd)));
        double CEA = -sqrt((nd - 1) * (alpha + nd - 1) * (beta + nd - 1) * (Nd - nd + 1) * (alpha + beta + nd - 1) *
                           (alpha + beta + 2 * nd + 1) * (alpha + beta + Nd + nd - 1) * pow(alpha + beta + 2 * nd, 2) /
                           (nd * (alpha + beta + nd) * (alpha + nd) * (beta + nd) * (Nd - nd) * (alpha + beta + 2 * nd - 3) *
                            (alpha + beta + Nd + nd) * pow(alpha + beta + 2 * nd - 2, 2)));

        for (int x = N / 2 - 2; x < N / 2; ++x) {
            double xd = (double)x;
            double B = xd - ((alpha - beta + 2 * Nd - 2) / 4) -
                       ((beta * beta - alpha * alpha) * (alpha + beta + 2 * Nd)) /
                       (4 * (alpha + beta + 2 * nd - 2) * (alpha + beta + 2 * nd));
            R[n][x] = (B * DA) * R[n - 1][x] + CEA * R[n - 2][x];
        }

        int f = 0;
        for (int x = N / 2 - 1; x > 1; --x) {
            double xd = (double)x;
            double B = sqrt((Nd - xd) * (beta + xd) * (Nd + alpha - xd) * xd);
            double AA = (-2 * xd * xd + (2 * Nd + alpha - beta + 2) * xd + (beta - 1) * Nd - alpha - 1) / B;
            double CC = -sqrt((beta + xd - 1) * (Nd - xd + 1) * (xd - 1) * (Nd + alpha - xd + 1)) / B;
            double BB = -nd * (alpha + beta + nd + 1) / B;
            R[n][x - 2] = (R[n][x] - (AA + BB) * R[n][x - 1]) / CC;
            if (abs(R[n][x - 2]) > abs(R[n][x - 1]) && abs(R[n][x - 2]) < 1e-6) {
                if (f == 0) f = 1;
                else {
                    R[n][x - 2] = 0;
                    break;
                }
            }
        }
    }

    R[0][N - 1] = exp((lgamma(beta + N) + lgamma(1 + alpha) - (lgamma(1 + beta) + lgamma(alpha + N))) / 2) * R[0][0];

    for (int n = 0; n < M - 1; ++n) {
        double nd = (double)n;
        R[n + 1][N - 1] = sqrt((2 * nd + 3 + alpha + beta) * (nd + alpha + beta + 1) * (nd + 1 + alpha) * (-nd - 1 + Nd) /
                              ((nd + beta + 1) * (nd + 1 + alpha + beta + Nd) * (2 * nd + alpha + beta + 1) * (nd + 1))) * R[n][N - 1];
    }

    for (int n = 0; n < M; ++n) {
        double nd = (double)n;
        R[n][N - 2] = ((-nd * nd + (-alpha - beta - 1) * nd + (Nd - 1) * (alpha + 1)) /
                      sqrt((beta + N - 1) * (Nd - 1) * (alpha + 1))) * R[n][N - 1];
    }

    for (int x = N - 1; x >= N / 2; --x) {
        double xd = (double)x;
        double B = sqrt((Nd - xd) * (beta + xd) * (Nd + alpha - xd) * xd);
        double AA = (-2 * xd * xd + (2 * Nd + alpha - beta + 2) * xd + (beta - 1) * Nd - alpha - 1) / B;
        double CC = -sqrt((beta + xd - 1) * (Nd - xd + 1) * (xd - 1) * (Nd + alpha - xd + 1)) / B;
        for (int n = 0; n < M; ++n) {
            double nd = (double)n;
            double BB = -nd * (alpha + beta + nd + 1) / B;
            R[n][x - 2] = (R[n][x] - (AA + BB) * R[n][x - 1]) / CC;
        }
    }

    for (int n = M; n < rows; ++n) {
        double nd = (double)n;
        double DA = sqrt((alpha + beta + 2 * nd + 1) * (alpha + beta + 2 * nd - 1) * pow(alpha + beta + 2 * nd, 2) /
                         (nd * (alpha + beta + nd) * (Nd - nd) * (alpha + nd) * (beta + nd) * (alpha + beta + Nd + nd)));
        double CEA = -sqrt((nd - 1) * (alpha + nd - 1) * (beta + nd - 1) * (Nd - nd + 1) * (alpha + beta + nd - 1) *
                           (alpha + beta + 2 * nd + 1) * (alpha + beta + Nd + nd - 1) * pow(alpha + beta + 2 * nd, 2) /
                           (nd * (alpha + beta + nd) * (alpha + nd) * (beta + nd) * (Nd - nd) * (alpha + beta + 2 * nd - 3) *
                            (alpha + beta + Nd + nd) * pow(alpha + beta + 2 * nd - 2, 2)));
        for (int x = N / 2; x < N / 2 + 2; ++x) {
            double xd = (double)x;
            double B = xd - ((alpha - beta + 2 * Nd - 2) / 4) -
                       ((beta * beta - alpha * alpha) * (alpha + beta + 2 * Nd)) /
                       (4 * (alpha + beta + 2 * nd - 2) * (alpha + beta + 2 * nd));
            R[n][x] = (B * DA) * R[n - 1][x] + CEA * R[n - 2][x];
        }

        int f = 0;
        for (int x = N / 2; x < N; ++x) {
            double xd = (double)x;
            double B = sqrt((Nd - xd) * (beta + xd) * (Nd + alpha - xd) * xd);
            double AA = (-2 * xd * xd + (2 * Nd + alpha - beta + 2) * xd + (beta - 1) * Nd - alpha - 1) / B;
            double CC = -sqrt((beta + xd - 1) * (Nd - xd + 1) * (xd - 1) * (Nd + alpha - xd + 1)) / B;
            double BB = -nd * (alpha + beta + nd + 1) / B;
            R[n][x] = (AA + BB) * R[n][x - 1] + CC * R[n][x - 2];
            if (abs(R[n][x]) > abs(R[n][x - 1]) && abs(R[n][x]) < 1e-6) {
                if (f == 0) f = 1;
                else {
                    R[n][x - 1] = 0;
                    break;
                }
            }
        }
    }
	// Trim to Order × N in case of Order less than N
	if (Order < N) {
		R.resize(Order);  // Trim rows to Order
		return R;
	}

	return R;
}


int main() {
	int Poly_Size = 1000;
	double Parameter_alpha = 100;
	double Parameter_beta = 50;
	int Ord = 1000;
	
	// Start timing
	auto start_time = chrono::high_resolution_clock::now();
	vector<vector<double>> H = DHP_(Poly_Size, Parameter_alpha, Parameter_beta, Ord);
	
	// End timing
	auto end_time = chrono::high_resolution_clock::now();

	chrono::duration<double> elapsed_time = end_time - start_time;
	cout << "Computation time: " << elapsed_time.count() << " seconds" << endl;

    std::cout << "Matrix size: " << H.size() << " rows × "
	    << (H.empty() ? 0 : H[0].size()) << " columns" << std::endl;

	return 0;
}