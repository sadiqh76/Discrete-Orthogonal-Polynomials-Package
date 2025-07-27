#include <iostream>
#include <vector>
#include <cmath>
#include <chrono>
#include <stdexcept>
#include <algorithm> 

using namespace std;

vector<vector<double>> DRP_(int N, double a, double alpha, double beta) {
    // Constraints
    if (alpha > a) {
        throw invalid_argument("alpha must be less than or equal to a");
    }
    if (beta > alpha) {
        throw invalid_argument("beta must be less than or equal to alpha");
    }

    int N_max = N;
    double b = a + N;
    double thres = 1e-5;

    vector<vector<double>> R(N, vector<double>(N, 0.0));
    int cor = 0;

    double FF1 = (
        lgamma(alpha + beta + 2)
        + lgamma(2 * a + N)
        + lgamma(beta + N)
        + lgamma(2 * a + 2 * N + alpha)
        - (
            lgamma(2 * a + 2 * N - 1)
            + lgamma(beta + 1)
            + lgamma(alpha + beta + N + 1)
            + lgamma(2 * a + N + alpha + 1)
        )
    );
    double F1 = exp(FF1 / 2);
    R[0 + cor][N - 1 + cor] = F1;

    for (int s = a + N - 2; s >= a; --s) {
        double numerator = (2 * s + 1) * (a - beta + s + 1) * (b + s + 1) * (b + alpha - s - 1) * (a - s - 1);
        double denominator = (a + s + 1) * (b + alpha + s + 1) * (a - beta - s - 1) * (2 * s + 3) * (b - s - 1);
        R[0 + cor][s - a + cor] = sqrt(numerator / denominator) * R[0 + cor][s - a + 1 + cor];
    }

    for (int s = a; s < a + N; ++s) {
        double numerator = -(
            ((-a + b - 1) * alpha + b * b - s * s - a - s - 1) * beta
            + (a * a - s * s + b - s - 1) * alpha
            + a * a + b * b - 2 * (s * s + s) - 1
        );
        double denominator = (a - b + 1) * (a + b - beta - 1) * (alpha + 1) * (beta + 1) * (a - b - alpha - beta - 1) * (a + b + alpha + 1);
        R[1 + cor][s - a + cor] = numerator * sqrt((alpha + beta + 3) / denominator) * R[0 + cor][s - a + cor];
    }

    int s = N - 1;
    for (int n = 1; n < N_max - 1; ++n) {
        double numerator = (N - n - 1) * (alpha + beta + 2 * n + 3) * (alpha + beta + n + 1) * (alpha + n + 1) * (2 * a + N - beta - n - 1);
        double denominator = (alpha + beta + 2 * n + 1) * (beta + n + 1) * (N + alpha + beta + n + 1) * (n + 1) * (2 * a + N + alpha + n + 1);
        R[n + 1 + cor][s + cor] = sqrt(numerator / denominator) * R[n + cor][s + cor];

        numerator = (N - n - 1) * (alpha + beta + 2 * n + 3) * ((alpha + beta + n + 1) * (beta + n + 1)) * (2 * a + N + alpha + n + 1);
        denominator = (2 * a + N - beta - n - 1) * (alpha + beta + 2 * n + 1) * (alpha + n + 1) * (N + alpha + beta + n + 1) * (n + 1);
        R[n + 1 + cor][cor] = -sqrt(numerator / denominator) * R[n + cor][cor];
    }

    auto max_N_a_1 = max_element(R.begin(), R.end(), [](const vector<double>& a, const vector<double>& b) {
        return a.back() < b.back();
    });
    int maxind_N_a_1 = distance(R.begin(), max_N_a_1) + 1;

    auto max_0 = max_element(R.begin(), R.end(), [](const vector<double>& a, const vector<double>& b) {
        return a[0] < b[0];
    });
    int maxind_0 = distance(R.begin(), max_0) + 1;

    if (maxind_N_a_1 == 1 && N > 4) {
        cout << "Unable to compute initial values for a=" << a << ", b=" << b << ", alpha=" << alpha << ", beta=" << beta << endl;
        cout << "Please terminate" << endl;
    }

    if (maxind_0 == 1) {
        maxind_0 = maxind_N_a_1;
    }

    maxind_N_a_1 = min(maxind_N_a_1, N_max - 1);
    maxind_0 = min(maxind_0, N_max - 1);
    if (maxind_N_a_1 <= 1) maxind_N_a_1 = 2;

    for (int s = a + 1; s < a + N - 1; ++s) {
        int mis = static_cast<int>(round(maxind_0 + (s - a) * (maxind_N_a_1 - maxind_0) / (N - 1)));
        for (int n = 2; n <= mis; ++n) {
            double A = (n * (alpha + beta + n)) / ((alpha + beta + 2 * n - 1) * (alpha + beta + 2 * n));
            double B = s * (s + 1)
                - 0.25 * (a * a + b * b + (a - beta) * (a - beta) + (b + alpha) * (b + alpha) - 2)
                + 0.125 * ((alpha + beta + 2 * n - 2) * (alpha + beta + 2 * n))
                - (((beta * beta - alpha * alpha) * ((b + alpha / 2) * (b + alpha / 2) - (a - beta / 2) * (a - beta / 2))) / (2 * (alpha + beta + 2 * n - 2) * (alpha + beta + 2 * n)));
            double C = -(alpha + n - 1) * (beta + n - 1)
                * ((a + b + (alpha - beta) / 2) * (a + b + (alpha - beta) / 2) - (n - 1 + (alpha + beta) / 2) * (n - 1 + (alpha + beta) / 2))
                * ((b - a + (alpha + beta) / 2) * (b - a + (alpha + beta) / 2) - (n - 1 + (alpha + beta) / 2) * (n - 1 + (alpha + beta) / 2))
                / ((alpha + beta + 2 * n - 2) * (alpha + beta + 2 * n - 1));
            double B1 = (alpha + beta + 2 * n + 1) * n * (alpha + beta + n)
                / ((alpha + beta + 2 * n - 1) * (-b + a + n) * (a + b - beta - n) * (alpha + n) * (beta + n) * (-b + a - alpha - beta - n) * (a + b + alpha + n));
            double C1 = (n - 1) * n * (alpha + beta + n - 1) * (alpha + beta + n) * (alpha + beta + 2 * n + 1)
                / ((alpha + n - 1) * (alpha + n) * (beta + n - 1) * (beta + n)
                    * (-b + a - alpha - beta - n + 1) * (-b + a - alpha - beta - n)
                    * (a + b + alpha + n - 1) * (a + b + alpha + n)
                    * (alpha + beta + 2 * n - 3)
                    * (-b + a + n) * (-b + a + n - 1)
                    * (a + b - beta - n) * (a + b - beta - n + 1));
            double P1 = (B / A) * sqrt(B1);
            double P2 = (C / A) * sqrt(C1);
            R[n + cor][s - a + cor] = P1 * R[n - 1 + cor][s - a + cor] + P2 * R[n - 2 + cor][s - a + cor];
        }
    }

    for (int s = a + 1; s < a + N - 1; ++s) {
        int mis = static_cast<int>(round(maxind_0 + (s - a) * (maxind_N_a_1 - maxind_0) / (N - 1)));
        for (int n = mis + 1; n < N_max; ++n) {
            double A = (n * (alpha + beta + n)) / ((alpha + beta + 2 * n - 1) * (alpha + beta + 2 * n));
            double B = s * (s + 1)
                - 0.25 * (a * a + b * b + (a - beta) * (a - beta) + (b + alpha) * (b + alpha) - 2)
                + 0.125 * ((alpha + beta + 2 * n - 2) * (alpha + beta + 2 * n))
                - (((beta * beta - alpha * alpha) * ((b + alpha / 2) * (b + alpha / 2) - (a - beta / 2) * (a - beta / 2))) / (2 * (alpha + beta + 2 * n - 2) * (alpha + beta + 2 * n)));
            double C = -(alpha + n - 1) * (beta + n - 1)
                * ((a + b + (alpha - beta) / 2) * (a + b + (alpha - beta) / 2) - (n - 1 + (alpha + beta) / 2) * (n - 1 + (alpha + beta) / 2))
                * ((b - a + (alpha + beta) / 2) * (b - a + (alpha + beta) / 2) - (n - 1 + (alpha + beta) / 2) * (n - 1 + (alpha + beta) / 2))
                / ((alpha + beta + 2 * n - 2) * (alpha + beta + 2 * n - 1));
            double B1 = (alpha + beta + 2 * n + 1) * n * (alpha + beta + n)
                / ((alpha + beta + 2 * n - 1) * (-b + a + n) * (a + b - beta - n) * (alpha + n) * (beta + n) * (-b + a - alpha - beta - n) * (a + b + alpha + n));
            double C1 = (n - 1) * n * (alpha + beta + n - 1) * (alpha + beta + n) * (alpha + beta + 2 * n + 1)
                / ((alpha + n - 1) * (alpha + n) * (beta + n - 1) * (beta + n)
                    * (-b + a - alpha - beta - n + 1) * (-b + a - alpha - beta - n)
                    * (a + b + alpha + n - 1) * (a + b + alpha + n)
                    * (alpha + beta + 2 * n - 3)
                    * (-b + a + n) * (-b + a + n - 1)
                    * (a + b - beta - n) * (a + b - beta - n + 1));
            double P1 = (B / A) * sqrt(B1);
            double P2 = (C / A) * sqrt(C1);
            R[n + cor][s - a + cor] = P1 * R[n - 1 + cor][s - a + cor] + P2 * R[n - 2 + cor][s - a + cor];
            if (abs(R[n + cor][s - a + cor]) <= thres &&
                abs(R[n + cor][s - a + cor]) > abs(R[n - 1 + cor][s - a + cor])) {
                R[n + cor][s - a + cor] = 0;
                break;
            }
        }
    }

    return R;
}



int main() {
	int N = 1000;
	double a = 200;
	double alpha = 100;
	double beta = 50;

	// Start timing
    auto start = chrono::high_resolution_clock::now();
	
	vector<vector<double>> R = DRP_(N, a, alpha, beta);

	// End timing
    auto end = chrono::high_resolution_clock::now();
    chrono::duration<double> duration = end - start;
    
	cout << "Computation time: " << duration.count() << " seconds" << endl;
	

	return 0;
}