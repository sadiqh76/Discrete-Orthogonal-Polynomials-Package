import numpy as np
import time
import math
import matplotlib.pyplot as plt

from scipy.special import gammaln


def DHP_(N, alpha, beta, Ord):
    """
    DHP_: Computes the Discrete Hahn Polynomials (DHP)

    Inputs:
       N      - Size of the signal (must be an even number)
       Ord    - Maximum order of the polynomials (1 <= Ord <= N)
       alpha  - Hahn polynomial parameter ?
       beta   - Hahn polynomial parameter ?

    Output:
       R - Matrix of size (Ord x N) containing the Hahn polynomial values
           computed using a recurrence relation in both n and x directions
    """
    H = np.zeros((N, N))
    M1 = int(np.floor(0.4 * N))
    M = min(M1, Ord)  # To compute HPL to the required order if Ord < M

    # Compute initial value h(0,0)
    H[0, 0] = np.exp((gammaln(alpha + beta + 2) + gammaln(N + alpha) - (gammaln(alpha + 1) + gammaln(N + alpha + beta + 1))) / 2)

    # Compute value of HLP h(n,0); n=1,2,3,...
    for n in range(M - 1):
        H[n + 1, 0] = -np.sqrt((2 * n + 3 + alpha + beta) * (n + beta + 1) * (n + alpha + beta + 1) * (-n - 1 + N) /
                               ((n + 1 + alpha + beta + N) * (2 * n + alpha + beta + 1) * (n + 1 + alpha) * (n + 1))) * H[n, 0]

    # Compute HPL h(n,1); n=0,1,2,...
    for n in range(M):
        H[n, 1] = (((N - n - 1) * beta - n ** 2 - (alpha + 1) * n + N - 1) *
                   np.sqrt((beta + 1) * (N - 1) / (alpha + N - 1)) / ((beta + 1) * (N - 1))) * H[n, 0]

    # Compute DHPCs values in part P1 using x-direction recurrence algorithm for x=2,3,... and n=0,1,..,M
    for x in range(2, int(N / 2)):
        B = np.sqrt((N - x) * (beta + x) * (N + alpha - x) * x)
        AA = (-2 * x ** 2 + (2 * N + alpha - beta + 2) * x + (beta - 1) * N - alpha - 1) / B
        CC = -np.sqrt((beta + x - 1) * (N - x + 1) * (x - 1) * (N + alpha - x + 1)) / B
        for n in range(M):
            BB = -n * (alpha + beta + n + 1) / B
            H[n, x] = (AA + BB) * H[n, x - 1] + CC * H[n, x - 2]

    # Compute DHPCs values in part P3 using both n and x direction in the range n=M, M+1,...,Ord and x=N/2-1,N/2-2, ..., adaptive location.
    for n in range(M, Ord):
        DA = np.sqrt((alpha + beta + 2 * n + 1) * (alpha + beta + 2 * n - 1) * (alpha + beta + 2 * n) ** 2 /
                     (n * (alpha + beta + n) * (N - n) * (alpha + n) * (beta + n) * (alpha + beta + N + n)))
        CEA = -np.sqrt((n - 1) * (alpha + n - 1) * (beta + n - 1) * (N - n + 1) * (alpha + beta + n - 1) *
                       (alpha + beta + 2 * n + 1) * (alpha + beta + N + n - 1) * (alpha + beta + 2 * n) ** 2 /
                       (n * (alpha + beta + n) * (alpha + n) * (beta + n) * (N - n) * (alpha + beta + 2 * n - 3) *
                        (alpha + beta + N + n) * (alpha + beta + 2 * n - 2) ** 2))
        for x in range(int(N / 2) - 2, int(N / 2)):
            B = x - ((alpha - beta + 2 * N - 2) / 4) - ((beta ** 2 - alpha ** 2) * (alpha + beta + 2 * N)) / (4 * (alpha + beta + 2 * n - 2) * (alpha + beta + 2 * n))
            H[n, x] = (B * DA) * H[n - 1, x] + CEA * H[n - 2, x]
        f = 0
        for x in range(int(N / 2) - 1, 1, -1):
            B = np.sqrt((N - x) * (beta + x) * (N + alpha - x) * x)
            AA = (-2 * x ** 2 + (2 * N + alpha - beta + 2) * x + (beta - 1) * N - alpha - 1) / B
            CC = -np.sqrt((beta + x - 1) * (N - x + 1) * (x - 1) * (N + alpha - x + 1)) / B
            BB = -n * (alpha + beta + n + 1) / B
            H[n, x - 2] = (H[n, x + 0] - (AA + BB) * H[n, x - 1]) / CC
            if abs(H[n, x - 2]) > abs(H[n, x - 1]) and abs(H[n, x - 2]) < 1e-6:
                if f == 0:
                    f = 1
                else:
                    H[n, x - 2] = 0
                    break

    # Compute h(0, N-1)
    H[0, N - 1] = np.exp((gammaln(beta + N) + gammaln(1 + alpha) - (gammaln(1 + beta) + gammaln(alpha + N))) / 2) * H[0, 0]

    # Compute value of HLP h(n, N-1); n=1,2,3,...
    for n in range(M - 1):
        H[n + 1, N - 1] = np.sqrt((2 * n + 3 + alpha + beta) * (n + alpha + beta + 1) * (n + 1 + alpha) * (-n - 1 + N) /
                                  ((n + beta + 1) * (n + 1 + alpha + beta + N) * (2 * n + alpha + beta + 1) * (n + 1))) * H[n, N - 1]

    # Compute DHPCs values h(n, N-2); n=0,1,2,...
    for n in range(M):
        H[n, N - 2] = ((-n ** 2 + (-alpha - beta - 1) * n + (N - 1) * (alpha + 1)) /
                       np.sqrt((beta + N - 1) * (N - 1) * (alpha + 1))) * H[n, N - 1]

    # Compute DHPCs values in part P2 using x-direction recurrence algorithm for x=N-1, N-2, ..., N/2
    for x in range(N - 1, int(N / 2) - 1, -1):
        B = np.sqrt((N - x) * (beta + x) * (N + alpha - x) * x)
        AA = (-2 * x ** 2 + (2 * N + alpha - beta + 2) * x + (beta - 1) * N - alpha - 1) / B
        CC = -np.sqrt((beta + x - 1) * (N - x + 1) * (x - 1) * (N + alpha - x + 1)) / B
        for n in range(M):
            BB = -n * (alpha + beta + n + 1) / B
            H[n, x - 2] = (H[n, x + 0] - (AA + BB) * H[n, x - 1]) / CC

    # Compute DHPCs values in part P4 using both n and x direction in the range n=M, M+1,...,Ord and x=N/2, N/2+1, ..., adaptive location.
    for n in range(M, Ord):
        DA = np.sqrt((alpha + beta + 2 * n + 1) * (alpha + beta + 2 * n - 1) * (alpha + beta + 2 * n) ** 2 /
                     (n * (alpha + beta + n) * (N - n) * (alpha + n) * (beta + n) * (alpha + beta + N + n)))
        CEA = -np.sqrt((n - 1) * (alpha + n - 1) * (beta + n - 1) * (N - n + 1) * (alpha + beta + n - 1) *
                       (alpha + beta + 2 * n + 1) * (alpha + beta + N + n - 1) * (alpha + beta + 2 * n) ** 2 /
                       (n * (alpha + beta + n) * (alpha + n) * (beta + n) * (N - n) * (alpha + beta + 2 * n - 3) *
                        (alpha + beta + N + n) * (alpha + beta + 2 * n - 2) ** 2))
        for x in range(int(N / 2), int(N / 2) + 2):
            B = x - ((alpha - beta + 2 * N - 2) / 4) - ((beta ** 2 - alpha ** 2) * (alpha + beta + 2 * N)) / (4 * (alpha + beta + 2 * n - 2) * (alpha + beta + 2 * n))
            H[n, x] = (B * DA) * H[n - 1, x] + CEA * H[n - 2, x]
        f = 0
        for x in range(int(N / 2), N):
            B = np.sqrt((N - x) * (beta + x) * (N + alpha - x) * x)
            AA = (-2 * x ** 2 + (2 * N + alpha - beta + 2) * x + (beta - 1) * N - alpha - 1) / B
            CC = -np.sqrt((beta + x - 1) * (N - x + 1) * (x - 1) * (N + alpha - x + 1)) / B
            BB = -n * (alpha + beta + n + 1) / B
            H[n, x] = (AA + BB) * H[n, x - 1] + CC * H[n, x - 2]
            if abs(H[n, x]) > abs(H[n, x - 1]) and abs(H[n, x]) < 1e-6:
                if f == 0:
                    f = 1
                else:
                    H[n, x - 1] = 0
                    break
    R = H[:Ord, :N]

    return R

# Example usage:
# Parameter setting
N = 500
Parameter_alpha = 100
Parameter_beta = 50
Ord = N-100

# Call to generate Hahn Polynomials Matrix
R_matrix = DHP_(N, Parameter_alpha, Parameter_beta, Ord)
