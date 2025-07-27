import numpy as np
import time
# import numba

# @numba.jit(nopython=True, parallel=True)
def DKP_(N, p, Ord):
    """
    DKP_: Computes the Discrete Krawtchouk Polynomials (DKP)

     Inputs:
       N   - Size of the signal (maximum x and n indices)
       Ord - Maximum order of the polynomials (1 <= Ord <= N)
       p   - Probability parameter (0 < p < 1)

     Output:
       R - Matrix of size (Ord x N) containing the DKP values
    """
    K = np.zeros((N, N), dtype=float)

    # Adjust p to ensure it lies in [0, 0.5] for numerical stability
    p_signal = 0
    if p > 0.5:
        p = 1 - p
        p_signal = 1

    # Initialize the first row
    K[0, 0] = (1.0 - p) ** ((N - 1.0) / 2)
    for x in range(1, N):
        K[0, x] = (np.sqrt(p * (N - x)) / np.sqrt(x * (1.0 - p))) * K[0, x - 1]

    # Initialize the second row
    for x in range(1, N - 1):
        K[1, x] = ((-x + p * (N - 1.0)) / (p * (N - 1))) * np.sqrt((N - 1.0) * p / (1.0 - p)) * K[0, x]

    # Compute the remaining rows
    for n in range(2, int(np.ceil(N / 2))):
        for x in range(n, N - n):
            an = ((N - 2 * n + 1) * p + n - x - 1) / np.sqrt(p * n * (1 - p) * (N - n))
            bn = np.sqrt((n - 1) * (N - n + 1) / (n * (N - n)))
            K[n, x] = an * K[n - 1, x] - bn * K[n - 2, x]

    # Symmetrize the matrix
    for x in range(int(np.floor(N / 2))):
        for n in range(x + 1, N - x):
            K[n, x] = K[x, n]

    # Handle the lower triangular part
    for n in range(1, N + 1):
        for x in range(N - n + 1, N + 1):
            K[n - 1, x - 1] = (-1) ** (N + n + x - 1) * K[N - n, N - x]

    # Flip and adjust signs if p_signal is 1
    if p_signal == 1:
        K = np.fliplr(K)
        for i in range(N):
            K[i, :] *= (-1) ** (i + 1)

    # Extract the required portion
    R = K[:Ord, :N]
    return R

# Example usage:
# Parameter setting
N = 1000
p = 0.5
Order = N

# Call to generate Krawtchouk Polynomials Matrix
R_matrix = DKP_(N, p, Order)
