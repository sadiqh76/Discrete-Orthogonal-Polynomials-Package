import numpy as np
import time
from math import lgamma

def DRP_(N, a, alpha, beta):
    """
    DRP_: Computes the Discrete Racah Polynomials (DRP)

     Inputs:
       N   - Size of the polynomial matrix (maximum degree + 1)
       a   - Racah polynomial parameter (positive integer)
       ph  - Racah parameter phi (must satisfy ph <= a)
       bt  - Racah parameter beta (must satisfy bt <= ph)

     Output:
       R   - Matrix of size (N x N) containing the DRP values

     Notes:
       - The function uses recurrence relations and gamma functions
         to construct the matrix of Racah polynomials.
       - Numerical stability is managed via thresholds and dynamic index bounds.
       - Ensure the input parameters satisfy: 0 < bt < ph < a.
    """
    # Constraints
    assert alpha <= a, "ph must be less than or equal to a"
    assert beta <= alpha, "bt must be less than or equal to ph"

    N_max = N
    b = a + N
    thres = 1e-5

    # Initialize the result matrix R
    R = np.zeros((N, N))
    cor = 0  # Correction for 0-based indexing in Python

    # Compute initial value
    FF1 = (lgamma(alpha + beta + 2) + lgamma(2 * a + N) + lgamma(beta + N)
        + lgamma(2 * a + 2 * N + alpha)  - ( lgamma(2 * a + 2 * N - 1)
        + lgamma(beta + 1) + lgamma(alpha + beta + N + 1) + lgamma(2 * a + N + alpha + 1)))
    F1 = np.exp(FF1 / 2)
    R[0 + cor, N - 1 + cor] = F1

    for s in range(a + N - 2, a - 1, -1):
        R[0 + cor, s - a + cor] = np.sqrt((2 * s + 1) * (a - beta + s + 1) * (b + s + 1) * (b + alpha - s - 1)
           * (a - s - 1) / ((a + s + 1) * (b + alpha + s + 1) * (a - beta - s - 1) * (2 * s + 3) * (b - s - 1)))\
            * R[0 + cor, s - a + 1 + cor]

    for s in range(a, a + N):
        R[1 + cor, s - a + cor] = (
            -(
                    ((-a + b - 1) * alpha + b ** 2 - s ** 2 - a - s - 1) * beta
                    + (a**2 - s**2 + b - s - 1) * alpha + a ** 2 + b ** 2 - 2 * (s**2 + s) - 1)
            * np.sqrt((alpha + beta + 3) / ((a - b + 1) * (a + b - beta - 1) * (alpha + 1) * (beta + 1)
                    * (a - b - alpha - beta - 1) * (a + b + alpha + 1))) * R[0 + cor, s - a + cor]
        )

    s = N - 1
    for n in range(1, N_max - 1):
        R[n + 1 + cor, s + cor] = np.sqrt(  (N - n - 1) * (alpha + beta + 2 * n + 3)
            * (alpha + beta + n + 1) * (alpha + n + 1) * (2 * a + N - beta - n - 1)
            / ((alpha + beta + 2 * n + 1) * (beta + n + 1) * (N + alpha + beta + n + 1)
            * (n + 1) * (2 * a + N + alpha + n + 1))) * R[n + cor, s + cor]

        R[n + 1 + cor, cor] = (
            -np.sqrt((N - n - 1) * (alpha + beta + 2 * n + 3) * ((alpha + beta + n + 1) * (beta + n + 1))
                * (2 * a + N + alpha + n + 1) / ( (2 * a + N - beta - n - 1) * (alpha + beta + 2 * n + 1)
                    * (alpha + n + 1) * (N + alpha + beta + n + 1) * (n + 1)))* R[n + cor, cor]
        )

    # Find max indices
    maxind_N_a_1 = np.argmax(R[:, -1]) + 1
    maxind_0 = np.argmax(R[:, 0]) + 1

    if maxind_N_a_1 == 1 and N > 4:
        print(f"Unable to compute initial values for a={a}, b={b}, ph={alpha}, bt={beta}")
        print("Please terminate")

    if maxind_0 == 1:
        maxind_0 = maxind_N_a_1

    maxind_N_a_1 = min(maxind_N_a_1, N_max - 1)
    maxind_0 = min(maxind_0, N_max - 1)
    if maxind_N_a_1 <= 1:
        maxind_N_a_1 = 2

    # Part 1
    for s in range(a + 1, a + N - 1):
        mis = round(maxind_0 + (s - a) * (maxind_N_a_1 - maxind_0) / (N - 1))
        for n in range(2, int(mis) + 1):
            A = (n * (alpha + beta + n)) / ((alpha + beta + 2 * n - 1) * (alpha + beta + 2 * n))
            B = (s * (s + 1)
                - (1 / 4) * (a ** 2 + b ** 2 + (a - beta) ** 2 + (b + alpha) ** 2 - 2)
                + (1 / 8) * ((alpha + beta + 2 * n - 2) * (alpha + beta + 2 * n))
                - (
                    ((beta ** 2 - alpha ** 2) * ((b + alpha / 2) ** 2 - (a - beta / 2) ** 2))
                    / (2 * (alpha + beta + 2 * n - 2) * (alpha + beta + 2 * n))
                )
            )
            C = (-(alpha + n - 1)
                * (beta + n - 1)
                * ((a + b + (alpha - beta) / 2) ** 2 - (n - 1 + (alpha + beta) / 2) ** 2)
                * ((b - a + (alpha + beta) / 2) ** 2 - (n - 1 + (alpha + beta) / 2) ** 2)
                / ((alpha + beta + 2 * n - 2) * (alpha + beta + 2 * n - 1))
            )
            B1 = (alpha + beta + 2 * n + 1) * n * (alpha + beta + n) / (
                (alpha + beta + 2 * n - 1)
                * (-b + a + n)
                * (a + b - beta - n)
                * (alpha + n)
                * (beta + n)
                * (-b + a - alpha - beta - n)
                * (a + b + alpha + n)
            )
            C1 = ((n - 1) * n * (alpha + beta + n - 1) * (alpha + beta + n) * (alpha + beta + 2 * n + 1)
                / ((alpha + n - 1) * (alpha + n) * (beta + n - 1) * (beta + n) * (-b + a - alpha - beta - n + 1)
                    * (-b + a - alpha - beta - n) * (a + b + alpha + n - 1) * (a + b + alpha + n)
                    * (alpha + beta + 2 * n - 3) * (-b + a + n) * (-b + a + n - 1) * (a + b - beta - n)
                    * (a + b - beta - n + 1))
            )
            P1 = (B / A) * np.sqrt(B1)
            P2 = (C / A) * np.sqrt(C1)
            R[n + cor, s - a + cor] = P1 * R[n - 1 + cor, s - a + cor] + P2 * R[n - 2 + cor, s - a + cor]

    # Part 2
    for s in range(a + 1, a + N - 1):
        mis = round(maxind_0 + (s - a) * (maxind_N_a_1 - maxind_0) / (N - 1))
        for n in range(int(mis) + 1, N_max):
            A = (n * (alpha + beta + n)) / ((alpha + beta + 2 * n - 1) * (alpha + beta + 2 * n))
            B = (
                s * (s + 1)
                - (1 / 4) * (a ** 2 + b ** 2 + (a - beta) ** 2 + (b + alpha) ** 2 - 2)
                + (1 / 8) * ((alpha + beta + 2 * n - 2) * (alpha + beta + 2 * n))
                - (
                    ((beta ** 2 - alpha ** 2) * ((b + alpha / 2) ** 2 - (a - beta / 2) ** 2))
                    / (2 * (alpha + beta + 2 * n - 2) * (alpha + beta + 2 * n))
                )
            )
            C = (-(alpha + n - 1) * (beta + n - 1)
                * ((a + b + (alpha - beta) / 2) ** 2 - (n - 1 + (alpha + beta) / 2) ** 2)
                * ((b - a + (alpha + beta) / 2) ** 2 - (n - 1 + (alpha + beta) / 2) ** 2)
                / ((alpha + beta + 2 * n - 2) * (alpha + beta + 2 * n - 1))
            )
            B1 = (alpha + beta + 2 * n + 1) * n * (alpha + beta + n) / ( (alpha + beta + 2 * n - 1)
                * (-b + a + n) * (a + b - beta - n) * (alpha + n) * (beta + n) * (-b + a - alpha - beta - n)
                * (a + b + alpha + n))
            C1 = ((n - 1) * n * (alpha + beta + n - 1) * (alpha + beta + n)  * (alpha + beta + 2 * n + 1)
                / ( (alpha + n - 1) * (alpha + n) * (beta + n - 1) * (beta + n) * (-b + a - alpha - beta - n + 1)
                    * (-b + a - alpha - beta - n) * (a + b + alpha + n - 1) * (a + b + alpha + n)
                    * (alpha + beta + 2 * n - 3) * (-b + a + n) * (-b + a + n - 1) * (a + b - beta - n)
                    * (a + b - beta - n + 1))
            )
            P1 = (B / A) * np.sqrt(B1)
            P2 = (C / A) * np.sqrt(C1)
            R[n + cor, s - a + cor] = P1 * R[n - 1 + cor, s - a + cor] + P2 * R[n - 2 + cor, s - a + cor]
            if abs(R[n + cor, s - a + cor ]) <= thres and (
                abs(R[n + cor, s - a + cor]) > abs(R[n - 1 + cor, s - a + cor])
            ):
                R[n + cor, s - a + cor] = 0
                break

    return R


# Example usage:
# Parameter setting
N = 1600
a=int(N/2)
alpha=int(N/4)
beta=int(N/8)

# Call to generate Racah Polynomials Matrix
R_matrix = DRP_(N,a,alpha,beta)
