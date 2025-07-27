import numpy as np
from math import lgamma
import time
import pandas as pd
from sympy import symbols, hyper, gamma, factorial, sqrt, rf

def DMP_(N, Ord, beta, c):
    """
    DMP_: Computes the Discrete Meixner Polynomials (DMP)

     Inputs:
       N     - Size of the signal (number of samples)
       Ord   - Maximum order of the polynomials (1 <= Ord <= N)
       c     - The Meixner parameter c (0 < c < 1)
       beta  - The Meixner parameter beta

     Output:
       R - Matrix of size (N x Ord) containing the DMP values
    """

    # Define constants and preallocate
    cor = 0 # Python indexing for matrices starts from 0 (zero-based indexing)
    MN = np.zeros((N, N))

    # Define symbolic function
    nn, xx, bb, cc = symbols('nn xx bb cc')
    r = rf(bb, nn) * hyper([-nn, -xx], [bb], 1 - 1/cc) * \
        sqrt(cc**xx * gamma(bb + xx) * cc**nn * (1 - cc)**bb / \
             (factorial(xx) * gamma(bb) * factorial(nn) * rf(bb, nn)))

    # Define key indices
    n1 = N // 2
    x1 = N // 2
    n0 = N // 2 - 1
    x0 = N // 2 - 1

    # Initialize key elements
    MN[n1 + cor, x1 + cor] = r.subs({nn: x1, xx: n1, bb: beta, cc: c}).evalf()
    MN[n1 + cor, x0 + cor] = r.subs({nn: n1, xx: x0, bb: beta, cc: c}).evalf()
    MN[n0 + cor, x1 + cor] = MN[n1 + cor, x0 + cor]
    MN[n0 + cor, x0 + cor] = r.subs({nn: n0, xx: x0, bb: beta, cc: c}).evalf()

    # Vectorized computations for "Initials" section
    for x in range(x0, x1 + 1):
        for n in range(n1 + 1, N - 0):
            A = c / (c - 1)
            B = (x - x * c - n + 1 - n * c + c - beta * c) / (1 - c)
            C = (n - 1) * (n - 2 + beta) / (1 - c)
            D = np.sqrt(c / (n * (beta + n - 1)))
            E = np.sqrt(c**2 / (n * (n - 1) * (beta + n - 2) * (beta + n - 1)))
            MN[n + cor, x + cor] = (B * D) / A * MN[n - 1 + cor, x + cor] + (C * E) / A * MN[n - 2 + cor, x + cor]
            MN[x + cor, n + cor] = MN[n + cor, x + cor]

    # Optimize A1 using vectorized computations
    for n in range(n0 - 1, N - 0):
        for x in range(x1, n - 0):
            c1 = (n * (c - 1) + x + (x + beta) * c) * np.sqrt(1 / (c * (x + 1) * (x + beta)))
            c2 = -np.sqrt((x * (x + beta - 1)) / ((x + 1) * (x + beta)))
            MN[n + cor, x + 1 + cor] = c1 * MN[n + cor, x + cor] + c2 * MN[n + cor, x - 1 + cor]
            MN[x + 1 + cor, n + cor] = MN[n + cor, x + 1 + cor]

    # Optimize A2
    for n in range(n1 - 1, N - 0):
        x = x0
        while x >= 1:
            c1 = (n * (c - 1) + x + (x + beta) * c) * np.sqrt(1 / (c * (x + 1) * (x + beta)))
            c2 = -np.sqrt((x * (x + beta - 1)) / ((x + 1) * (x + beta)))
            MN[n + cor, x - 1 + cor] = -c1 / c2 * MN[n + cor, x + cor] + 1 / c2 * MN[n + cor, x + 1 + cor]
            MN[x - 1 + cor, n + cor] = MN[n + cor, x - 1 + cor]
            if abs(MN[n + cor, x - 1 + cor]) < 1e-7 and abs(MN[n + cor, x + cor]) < 1e-5:
                break
            x -= 1

    # Optimize A3
    for x in range(n1 - 2, -1, -1):
        n = x + 2
        while n >= 2:
            A = c / (c - 1)
            B = (x - x * c - n + 1 - n * c + c - beta * c) / (1 - c)
            C = (n - 1) * (n - 2 + beta) / (1 - c)
            D = np.sqrt(c / (n * (beta + n - 1)))
            E = np.sqrt(c**2 / (n * (n - 1) * (beta + n - 2) * (beta + n - 1)))
            MN[n - 2 + cor, x + cor] = A / (C * E) * MN[n + cor, x + cor] - (B * D) / (C * E) * MN[n - 1 + cor, x + cor]
            MN[x + cor, n - 2 + cor] = MN[n - 2 + cor, x + cor]
            if abs(MN[n - 2 + cor, x + cor]) < 1e-6 and abs(MN[n - 1 + cor, x + cor]) < 1e-4:
                break
            n -= 1
    R=MN[:Ord, :N]
    return R


# Example usage:
# Parameter setting
N = 1000
beta = N/4
c=0.3
Order = N

# Call to generate Meixner Polynomials Matrix
R_matrix = DMP_(N, Order, beta, c)

