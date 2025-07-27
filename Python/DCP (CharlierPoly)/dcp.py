import numpy as np
from math import lgamma


def DCP_(N, Ord, Parameter_p):
    """
    DCP_: Computes the Discrete Charlier Polynomials (DCP)

    Inputs:
      N          - Size of the signal (maximum n index)
      Ord        - Maximum order of the polynomials (1 <= Ord <= N)
      Parameter_p - The parameter p for the Charlier polynomials

     Output:
       CH - Matrix of size (Ord x N) containing the DCP values
    """
    # Initialize variables
    CH = np.zeros((N, N), dtype=float)

    # Force `Parameter_p` to integer value and compute initial condition
    aa = int(np.floor(Parameter_p))
    Parameter_p = aa
    cor = 0 # Python indexing for matrices starts from 0 (zero-based indexing)
    init = np.sqrt(np.exp(-aa + (aa - 1) * np.log(aa) - lgamma(aa)))

    # Set initial values
    CH[aa - 1 + cor, 0 + cor] = init
    CH[aa + cor, 0 + cor] = init

    # Compute CH(0, x) using CH(0, 0)
    for n in range(aa - 1, 0, -1):
        CH[n - 1 + cor, 0 + cor] = CH[n + + cor, 0 + cor] / np.sqrt(Parameter_p / (n))

    for n in range(aa + 1, N - 1 + 1):
        CH[n + cor, 0 + cor] = CH[n - 1 + cor, 0 + cor] * np.sqrt(Parameter_p / (n))

    # Compute CH(1, x) using CH(0, x)
    for n in range(1, N):
        CH[n + cor, 1 + cor] = CH[n + cor, 0 + cor] * (Parameter_p - n) / np.sqrt(Parameter_p)

    # Compute higher-order Charlier polynomials
    for x in range(1, Ord - 1):
        for n in range(x + 1, N - 1 + 1):
            c1 = (Parameter_p - n + x) * np.sqrt(1 / (Parameter_p * (x + 1)))
            c2 = -np.sqrt(x / (x + 1))
            CH[n + + cor, x + 1 + cor] = c1 * CH[n + cor, x + cor] + c2 * CH[n + cor, x - 1 + cor]

    # Symmetry property: CH(x, n) = CH(n, x)
    for n in range(1, N):
        for x in range(n):
            CH[x + cor, n + cor] = CH[n + cor, x + cor]

    R=CH[:Ord, :N]
    return R



# Example usage:
# Parameter setting
N = 1500
poly_parameter = N/6
Order = N
# Call to generate Charlier Polynomials Matrix
R_matrix = DCP_(N, Order, poly_parameter)
