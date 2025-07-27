import numpy as np

def DTP_(N, Ord):
    """
    DTP_: Computes the Discrete Tchebichef Polynomials (DTP)

	 Inputs:
	   N   - Size of the signal (maximum x and n indices)
	   Ord - Maximum order of the polynomials (1 <= Ord <= N)

	 Output:
	   T - Matrix of size (Ord x N) containing the DKP values
    """

    # Initialize the matrix t with zeros
    T = np.zeros((Ord, N))

    # Set the first element
    T[0, 0] = 1.0 / np.sqrt(N)

    # Precompute the scaling factor ss
    ss = -np.sqrt((N - np.arange(1, Ord)) / (N + np.arange(1, Ord))) * \
         np.sqrt((2 * np.arange(1, Ord) + 1) / (2 * np.arange(1, Ord) - 1))

    # Calculate T(n, 0) and T(n, 1)
    for n in range(1, Ord):
        T[n, 0] = ss[n - 1] * T[n - 1, 0]

    for n in range(Ord):
        T[n, 1] = (1 + n * (n + 1) / (1 - N)) * T[n, 0]

    # Determine the limit for the first loop
    Ord1 = min(int(np.floor(N / 4)), Ord)

    # Calculate T(n, x) for x = 2:(N/2)-1
    for x in range(2, int(N / 2)):
        for n in range(Ord1):
            b1 = -n * (n + 1) - (2 * x - 1) * (x - N - 1) - x
            lamda1 = b1 / (x * (N - x))
            lamda2 = ((x - 1) * (x - N - 1)) / (x * (N - x))
            T[n, x] = lamda1 * T[n, x - 1] + lamda2 * T[n, x - 2]

    # Calculate T(n, x) for x = (N/2):(N-1)
    for n in range(Ord1, Ord):
        xx = int(np.floor(0.5 * N - np.sqrt((N * 0.5) ** 2 - (n / 2) ** 2))) % Ord
        for x in range(xx, int(N / 2)):
            a1 = (2 / n) * np.sqrt((4 * n ** 2 - 1) / (N ** 2 - n ** 2))
            a2 = ((1 - N) / n) * np.sqrt((4 * n ** 2 - 1) / (N ** 2 - n ** 2))
            a3 = ((1 - n) / n) * np.sqrt((2 * n + 1) / (2 * n - 3)) * np.sqrt((N ** 2 - (n - 1) ** 2) / (N ** 2 - n ** 2))
            T[n, x] = a1 * x * T[n - 1, x] + a2 * T[n - 1, x] + a3 * T[n - 2, x]

    # Calculate T(n, x) for x = (N/2)-1:-1:2
    for n in range(Ord - 1, Ord1 - 1, -1):
        xx = int(np.floor(0.5 * N - np.sqrt((N * 0.5) ** 2 - (n / 2) ** 2))) % Ord
        for x in range(xx + 1, 1, -1):
            b1 = -n * (n + 1) - (2 * x - 1) * (x - N - 1) - x
            lamda1 = b1 / (x * (N - x))
            lamda2 = ((x - 1) * (x - N - 1)) / (x * (N - x))
            T[n, x - 2] = (T[n, x] - lamda1 * T[n, x - 1]) / lamda2
            if abs(T[n, x - 2]) > abs(T[n, x - 1]):
                break

    # Calculate T(n, x) for x = (N/2):N-1
    for x in range(int(N / 2), N):
        for n in range(Ord):
            T[n, x] = T[n, N - 1 - x] / (-1) ** n

    return T

# Example usage:
# Parameter setting
N = 1000
Order = N

# Call to generate Tchebichef Polynomials Matrix
R = DTP_(N, Order)
