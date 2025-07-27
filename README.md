# Discrete Orthogonal Polynomial Matrix Generators

This repository provides implementations for generating six families of discrete orthogonal polynomials using three programming languages: **C++**, **Python**, and **MATLAB**. These polynomials are widely used in signal processing, pattern recognition, and numerical computation.

---

## 📚 Implemented Discrete Orthogonal Polynomial Families

- **DHP** – Discrete Hahn Polynomials, implemented based on the algorithm described in [1].
- **DMP** – Discrete Meixner Polynomials, implemented based on the algorithm described in [2].
- **DCP** – Discrete Charlier Polynomials, implemented based on the algorithm described in [3].
- **DKP** – Discrete Krawtchouk Polynomials, implemented based on the algorithm described in [4].
- **DTP** – Discrete Tchebichef Polynomials, implemented based on the algorithm described in [5].
- **DRP** – Discrete Racah Polynomials, implemented based on the algorithm described in [6].

---
## 📁 Folder Structure by Language
The repository is organized into three main language folders: C++, Python, and MATLAB. Each language folder contains subfolders representing different families of discrete orthogonal polynomials implemented in that language. The folder names use abbreviations for the polynomial families along with their full names in parentheses for clarity.

This structure allows easy navigation and access to specific polynomial implementations in the language of your choice:

```text
.
├── cpp/
│   ├── DHP (HahnPoly)/          # Discrete Hahn Polynomials in C++
│   ├── DMP (MeixnerPoly)/       # Discrete Meixner Polynomials in C++
│   ├── DCP (CharlierPoly)/      # Discrete Charlier Polynomials in C++
│   ├── DKP (KrawtchoukPoly)/    # Discrete Krawtchouk Polynomials in C++
│   ├── DTP (TchebichefPoly)/    # Discrete Tchebichef Polynomials in C++
│   └── DRP (RacahPoly)/         # Discrete Racah Polynomials in C++

├── python/
│   ├── DHP (HahnPoly)/          # Discrete Hahn Polynomials in Python
│   ├── DMP (MeixnerPoly)/       # Discrete Meixner Polynomials in Python
│   ├── DCP (CharlierPoly)/      # Discrete Charlier Polynomials in Python
│   ├── DKP (KrawtchoukPoly)/    # Discrete Krawtchouk Polynomials in Python
│   ├── DTP (TchebichefPoly)/    # Discrete Tchebichef Polynomials in Python
│   └── DRP (RacahPoly)/         # Discrete Racah Polynomials in Python

├── matlab/
│   ├── DHP (HahnPoly)/          # Discrete Hahn Polynomials in MATLAB
│   ├── DMP (MeixnerPoly)/       # Discrete Meixner Polynomials in MATLAB
│   ├── DCP (CharlierPoly)/      # Discrete Charlier Polynomials in MATLAB
│   ├── DKP (KrawtchoukPoly)/    # Discrete Krawtchouk Polynomials in MATLAB
│   ├── DTP (TchebichefPoly)/    # Discrete Tchebichef Polynomials in MATLAB
│   └── DRP (RacahPoly)/         # Discrete Racah Polynomials in MATLAB
```
## 💻 Language Implementations

### 1. 🔷 C++

#### 1.1. Contents

This section describes the C++ source files included in the repository. Each `.cpp` file implements a specific discrete orthogonal polynomial family generator function. The primary polynomial computation function in each file is named with an underscore suffix (e.g., `DCP_()`), and a `main()` function is provided for demonstration, testing, and performance measurement.

- `DHP.cpp` contains the function `DHP_()` – Discrete Hahn polynomial generator 
- `DMP.cpp` contains the function `DMP_()` – Discrete Meixner polynomial generator  
- `DCP.cpp` contains the function `DCP_()` – Discrete Charlier polynomial generator  
- `DKP.cpp` contains the function `DKP_()` – Discrete Krawtchouk polynomial generator  
- `DTP.cpp` contains the function `DTP_()` – Discrete Tchebichef polynomial generator  
- `DRP.cpp` contains the function `DRP_()` – Discrete Racah polynomial generator  
- Each file includes a `main()` function used for demonstration and testing purposes, including performance timing
---
#### 1.2. Build Instructions

Compile with any C++17-compliant compiler. Example using `g++`:

```bash
g++ -std=c++17 -O2 poly_file_name.cpp -o generated_file_name
````

---

#### 1.3. ⏱ Sample Benchmark

The following is an example demonstrating how to benchmark the computation time of the `DRP_()` function, which generates the Discrete Racah Polynomials matrix. 

**Note:** The input parameters to the function (such as `N`, `Order`, `alpha`, `beta`, and `a`) must be carefully chosen according to the specific polynomial family and use case, as they directly affect the output and performance.

```cpp
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
```
---

#### 1.4. 📋 Requirements


- C++17 or later compliant compiler (e.g., `g++` 7+, `clang++` 6+)
- Standard C++ libraries:
  - `<iostream>`, `<vector>`, `<cmath>`, `<chrono>`, `<fstream>`, `<algorithm>`, `<stdexcept>`
- Optional third-party libraries (depending on polynomial implementation):
  - SymEngine — symbolic computation (required for DMP)
  - GNU Scientific Library (GSL) — special functions (required for DMP)
  - Eigen — matrix operations (required for DMP)
  - BLAS/LAPACK — for optimized numerical routines (used by Armadillo/Eigen)

---


### 2. 🐍 Python

This repository offers efficient Python implementations for generating six families of discrete orthogonal polynomials widely used in signal processing, numerical analysis, and pattern recognition.

#### 2.1. Contents

This section describes the Python source files included in the repository. Each `.py` file implements a specific discrete orthogonal polynomial family generator function. The primary polynomial computation function in each file is named according to the polynomial family (e.g., `dcp()`), and each file contains example code for demonstration, testing, and performance measurement.

- `dhp.py` contains the function `dhp_()` – Discrete Hahn polynomial generator  
- `dmp.py` contains the function `dmp_()` – Discrete Meixner polynomial generator  
- `dcp.py` contains the function `dcp_()` – Discrete Charlier polynomial generator  
- `dkp.py` contains the function `dkp_()` – Discrete Krawtchouk polynomial generator  
- `dtp.py` contains the function `dtp_()` – Discrete Tchebichef polynomial generator  
- `drp.py` contains the function `drp_()` – Discrete Racah polynomial generator  
- Each file includes example usage code and test routines, including performance timing

#### 2.2. Run Instructions

Use Python 3.8.5 or newer Example usage:

```bash
cd "DRP (HahnPoly)"
python dhp.py
````

---
#### 2.3. ⏱ Sample Benchmark

The following is an example demonstrating how to benchmark the computation time of the `drp()` function, which generates the Discrete Racah Polynomials matrix.

**Note:** The input parameters (such as `N`, `a`, `alpha`, and `beta`) must be chosen carefully according to the specific polynomial family and use case, as they directly affect the output and performance.

```python

# Parameter setting
N = 1600
a = int(N / 2)
alpha = int(N / 4)
beta = int(N / 8)

# Call to generate Discrete Racah Polynomials matrix
R_matrix = drp(N, a, alpha, beta)

```
---

#### 2.4. 📋 Requirements

- Python 3.8.5 or later  
- Required Python packages (depending on polynomial implementation):
  - `numpy`
  - `math` (standard library)
  - `time` (standard library)
  - `matplotlib` (for DHP)
  - `pandas` (for DMP)
  - `sympy` (symbolic computations, required for DMP)

---

**Note:**  
Make sure to install all required packages before running the scripts, e.g.:

```bash
pip install numpy matplotlib pandas sympy

```
---

### 3. 📊 MATLAB

This repository provides MATLAB implementations for generating six families of discrete orthogonal polynomials. These polynomials are widely used in signal processing, numerical analysis, and pattern recognition applications.

#### 3.1. Content
This section describes the Matlab source files included in the repository. Each `.m` file implements a specific discrete orthogonal polynomial family generator function. The main polynomial computation function in each file is named according to the polynomial family with an underscore suffix (e.g., `dcp_()`), and each file contains example scripts for demonstration, testing, and performance measurement.

| File Name     | Description                                 |
|---------------|---------------------------------------------|
| `DHP_.m`      | Computes Discrete Hahn Polynomial matrix.   |
| `DMP_.m`      | Computes Discrete Meixner Polynomial matrix.|
| `DCP_.m`      | Computes Discrete Chebyshev Polynomial matrix.|
| `DKP_.m`      | Computes Discrete Krawtchouk Polynomial matrix.|
| `DTP_.m`      | Computes Discrete Tchebichef Polynomial matrix.|
| `DRP_.m`      | Computes Discrete Racah Polynomial matrix.   |

---

#### 3.2. Run Instructions

Each function follows this typical structure:

```matlab
M = <POLY>_(N, Ord, param1, param2, ...)
```

- `N`     – Number of discrete points.
- `Ord`   – Maximum order of the polynomials (rows in the output matrix).
- `paramX` – Family-specific parameters.

#### 3.3. ⏱ Sample Benchmark

To run a polynomial generator, navigate to the directory containing the `.m` files and run the corresponding function from the MATLAB command window. For example:

```matlab
% Discrete Racah Polynomial Example
N = 4000;
a = N / 2;
alpha = N / 4;
beta = N / 8;

H = DRP_(N, a, alpha, beta);
````
---

#### 3.4. 📋 Requirements

- Matlab 2019b or later  
---

## 📚 References
[1] B. M. Mahmmod, S. H. Abdulhussain, T. Suk, and A. Hussain, “Fast Computation of Hahn Polynomials for High Order Moments,” IEEE Access, vol. 10, pp. 48719–48732, 2022.

[2] S. H. Abdulhussain and B. M. Mahmmod, “Fast and efficient recursive algorithm of Meixner polynomials,” J. Real-Time Image Process., vol. 18, no. 6, pp. 2225–2237, Dec. 2021.

[3] A. M. Abdul-Hadi, S. H. Abdulhussain, and B. M. Mahmmod, “On the computational aspects of Charlier polynomials,” Cogent Eng., vol. 7, no. 1, p. 1763553, Jan. 2020.

[4] S. H. Abdulhussain, A. R. Ramli, S. A. R. Al-Haddad, B. M. Mahmmod, and W. A. Jassim, “Fast Recursive Computation of Krawtchouk Polynomials,” J. Math. Imaging Vis., vol. 60, no. 3, pp. 285–303, Mar. 2018.

[5] S. H. Abdulhussain, A. R. Ramli, S. A. R. Al-Haddad, B. M. Mahmmod, and W. A. Jassim, “On Computational Aspects of Tchebichef Polynomials for Higher Polynomial Order,” IEEE Access, vol. 5, no. 1, pp. 2470–2478, 2017.

[6] B. M. Mahmmod, S. H. Abdulhussain, T. Suk, M. Alsabah, and A. Hussain, “Accelerated and Improved Stabilization for High Order Moments of Racah Polynomials,” IEEE Access, vol. 11, pp. 110502–110521, 2023.


## 📄 License
This project is licensed under the MIT License. See the [LICENSE](./LICENSE) file for details.
