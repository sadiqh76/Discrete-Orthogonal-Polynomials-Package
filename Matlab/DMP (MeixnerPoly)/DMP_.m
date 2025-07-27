% DMP_: Computes the Discrete Meixner Polynomials (DMP)
%
% Inputs:
%   N     - Size of the signal (number of samples)
%   Ord   - Maximum order of the polynomials (1 <= Ord <= N)
%   c     - The Meixner parameter c (0 < c < 1)
%   beta  - The Meixner parameter beta 
%
% Output:
%   R - Matrix of size (N x Ord) containing the DMP values

function R = DMP_(N, Ord, c, beta)
    %% Define constants and preallocate output matrix
    cor = 1; % Offset for MATLAB 1-based indexing
    R = zeros(N, Ord);

    %% Precompute symbolic function and convert to numeric for efficiency
    syms r(nn,xx,bb,cc)
    r(nn,xx,bb,cc)=pochhammer(bb, nn)*hypergeom([-nn, -xx], [bb], 1 - 1/cc)* ...
        sqrt(cc^xx*gamma(bb+xx)*cc^nn*(1-cc)^bb/(factorial(xx)* ...
        gamma(bb)*factorial(nn)*pochhammer(bb, nn)));

    %% Define key central indices
    n1 = N / 2; x1 = N / 2;
    n0 = N / 2 - 1; x0 = N / 2 - 1;

    %% Initialize central values
    R(n1 + cor, x1 + cor) = r(x1, n1, beta, c);
    R(n1 + cor, x0 + cor) = r(n1, x0, beta, c);
    R(n0 + cor, x1 + cor) = R(n1 + cor, x0 + cor);
    R(n0 + cor, x0 + cor) = r(n0, x0, beta, c);

    %% Computing Meixner coefficients using Recurrence Relations
    for x = x0:x1
        for n = n1 + 1:N - 1
            A = c / (c - 1);
            B = (x - x * c - n + 1 - n * c + c - beta * c) / (1 - c);
            C = (n - 1) * (n - 2 + beta) / (1 - c);
            D = sqrt(c / (n * (beta + n - 1)));
            E = sqrt(c^2 / (n * (n - 1) * (beta + n - 2) * (beta + n - 1)));
            R(n + cor, x + cor) = (B * D) / A * R(n -1 + cor, x + cor) + (C * E) / A * R(n - 2 + cor, x + cor);
            R(x + cor, n + cor) = R(n + cor, x + cor);
        end
    end

    for n = n0 - 1:N - 1
        for x = x1:n - 1
            c1 = (n * (c - 1) + x + (x + beta) * c) * sqrt(1 / (c * (x + 1) * (x + beta)));
            c2 = -sqrt((x * (x + beta - 1)) / ((x + 1) * (x + beta)));
            R(n + cor, x + 1 + cor) = c1 * R(n + cor, x + cor) + c2 * R(n + cor, x - 1 + cor);
            R(x + 1 + cor, n + cor) = R(n + cor, x + 1 + cor);
        end
    end

    for n = n1 - 1:N - 1
        x = x0;
        while x >= 1
            c1 = (n * (c - 1) + x + (x + beta) * c) * sqrt(1 / (c * (x + 1) * (x + beta)));
            c2 = -sqrt((x * (x + beta - 1)) / ((x + 1) * (x + beta)));
            R(n + cor, x - 1 + cor) = -c1 / c2 * R(n + cor, x + cor) + 1 / c2 * R(n + cor, x + 1 + cor);
            R(x - 1 + cor, n + cor) = R(n + cor, x - 1 + cor);
            if abs(R(n + cor, x - 1 + cor)) < 1e-7 && abs(R(n + cor, x + cor)) < 1e-5
                break;
            end
            x = x - 1;
        end
    end

    for x = n1 - 2:-1:0
        n = x + 2;
        while n >= 2
            A = c / (c - 1);
            B = (x - x * c - n + 1 - n * c + c - beta * c) / (1 - c);
            C = (n - 1) * (n - 2 + beta) / (1 - c);
            D = sqrt(c / (n * (beta + n - 1)));
            E = sqrt(c^2 / (n * (n - 1) * (beta + n - 2) * (beta + n - 1)));
            R(n - 2 + cor, x + cor) = A / (C * E) * R(n + cor, x + cor) - (B * D) / (C * E) * R(n - 1 + cor, x + cor);
            R(x + cor, n - 2 + cor) = R(n - 2 + cor, x + cor);
            if abs(R(n - 2 + cor, x + cor)) < 1e-6 && abs(R(n - 1 + cor, x + cor)) < 1e-4
                break;
            end
            n = n - 1;
        end
    end

    %% Return only required order
    if Ord < N
        R=R(1:Ord,:);
    end