function t = DTP_(N, Ord)
	% DTP_: Computes the Discrete Tchebichef Polynomials (DTP)
	%
	% Inputs:
	%   N   - Size of the signal (maximum x and n indices)
	%   Ord - Maximum order of the polynomials (1 <= Ord <= N)
	%
	% Output:
	%   R - Matrix of size (Ord x N) containing the DKP values

    t = zeros(Ord, N);
    t(0+1, 0+1) = 1.0 / sqrt(N);
    
    ss=-sqrt((N-(1:Ord-1)) ./ (N+(1:Ord-1))) .* sqrt((2*(1:Ord-1)+1) ./ (2*(1:Ord-1)-1));
    % Calculate t(n, 0) and t(n, 1)
    for n = 1:Ord-1
        t(n+1, 0+1) = ss(n) * t(n, 1);
    end
    
    for n = 0:Ord-1
        t(n+1, 1+1) = (1 + n*(n+1) / (1-N)) * t(n+1, 1);
    end
    
    Ord1 = min(floor(N/4), Ord);
    
    % Calculate t(n, x) for x = 2:(N/2)-1
    for x = 2:(N/2)-1
        for n = 0:Ord1-1
            b1 = -n*(n+1) - (2*x-1)*(x-N-1) - x;
            lamda1 = b1 / (x*(N-x));
            lamda2 = ((x-1)*(x-N-1)) / (x*(N-x));
            t(n+1, x+1) = lamda1 * t(n+1, x) + lamda2 * t(n+1, x-1);
        end
    end
    
    % Calculate t(n, x) for x = (N/2):(N-1)
    for n = Ord1:Ord-1
        xx = mod(floor(0.5*N - sqrt((N*0.5)^2 - (n/2)^2)), Ord);
        for x = xx:N/2-1
            a1 = (2/n) * sqrt((4*n^2 - 1) / (N^2 - n^2));
            a2 = ((1-N)/(n)) * sqrt((4*n^2 - 1) / (N^2 - n^2));
            a3 = ((1-n)/(n)) * sqrt((2*n+1) / (2*n-3)) * sqrt((N^2 - (n-1)^2) / (N^2 - n^2));
            t(n+1, x+1) = a1 * x * t(n-1+1, x+1) + a2 * t(n-1+1, x+1) + a3 * t(n-2+1, x+1);
        end
    end
    
    % Calculate t(n, x) for x = (N/2)-1:-1:2
    for n = Ord-1:-1:Ord1
        xx = mod(floor(0.5*N - sqrt((N*0.5)^2 - (n/2)^2)), Ord);
        for x = xx+1:-1:2
            b1 = -n*(n+1) - (2*x-1)*(x-N-1) - x;
            lamda1 = b1 / (x*(N-x));
            lamda2 = ((x-1)*(x-N-1)) / (x*(N-x));
            t(n+1, x-1) = (t(n+1, x+1) - lamda1 * t(n+1, x)) / lamda2;
            if abs(t(n+1, x-1)) > abs(t(n+1, x))
                break
            end
        end
    end
    
    % Calculate t(n, x) for x = (N/2):N-1
    for x = (N/2):N-1
        for n = 0:Ord-1
            t(n+1, x+1) = t(n+1, N-1-x+1) / (-1)^n;
        end
    end
end