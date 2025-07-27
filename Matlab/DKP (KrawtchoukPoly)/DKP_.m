function R = DKP_(N,Ord, p)

% DKP_: Computes the Discrete Krawtchouk Polynomials (DKP)
%
% Inputs:
%   N   - Size of the signal (maximum x and n indices)
%   Ord - Maximum order of the polynomials (1 <= Ord <= N)
%   p   - Probability parameter (0 < p < 1)
%
% Output:
%   R - Matrix of size (Ord x N) containing the DKP values


    % Preallocate polynomial matrix
    R = double(zeros(N,Ord));
    
    % Adjust p to ensure it lies in [0, 0.5] for numerical stability
    p_signal=0;
    if p > 0.5
        p=1-p;
        p_signal=1;
    end
    
    %% Initial condition for n=0, x=0
    R(0+1,0+1)=((1.0-p)^((N-1.0)/2));
    
    %% Compute K(0,x) for x = 1 to N-1 
    for x=1:N-1
        R(0+1,x+1)= (sqrt((p*(N-x)))/sqrt((x*(1.0-p))))*R(0+1,x);
    end
    
    %% Compute K(1,x) for x = 1 to N-2 
    for x=1:N-2
        R(1+1,x+1)=((-x+p*(N-1.0))/(p*(N-1)))*sqrt(((N-1.0)*p)/(1.0-p))*R(0+1,x+1);
    end
    
    %% Compute K(n,x) using recurrence for n >= 2 and central range of x
    for n = 2:(ceil(N/2))
        for x = n:N-n-1
            an=((N-2*n+1)*p+n-x-1)/sqrt(p*n*(1-p)*(N-n));
            bn=sqrt((n-1)*(N-n+1)/(n*(N-n)));
            R(n+1,x+1)=an*R(n-1+1,x+1)-bn*R(n-2+1,x+1);
        end
    end
    
    %% Use symmetry to compute upper triangle from lower triangle
    for x=0:floor(N/2)-1
        for n=x+1:N-1-x
            R(n+1,x+1)=R(x+1,n+1);
        end
    end
    
    %% Compute K(n,x) for remaining positions using anti-symmetry
    for n=1:N-0
        for x=N-n+2:N-0
            R(n,x)=(-1)^(N+n+x-1)*R(N-n+1,N-x+1);
        end
    end
    
    %% Adjust output if p was originally > 0.5
    if p_signal==1
        R=fliplr(R);
        for i=1:N
            R(i,:)=(-1)^(i+1) .* R(i,:);
        end
    end
    
    
    %% Return only required order
    if Ord<N
		R = R(1:Ord, :);  
    end
    
