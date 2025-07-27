function R = DHP_(N, Ord, alpha, beta)

% DHP_: Computes the Discrete Hahn Polynomials (DHP)
%
% Inputs:
%   N      - Size of the signal (must be an even number)
%   Ord    - Maximum order of the polynomials (1 <= Ord <= N)
%   alpha  - Hahn polynomial parameter ?
%   beta   - Hahn polynomial parameter ?
%
% Output:
%   R - Matrix of size (Ord x N) containing the Hahn polynomial values
%       computed using a recurrence relation in both n and x directions

R=zeros(N,N); % Preallocate 2D matrix to store Hahn polynomial values
M1=floor(0.4*N); % Heuristic threshold for adaptive computation
M=min(M1,Ord); % To compute HPL to the required order if Ord < M;

%% Compute initial value h(0,0). Note that in Matlab, the indexing starts from 1 not from 0.
R(0+1,0+1)=exp((gammaln(alpha+beta+2)+gammaln(N+alpha)-(gammaln(alpha+1)+gammaln(N+alpha+beta+1)))/2);

%% compute value of HLP h(n,0); n=1,2,3,...
for n=0:M-2
    R(n+2,0+1)=-sqrt((2*n+3+alpha+beta)*(n+beta+1)*(n+alpha+beta+1)*(-n-1+N)/((n+1+alpha+beta+N)*(2*n+alpha+beta+1)*(n+1+alpha)*(n+1))) * R(n+1,0+1);
end

%% Compute HPL h(n,1); n=0,1,2,...
for n=0:M-1
    R(n+1,1+1)= (((N-n-1)*beta-n^2-(alpha+1)*n+N-1)*sqrt((beta+1)*(N-1)/(alpha+N-1))/((beta+1)*(N-1))) * R(n+1,0+1);
end

%% Compute DHPCs values in part P1 using x-direction recurrence algorithm for x=2,3,... and n=0,1,..,M
for x=2:(N/2)-1
    B=sqrt((N-x)*(beta+x)*(N+alpha-x)*x);
    AA=(-2*x^2+(2*N+alpha-beta+2)*x+(beta-1)*N-alpha-1)/B;
    CC=- sqrt((beta+x-1)*(N-x+1)*(x-1)*(N+alpha-x+1))/B;
    for n=0:M-1
        BB=-n*(alpha+beta+n+1)/B;
        R(n+1,x+1)= (AA+BB) * R(n+1,x) + CC * R(n+1,x-1);
    end
end

%% Compute DHPCs values in part P3 using both n and x direction in the range n=M, M+1,...,Ord and x=N/2-1,N/2-2, ... , adaptive location.
%  The adaptive location is computed using criteria when using the
%  x-direction recurrence relation.
for n=M:Ord-1
    DA=sqrt((alpha+beta+2*n+1)*(alpha+beta+2*n-1)*(alpha+beta+2*n)^2/(n*(alpha+beta+n)*(N-n)*(alpha+n)*(beta+n)*(alpha+beta+N+n)));
    CEA=-sqrt((n-1)*(alpha+n-1)*(beta+n-1)*(N-n+1)*(alpha+beta+n-1)*(alpha+beta+2*n+1)*(alpha+beta+N+n-1)*(alpha+beta+2*n)^2/(n*(alpha+beta+n)*(alpha+n)*(beta+n)*(N-n)*(alpha+beta+2*n-3)*(alpha+beta+N+n)*(alpha+beta+2*n-2)^2));
    for x=N/2-2:N/2-1
        B = x - ((alpha-beta+2*N-2)/(4)) - ( (beta^2 - alpha^2)*(alpha+beta+2*N) )/( 4*(alpha+beta+2*n-2)*(alpha+beta+2*n) );
        R(n+1,x+1)= (B*DA) * R(n,x+1) + (CEA)*R(n-1,x+1);
    end
    f=0;
    for x=N/2-2+1:-1:2
        B=sqrt((N-x)*(beta+x)*(N+alpha-x)*x);
        AA=(-2*x^2+(2*N+alpha-beta+2)*x+(beta-1)*N-alpha-1)/B;
        CC=- sqrt((beta+x-1)*(N-x+1)*(x-1)*(N+alpha-x+1))/B;
        BB=-n*(alpha+beta+n+1)/B;
        R(n+1,x-1)=  ( R(n+1,x+1) - ((AA+BB)) * R(n+1,x))/CC;
        if abs(R(n+1,x-1))>abs(R(n+1,x)) && abs(R(n+1,x-1))<1e-6
            if f==0
                f=1;
            else
                R(n+1,x-1)=0;
                break
            end
        end
    end
end

%% x= N/2, ..., N-1
R(0+1, N-1 +1)=exp((gammaln(beta+N)+gammaln(1+alpha)-(gammaln(1+beta)+gammaln(alpha+N)))/2) * R(0 +1,0 +1);
%% compute value of HLP h(n,N-1); n=1,2,3,...
for n=0:M-2
    R(n+2, N-1 +1)=sqrt((2*n+3+alpha+beta)*(n+alpha+beta+1)*(n+1+alpha)*(-n-1+N)/((n+beta+1)*(n+1+alpha+beta+N)*(2*n+alpha+beta+1)*(n+1))) * R(n+1, N-1 +1);
end

%% Compute DHPCs values h(n,1); n=0,1,2,...
for n=0:M-1
    R(n+1,N-2 +1)=((-n^2+(-alpha-beta-1)*n+(N-1)*(alpha+1))/sqrt((beta+N-1)*(N-1)*(alpha+1))) * R(n+1,N-1 +1);
end

%% Compute DHPCs values in part P2 using x-direction recurrence algorithm for x=2,3,... and n=0,1,..,M
for x=N-1:-1:(N/2)
    B=sqrt((N-x)*(beta+x)*(N+alpha-x)*x);
    AA=(-2*x^2+(2*N+alpha-beta+2)*x+(beta-1)*N-alpha-1)/B;
    CC=- sqrt((beta+x-1)*(N-x+1)*(x-1)*(N+alpha-x+1))/B;
    for n=0:M-1
        BB=-n*(alpha+beta+n+1)/B;
        R(n+1,x-1)=  ( R(n+1,x+1) - ((AA+BB)) * R(n+1,x))/CC;
    end
end

%% Compute DHPCs values in part P4 using both n and x direction in the range n=M, M+1,...,Ord and x=N/2-1,N/2-2, ... , adaptive location.
%  The adaptive location is computed using criteria when using the
%  x-direction recurrence relation.
for n=M:Ord-1
    DA=sqrt((alpha+beta+2*n+1)*(alpha+beta+2*n-1)*(alpha+beta+2*n)^2/(n*(alpha+beta+n)*(N-n)*(alpha+n)*(beta+n)*(alpha+beta+N+n)));
    CEA=-sqrt((n-1)*(alpha+n-1)*(beta+n-1)*(N-n+1)*(alpha+beta+n-1)*(alpha+beta+2*n+1)*(alpha+beta+N+n-1)*(alpha+beta+2*n)^2/(n*(alpha+beta+n)*(alpha+n)*(beta+n)*(N-n)*(alpha+beta+2*n-3)*(alpha+beta+N+n)*(alpha+beta+2*n-2)^2));
    for x=N/2:N/2+1
        B = x - ((alpha-beta+2*N-2)/(4)) - ( (beta^2 - alpha^2)*(alpha+beta+2*N) )/( 4*(alpha+beta+2*n-2)*(alpha+beta+2*n) );
        R(n+1,x+1)= (B*DA) * R(n,x+1) + (CEA)*R(n-1,x+1);
    end
    
    f=0;
    for x=N/2:N-1
        B=sqrt((N-x)*(beta+x)*(N+alpha-x)*x);
        AA=(-2*x^2+(2*N+alpha-beta+2)*x+(beta-1)*N-alpha-1)/B;
        CC=- sqrt((beta+x-1)*(N-x+1)*(x-1)*(N+alpha-x+1))/B;
        BB=-n*(alpha+beta+n+1)/B;
        R(n+1,x+1)= (AA+BB) * R(n+1,x) + CC * R(n+1,x-1);
        if abs(R(n+1,x+1))>abs(R(n+1,x)) && abs(R(n+1,x+1))<1e-6
            if f==0
                f=1;
            else
                R(n+1,x-1)=0;
                break
            end
        end
    end
end

%% Return only required order
if Ord<N
    R = R(1:Ord, :);
end