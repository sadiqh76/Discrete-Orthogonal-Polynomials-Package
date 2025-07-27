% DRP_: Computes the Discrete Racah Polynomials (DRP)
%
% Inputs:
%   N   - Size of the polynomial matrix (maximum degree + 1)
%   a   - Racah polynomial parameter (positive integer)
%   ph  - Racah parameter phi (must satisfy ph <= a)
%   bt  - Racah parameter beta (must satisfy bt <= ph)
%
% Output:
%   R   - Matrix of size (N x N) containing the DRP values
%
% Notes:
%   - The function uses recurrence relations and gamma functions
%     to construct the matrix of Racah polynomials.
%   - Numerical stability is managed via thresholds and dynamic index bounds.
%   - Ensure the input parameters satisfy: 0 < bt ? ph ? a.

function R=DRP_(N,a,alpha,beta)

N_max=N;    % Maximum number of polynomial orders
b=a+N;      % Upper limit for b in Racah polynomial parameters

thres=1e-5;  % Threshold for pruning small values (stability)

R=zeros(N,N); % Preallocate output matrix
cor = 1;     % Index correction for MATLAB's 1-based indexing

%% Initial condition using logarithmic gamma for numerical stability
FF1=gammaln(alpha+beta+2)+gammaln(2*a+N)+gammaln(beta+N)+gammaln(2*a+2*N+alpha)-...
    (gammaln(2*a+2*N-1)+gammaln(beta+1)+gammaln(alpha+beta+N+1)+gammaln(2*a+N+alpha+1));
F1=exp(FF1/2);
R(0 +cor,N-1 +cor)=F1;

%% First row (n=0) backwards in s
for s=a+N-2:-1:a
    R(0 +cor,s-a +cor)= ...
        sqrt((2*s+1)*(a-beta+s+1)*(b+s+1)*(b+alpha-s-1)*(a-s-1)/((a+s+1)*(b+alpha+s+1)*(a-beta-s-1)*(2*s+3)*(b-s-1))) * ...
        R(0 +cor,s-a+1 +cor);
end

%% Second row (n=1) backwards in s
for s=a:a+N-1
    R(1 +cor,s-a +cor) = ...
        -(((-a+b-1)*alpha+b^2-s^2-a-s-1)*beta+(a^2-s^2+b-s-1)*alpha+a^2+b^2-2*(s^2+s)-1) * ...
        sqrt((alpha+beta+3)/((a-b+1)*(a+b-beta-1)*(alpha+1)*(beta+1)*(a-b-alpha-beta-1)*(a+b+alpha+1))) * R(0 +cor,s-a +cor);
end

%% Fill last column (s = N-1) and first column (s = a) for all n
s=N-1;
for n=1:N_max-2
    R(n+1 +cor,s +cor) = ...
        sqrt((N-n-1)*(alpha+beta+2*n+3)*(alpha+beta+n+1)*(alpha+n+1)*(2*a+N-beta-n-1)/((alpha+beta+2*n+1)*(beta+n+1)*(N+alpha+beta+n+1)*(n+1)*(2*a+N+alpha+n+1))) * ...
        R(n +cor,s +cor);
    R(n+1 +cor,cor) = ...
        -sqrt((N-n-1)*(alpha+beta+2*n+3)*((alpha+beta+n+1)*(beta+n+1))*(2*a+N+alpha+n+1)/((2*a+N-beta-n-1)*(alpha+beta+2*n+1)*(alpha+n+1)*(N+alpha+beta+n+1)*(n+1))) * ...
        R(n +cor,cor);
end

%% Identify dominant rows for improved numerical propagation
[~, maxind_N_a_1]=max(R(:,end));
[~, maxind_0]=max(R(:,1));

%% Handle potential instability or failure
if maxind_N_a_1==1 && N>4
    disp(['Unable to compute initial values for a=',num2str(a),' b=',num2str(b),' ph=',num2str(alpha),' bt=',num2str(beta)])
    disp('Please terminate')
end

if maxind_0==1
    maxind_0=maxind_N_a_1;
end

%% Bound adjustments to avoid indexing errors
maxind_N_a_1=min(maxind_N_a_1,N_max-1);
maxind_0=min(maxind_0,N_max-1);
if maxind_N_a_1<=1
    maxind_N_a_1=2;
end

%% Part 1: Compute polynomials for small n (low degrees)
for s=a+1:a+N-2
    mis=round(maxind_0+(s-a)*(maxind_N_a_1-maxind_0)/(N-1));
    for n=2:mis
        A=(n*(alpha+beta+n))/((alpha+beta+2*n-1)*(alpha+beta+2*n));
        B=s.*(s+1) - ...
            (1/4)*(a^2+b^2+(a-beta)^2+(b+alpha)^2 - 2) + ...
            ((1/8)*((alpha+beta+2*n-2)*(alpha+beta+2*n))) - ...
            (((beta^2-alpha^2)*((b+alpha/2)^2-(a-beta/2)^2))/...
            (2*(alpha+beta+2*n-2)*(alpha+beta+2*n)));
        C = -(alpha+n-1)*(beta+n-1) * ...
            ((a+b+(alpha-beta)/2)^2-(n-1+(alpha+beta)/2)^2)*...
            ((b-a+(alpha+beta)/2)^2-(n-1+(alpha+beta)/2)^2) / ...
            ((alpha+beta+2*n-2)*(alpha+beta+2*n-1));
        B1=(alpha+beta+2*n+1)*n*(alpha+beta+n)/((alpha+beta+2*n-1)*(-b+a+n)*(a+b-beta-n)*(alpha+n)*(beta+n)*(-b+a-alpha-beta-n)*(a+b+alpha+n));
        C1=(n-1)*n*(alpha+beta+n-1)*(alpha+beta+n)*(alpha+beta+2*n+1)/((alpha+n-1)*(alpha+n)*(beta+n-1)*(beta+n)*(-b+a-alpha-beta-n+1)*...
            (-b+a-alpha-beta-n)*(a+b+alpha+n-1)*(a+b+alpha+n)*(alpha+beta+2*n-3)*(-b+a+n)*(-b+a+n-1)*(a+b-beta-n)*(a+b-beta-n+1));
        P1=(B/A)*sqrt(B1);
        P2=(C/A)*sqrt(C1);
        R(n +cor,s-a +cor)=P1 .* R(n-1 +cor,s-a +cor)+P2 .* R(n-2 +cor,s-a +cor);
        
    end
end

%% Part 2: Higher n with thresholding to ensure numerical stability
for s=a+1:a + N - 2
    mis=round(maxind_0+(s-a)*(maxind_N_a_1-maxind_0)/(N-1));
    for n=mis+1:N_max-1
        A=(n*(alpha+beta+n))/((alpha+beta+2*n-1)*(alpha+beta+2*n));
        B=s.*(s+1) - ...
            (1/4)*(a^2+b^2+(a-beta)^2+(b+alpha)^2 - 2) + ...
            ((1/8)*((alpha+beta+2*n-2)*(alpha+beta+2*n))) - ...
            (((beta^2-alpha^2)*((b+alpha/2)^2-(a-beta/2)^2))/...
            (2*(alpha+beta+2*n-2)*(alpha+beta+2*n)));
        C = -(alpha+n-1)*(beta+n-1) * ...
            ((a+b+(alpha-beta)/2)^2-(n-1+(alpha+beta)/2)^2)*...
            ((b-a+(alpha+beta)/2)^2-(n-1+(alpha+beta)/2)^2) / ...
            ((alpha+beta+2*n-2)*(alpha+beta+2*n-1));
        B1=(alpha+beta+2*n+1)*n*(alpha+beta+n)/((alpha+beta+2*n-1)*(-b+a+n)*(a+b-beta-n)*(alpha+n)*(beta+n)*(-b+a-alpha-beta-n)*(a+b+alpha+n));
        C1=(n-1)*n*(alpha+beta+n-1)*(alpha+beta+n)*(alpha+beta+2*n+1)/((alpha+n-1)*(alpha+n)*(beta+n-1)*(beta+n)*(-b+a-alpha-beta-n+1)*(-b+a-alpha-beta-n)*(a+b+alpha+n-1)*(a+b+alpha+n)*(alpha+beta+2*n-3)*(-b+a+n)*(-b+a+n-1)*(a+b-beta-n)*(a+b-beta-n+1));
        P1=(B/A)*sqrt(B1);
        P2=(C/A)*sqrt(C1);
        R(n +cor,s-a +cor)=P1 .* R(n-1 +cor,s-a +cor)+P2 .* R(n-2 +cor,s-a +cor);
        if abs(R(n +cor,s-a +cor))<=thres && (abs(R(n +cor,s-a +cor))>abs(R(n-1 +cor,s-a +cor)))
            R(n +cor, s-a +cor)=0;
            break;
        end
    end
end
