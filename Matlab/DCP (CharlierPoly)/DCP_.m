function R = DCP_(N, Ord, Parameter_p)
    % DCP_: Computes the Discrete Charlier Polynomials (DCP)
    % 
    % Inputs:
    %   N          - Size of the signal (maximum n index)
    %   Ord        - Maximum order of the polynomials (1 <= Ord <= N)
    %   Parameter_p - The parameter p for the Charlier polynomials
    %
    % Output:
    %   R - Matrix of size (Ord x N) containing the DCP values

    % Preallocate memory and initialize variables
	R = double(zeros(Ord,N));
	
	%% Initial Values
	aa=floor(Parameter_p);
	Parameter_p=aa;
	cor = 1;
	init=sqrt(exp(-aa + (aa-1)*log(aa)-gammaln(aa)));

	%% Set initial values for CH(aa-1, 0) and CH(aa, 0)
	R(aa-1+cor,0+cor)=init;
	R(aa+cor,0+cor)=init;

	%% Compute CH(0,x) using CH(0,0) for n < aa -1
	for n=aa-1:-1:1 
		R(n-1+cor,0+cor)=R(n+cor,0+cor)/sqrt(Parameter_p/n);
	end

	%% Compute CH(0,x) using CH(0,0) for n > aa
	for n=aa+1:N-1 
		R(n+cor,0+cor)=R(n-1+cor,0+cor)*sqrt(Parameter_p/n);
	end

	%% Compute CH(1,x) using CH(0,x)
	for n=1:N-1 
		R(n+cor,1+cor)=R(n+cor,0+cor)*(Parameter_p-n)/sqrt(Parameter_p);
	end

	%% Recurrence for Charlier polynomials CH(n, x)
	for x = 1:Ord-2
		for n=x+1:N-1
			c1=(Parameter_p-n+x) * sqrt(1/(Parameter_p*(x+1)));
			c2=-sqrt(x/(x+1));
			R(n+cor,x+1+cor)=c1*R(n+cor,x+cor)+c2*R(n+cor,x-1+cor);
		end
	end

	%% Symmetrize matrix: Fill lower triangle from upper triangle
	for n = 1:N-1
		for x=0:n-1
			R(x+cor,n+cor)= R(n+cor,x+cor);
		end
    end
	
    %% Return only required order
	if Ord<N
		R = R(1:Ord, :);  
    end
    
   
