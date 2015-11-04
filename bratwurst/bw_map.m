function [psi, dpsi, capacity, M, N] = bw_map(lambda_m,phi,eps_thick)
% function [psi, dpsi, capacity, M, N] = BW_MAP(lambda_m,phi,eps_thick)
% Gives a function handle to the exterior mapping function of the
% 'bratwurst' shaped set with parameters lambda_m, phi and eps.
% 
% Input parameters:
%   lambda_m      -- on unit circle
%   phi           -- in ]0, 2pi[
%   eps_thickness -- in [0, eps_max[
% 
% We use the notation of Theorem 3.1 in [1].
% 
% Output:
%   psi  -- function handle of exterior mapping function
%   dpsi -- function handle of d psi / dz
%   capacity -- capacity of 'bratwurst' shaped set
%   M    -- see [1, Theorem 3.1]
%   N    -- see [1, Theorem 3.1]
% 
% 
% [1] J. Liesen, T. Koch, The conformal 'bratwurst' maps and associated
%     Faber polynomials, Numer. Math. (2000) 86: 173-191.
% 

% Computation of maximal thickness of the 'bratwurst' set
eps_max = tan(phi/4) * (1 + tan(phi/8));

% parameter check
assert(phi > 0,'phi must be in ]0,2*pi[')
assert(phi < 2*pi,'phi must be in ]0,2*pi[')
assert(eps_thick >= 0,'eps_thick must be nonengative')
assert(eps_thick < eps_max,'eps_thick must be smaller than eps_max=%d',eps_max)
assert(abs(lambda_m) == 1,'lambda_m must be on the unit circle')

sigma = 1+eps_thick ;

% Computation of additional parameters
P = tan(phi/4) + 1/cos(phi/4);
N = 0.5 * ( P/sigma + sigma/P );  %% N_eps
M = ( sigma^2 - 1)/( 2 * sigma * tan(phi/4) ); %% M_eps

psi = @(z) ((z - lambda_m*N).*(z - lambda_m*M))./( (N-M).*z + lambda_m*(N*M-1) );

dpsi = @(z) ( (2*z - lambda_m*(M+N)).*( (N-M).*z + lambda_m * (N*M-1)) ...
    - (z-lambda_m*M).*(z-lambda_m*N).*(N-M) )./ ...
    ( ( (N-M).*z + lambda_m*(N*M-1) ).^2 ) ;

capacity = 1/(N - M);
% Romega = 1/N;

end
