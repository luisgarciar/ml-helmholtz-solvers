function F = bw_F_poly(lambda_m, phi, eps_thick, s, b, NMAX)
% function F = BW_F_POLY(lambda_m, phi, eps_thick, s, b, NMAX)
% Computation of Faber polynomials associated with a scaled and translated
% bratwurst set, using a 3-term recursion.
% 
% Let Omega_eps be the bratwurst set determined with lambda_m, phi and
% eps_thick (cf. [1]).  The Faber polynomials F_1, ..., F_N associated with
% s*Omega_eps + b are computed using a 3-term recurrence introduced in [2].
% 
% Input:
%   lambda_m, phi, eps_thich -- parameters for bratwurst set Omega_eps
%   s    -- scaling of Omega_eps
%   b    -- translation of s*Omega_eps
%   NMAX -- number of Faber polynomials to be computed (F_1, ..., F_NMAX)
% 
% 
% Output:
%   F = (NMAX+1)x1 cell array containing F_0, F_1, ..., F_NMAX
%       (note F_0(z)=1).
% 
% 
% [1] T. Koch, J. Liesen, The conformal 'bratwurst' maps and associated
%     Faber polynomials, Numer. Math. (2000) 86: 173-191.
% [2] J. Liesen, Faber Polynomials Corresponding to Rational Exterior
%     Mapping Functions, Constr. Approx. (2001) 17: 267-274.
% 

% Computation of maximal thickness of the 'bratwurst' set
eps_max = tan(phi/4) * (1 + tan(phi/8));

% Parameter check
assert(phi > 0,'phi must be in ]0,2*pi[')
assert(phi < 2*pi,'phi must be in ]0,2*pi[')
assert(eps_thick >= 0,'eps_thick must be nonengative')
assert(eps_thick < eps_max,'eps_thick must be smaller than eps_max=%d',...
    eps_max)
assert(abs(lambda_m) == 1,'lambda_m must be on the unit circle')

% Computation of additional parameters
P = tan(phi/4) + 1/(cos(phi/4));
N = 0.5 * ( P/(1+eps_thick) + (1+eps_thick)/P );
M = ( (1+eps_thick)^2 - 1)/( 2 * (1+eps_thick) * tan(phi/4) );


% Coefficients of psi(w) = (w^2 + mu1*w + mu0)/(nu1*w + nu0)
mu1 = ((N-M)*b/s) - lambda_m*(M+N);
mu0 = (lambda_m^2)*M*N + lambda_m*(M*N - 1)*b/s;
nu1 = (N - M)/s;
nu0 = (N*M - 1)*lambda_m/s;


% shifted Faber polynomials hatF = { hatF_0, hatF_1, ..., hatF_N }
F = cell(NMAX+1,1);
F{1} = 2; % hatF_0(z) = 2;
F{2} = [nu1, -mu1]; % hatF_1(z) = nu1*z-mu1;
for n=3:NMAX+1
    F{n} = polyadd(conv([nu1, -mu1],F{n-1}), conv([nu0, -mu0],F{n-2}));
end

% Faber polynomials F_k for s*Omega_eps + b:
% F = { F_0, F_1, ..., F_NMAX },  note that F_0(z) = 1.
% F{1} = [1];
for n=1:NMAX+1
    F{n}(end) = F{n}(end)-(-nu0/nu1)^(n-1);
end

end
