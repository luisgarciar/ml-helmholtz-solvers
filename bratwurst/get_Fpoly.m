function F = get_Fpoly(kmax, M, N)
% function Fk = GET_FPOLY(kmax, M, N)
% Compute the Faber polynomials F_0, F_1, ..., F_kmax for a
% certain scaled and shifted bw set (lambda = -1 fixed).
% 
% See [S., Approximating the inverse on bw set] for details.

mu1 = 2*N ;
mu0 = 1 ;
nu1 = 2*(N-M) ;
nu0 = 2*(1-M*N) ;

W2 = [nu1, -mu1] ;
% V = [2*M*N-2, 1] ;
V = [-nu0, mu0] ;

% Computed shifted Faber polynomials (3-term recurrence)
F = cell(kmax+1, 1) ;
F{1} = [2] ;  %% Fhat_0
if kmax >= 2
    F{2} = W2 ; %% Fhat_1
end
for kk = 2:kmax
    F{kk+1} = polyadd( conv(W2, F{kk}), -conv(V, F{kk-1}) ) ; %% Fhat_kk=W2*F{kk}-V*F{kk-1}
end

% Compute the Faber polynomials
for kk=1:kmax+1
    F{kk}(end) = F{kk}(end) - (( -nu0/nu1 )^(kk-1));
end

end
