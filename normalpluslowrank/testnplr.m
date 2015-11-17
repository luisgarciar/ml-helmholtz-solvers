%Script to test if a Helmholtz matrix A preconditioned by the Shifted
%Laplacian S is of the form inv(S)*A = N + F, with N normal and F a low
%rank matrix.

%1D Test
% Helmholtz matrix

f  = @dirac1d;
k  = 200;
np = ceil(5*k/pi)+1;
bc = 'som';
beta = 1-0.5*1i;
dim = 1;

[A,~,~] = helmholtzbc(dim,bc,np,k,f,0);
[M,~,~] = helmholtzbc(dim,bc,np,k*sqrt(beta),f,0);

S = M\A;

% eigv= eigs(S_s, length(S_s)-2);
% plot(eigv,'+'), axis([-1 2 -1 2]), shg

SCA   = A*A'-A'*A;    %Self-commutator of A
rkSCA  = sprank(SCA); %Rank of self-commutator of A

HA    = (A+A')/2;    %Hermitian part of A
sHA   = A-HA;         %Skew-Hermitian part of A
rksHA = sprank(sHA); %rank of skew hermitian part

SCS   = S*S'-S'*S;    %Self-commutator of S
rkSCS = sprank(SCS);  %rank of self commutator
rksHA = sprank(sHA);  %rank of skew Hermitian part

C = schur(full(S));
D = diag(diag(C));
T = C-D;




