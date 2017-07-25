function [fovAhat,minfov]=slapfov(A,Aeps,fvpts)
%% SLAPFOV
%  [fovSL,min]= SLapfov(A,Aeps) computes the field of values of the
%  matrix S=A*inv(Aeps), where A is the discrete Helmholtz operator 
%  and Aeps is the shifted Laplacian with shift eps
%   
%  INPUT
%  A: Helmholtz sparse matrix
%  Aeps: shifted Laplace matrix

%  OUTPUT
%  fovSL : Field of values and eigenvalues of A*inv(Aeps)
%  m: Inner numerical radius of A*inv(Aeps) (distance of the fov to zero)
% 
% Author: Luis Garcia Ramos, TU Berlin
%         Based on a routine of Nick Higham (Matrix Toolbox)
%         version 0.1 - Jul 2017
%         
%%
[L,U] = lu(Aeps);
LH = U'; UH = L';
N  = length(A);

Ahat     = @(x) A*(U\(L\x));
AhatH    = @(x) UH\(LH\(A'*x));
H        = @(x) 0.5*(feval(Ahat,x) + feval(AhatH,x)); %Hermitian part of Ahat

% %Compute first the maximum eigenvalue of A
% [eigvecA, eigvalA] =  eig(full(A));
% [~,ind] = max(diag(eigvalA));
% vmaxA   = eigvecA(:,ind);

%Compute first the maximum eigenvalue of A (sparsely)
%opts.p      = 30;
opts.p      = 60;
opts.tol    = 1e-4;
opts.isreal = 0;
[vmaxA,~]   = eigs(A,1,'lr',opts);


%Compute now the maximum eigenvalue of the Hermitian part of A
opts.p      = 30;
opts.tol    = 1e-6;
opts.isreal = 0;
opts.v0     = vmaxA;
[vmaxH,~]   = eigs(H,N,1,'LM',opts);
    
%Use the maximum eigenvalue of the Hermitian part of A
[fovAhat,~,~] = sfov(Ahat,AhatH,vmaxH,N,fvpts);

%Find the convex hull of the FOV and plot it
%reFOV = real(fovAhat); imFOV = imag(fovAhat);
%cvh = convhull(reFOV,imFOV);

%fovAhat = reFOV(cvh)+1i*imFOV(cvh);
minfov  = abs(min(fovAhat));



