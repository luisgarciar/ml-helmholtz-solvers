% Test of sFOV (function that plots the FOV of a matrix)
% We test it with a shifted Laplace preconditioned matrix
clear all
close all

ppw = 12;   %points per wavelength
bc = 'som'; %boundary conditions ('som' or 'dir')

%parameters of the shift eps = factoreps*k^poweps
poweps    = 1;
factoreps = 1;

%Choose which fovs are generated
plot_csl    = 'no';
plot_defcsl = 'no';

kk = [20];

k   = kk(1);
eps = factoreps*(k^poweps);

%if pollution = 'no' the number of points np= ceil(k^(3/2)) 
pollution = 'yes';

%otherwise the number of grid points is chosen
%with a fixed number of points per wavelength
np = ceil(ppw*k/(2*pi))-1; 
if strcmp(pollution,'no')
np = ceil(k^(3/2));
end

if (mod(np+1,2)==1)  %set an even number of interior points in 1D
    np = np+1; 
end
npc = (np-1)/2; 

%computing the FOV of A
A    = helmholtz2_ord1(k,0,np,np,bc);
Aeps = helmholtz2_ord1(k,eps,np,np,bc);
N    = length(A);

Ahatm = full(Aeps\A);

eigvA = eig(Ahatm);

%%

[LH,UH] = lu(Aeps');
[L,U] = lu(Aeps);    

%Compute first the maximum eigenvalue of A
[vmaxA, eigmaxA] =  eigs(A,1,'LR');

Ahat     = @(x) A*(U\(L\x));
AhatH    = @(x) UH\(LH\(A'*x));
H        = @(x) 0.5*(feval(Ahat,x) + feval(AhatH,x)); %Hermitian part of Ahat

opts.isreal = 0;
%opts.p = 20;
opts.v0 = vmaxA;
[vmaxH, eigmaxH] =  eigs(H,N,1,'LR',opts);

[fovAhat,m,M] = sfov(Ahat,AhatH,vmaxH,N,32);


%% Computing the FOV of the deflated matrix
Z = lininterpol(npc,2,bc);
[Q,R] = qr(full(Z));
r = size(Z,2);
Zperp = Q(:,r+1:end);

%Matrix form Ahat=A*inv(Aeps)
[L1,U1] = lu(A); 

Ahperp = Zperp'*Aeps*(U1\(L1\Zperp));
[Lperp,Uperp]   = lu(full(Ahperp));
[LperpH,UperpH] = lu(full(Ahperp'));

Ahperps = @(x)  Uperp\(Lperp\x);
AhperpHs = @(x) UperpH\(LperpH\x);

n = length(Uperp);

%Compute first the maximum eigenvalue of Ahperp
opts2.isreal = 0;
[vmaxAhperp, eigmaxAhperp] =  eigs(Ahperps,n,1,'LR',opts2);

H = @(x) 0.5*(feval(Ahperps,x) + feval(AhperpHs,x)); %Hermitian part of Ahat
[fovAhatperp,m,M] = sfov(Ahperps,AhperpHs,vmaxAhperp,n,32);

