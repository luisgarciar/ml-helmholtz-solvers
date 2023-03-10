function [fovAhat] = computefovAhat(kmult,opts)
%COMPUTEFOVAHAT Computes the field of values of 
%preconditioned Helmholtz matrices

%Input:
%kmult: vector of integers 
        %(wavenumbers, to be multiplied by pi)
%opts: Data structure with options for the computation
%opts.prec = 'csl' or 'adef': preconditioner
%opts.poweps    = power of the shifted Laplacian
%opts.factoreps = factor of the shifted Laplacian
%opts.dim  = dimension of the problem (1 or 2)
%opts.bc  = Boundary conditions, 'som' or 'dir' 
%opts.disc = Discretization options 
%            'q' for quasioptimal --> n = Ck^2
%            'pf' for pollution free --> n = Ck^1.5
%            'pw' for points per wavelength --> n = Ck  
%
%opts.fvpts = number of points in the field of values
%
%
%Output:
%fovAhat = Matrix of size (opts.fvpts,length(kmult)) containing
%          the boundaries of the field of values of the matrices
%          
%Version 0.1 -- Feb 2019
%
%%

dim = opts.dim;
poweps  = opts.poweps;
factoreps = opts.factoreps;
bc = opts.bc;
prec = opts.prec;
disc = opts.disc;
fvpts = opts.fvpts;

%Wavenumber
kk      =  kmult;

switch disc
    case 'q'
        pollexp = 2;
    case 'pf'
        pollexp = 1.5;
    case 'pw'
        pollexp = 1;
end

fovAhat = zeros(fvpts, length(kk));

numCores = feature('numcores');
pool = parpool(numCores);
for i=1:length(kk)
    k   = kk(i);
    eps = factoreps*k^poweps;
    npc = ceil(k^(pollexp)/2);
    npf = 2*npc+1;
   
    A    = helmholtzfem(k,npf,0,bc);           %Helmholtz matrix
    Aeps = helmholtzfem(k,npf,eps,bc);         %Shifted Laplace matrix
    M    = mass(npf,bc);
    R    = fwrestrictionfem(npf,dim,bc);
    Z    = R';             %Prolongation. Deflation subspace: columns of Z
    [L,U] = lu(Aeps);
    Ac = Z'*A*Z; %Coarse operator
    [Lc, Uc] = lu(Ac);
    LcH = Uc'; UcH = Lc';
    

    %Function handles for FOV of shifted Laplacian
    LH    = U'; UH = L';
    N     = length(A);
    
    %Let C = MAepsinv
    C    = @(x) M*(U\(L\x));
    CH   = @(x) UH\(LH\(M*x));    
    
    %Function handles for FOV of Deflated shifted Laplacian      
    %Deflated operator
    N        = length(A);
    Aepsinv  =  @(x) U\(L\x);     %Inverse of shifted Laplace
    AepsHinv =  @(x) UH\(LH\x);   %Inverse of Hermitian transpose of shifted Laplace
    Acinv    =  @(x) Uc\(Lc\x);   %Inverse of coarse Helmholtz
    AcHinv   =  @(x) UcH\(LcH\x); %Inverse of Hermitian transpose of coarse Helmholtz
    
    %Function handles Adef
    AP  =  @(x)  M*Aepsinv(x-A*Z*Acinv((Z'*x)));
    APH =  @(x)  AepsHinv(M*x)- Z*AcHinv(Z'*A'*AepsHinv(M*x));
    v0 = randn(length(A),1);
    
    switch prec
        case 'csl'
            Ahat = C; AhatH = CH;
            
        case 'adef'
            Ahat = AP;  AhatH = APH;
    end
    
    %Field of values of AHat in Euclidean inner product
    [fovAhat(:,i),~,~] = sfov(Ahat,AhatH,v0,N,fvpts);
    fovAhat(:,i) = 1+1i*eps*fovAhat(:,i);
    
   
end
delete(pool)
end