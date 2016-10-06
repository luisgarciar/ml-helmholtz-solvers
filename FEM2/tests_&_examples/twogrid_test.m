%% Test of a two-grid cycle for the Helmholtz problem 
%  with  small k (otherwise diverges)
%  

dim  = 1;
k    = 1;
eps  = 0.5*k^2;                          %Imaginary part of shift (for shifted Laplacian)
npf  = 100; npc = npf/2;                 %number of points in fine and coarse grids
h    = 1/npf; grid = h*(1:1:npf); 
f    = ones(npf,1); f(npf,1)=0.5; f=h*f; %right hand side

%Helmholtz matrix, restriction, interpolation, Galerkin coarse grid matrix
A = helmholtzfem(k,npf,0);       
restr  = fwrestrictionfem(npf,1);
interp = lininterpolfem(npc,1);
Ac = restr*A*interp;   

x = zeros(npf,1);
L = sparse(tril(A,-1)); U = sparse(triu(A,1)); D=sparse(diag(diag(A)));
N = -(L+U);

error = zeros(numit,1);
ex_sol = A\f;

for i= 1:numit
    %Presmoothing with damped Jacobi
    for j=1:npre 
      x = w*(D\(N*x)+D\f)+(1-w)*x;
    end
    
    res = f-A*x;
    
    %Restrict residual to coarse grid and solve coarse error equation exactly
    rc   = restr*res; 
    ec   = Ac\rc;
    
    %interpolate coarse error to fine grid and correct xf
    ef = interp*ec;
    x  = x + ef; 
    
    %Postsmoothing with damped Jacobi
    for j=1:npos 
      x = w*(D\(N*x)+D\f)+(1-w)*x;
    end
    
    error(i) = norm(ex_sol-x);
end

semilogy(error)






