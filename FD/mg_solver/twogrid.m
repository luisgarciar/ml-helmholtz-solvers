function [x_sol,relres] = twogrid(A,L,U,D,restr,interp,f,x_init,npre,npos,numcycles)

% L=sparse(tril(A));
% U=sparse(triu(A)); 
% D=sparse(diag(diag(A))); 
M = D+L; N = -U;

%Construct Galerkin coarse matrix
Ac        =  restr*A*interp;
x_sol     = x_init;
relres    = zeros(numcycles+1,1);
relres(1) = norm(f-A*x_sol);

for i=1:numcycles
    %Presmoothing
    for j=1:npre            
        x_sol = M\(N*x_sol+f);
    end   
    
    res     = f-A*x_sol; 
    
    %Restrict to coarse grid and solve coarse error equation exactly
    rc   = restr*res; 
    ec   = Ac\rc;
    
    %interpolate coarse error to fine grid and correct xf
    ef     = interp*ec;
    x_sol  = x_sol + ef; 
    
    %Post-smoothing
    for j=1:npos            
        x_sol = M\(N*x_sol+f);
    end   
    
    relres(i+1) = norm(f-A*x_sol);

end %end of cycle
relres = relres/relres(1); %relative residuals
end %end of function
