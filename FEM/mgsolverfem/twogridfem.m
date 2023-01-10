function [x_sol,relres] = twogridfem(A,restr,interp,f,x_init,npre,npos,w,smo,numcycles)
%Construct Galerkin coarse matrix
Ac        =  restr*A*interp;
x_sol     =  x_init;
relres    =  zeros(numcycles+1,1);
relres(1) =  norm(f-A*x_sol);

for i=1:numcycles
    %Presmoothing
    x_sol = smoother(A,f,x_sol,w,nu);
    res     = f-A*x_sol; 

    %Restrict to coarse grid and solve coarse error equation exactly
    rc   = restr*res; 
    ec   = Ac\rc;
    
    %interpolate coarse error to fine grid and correct xf
    ef     = interp*ec;
    x_sol  = x_sol + ef; 
    
    %Post-smoothing
    x_sol   = smoother(A,f,x_sol,mu,w);
    relres(i+1) = norm(f-A*x_sol);

end %end of cycle
relres = relres/relres(1); %relative residuals
end %end of function
