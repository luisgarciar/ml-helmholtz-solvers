function [v] = wJacobi(A,v,maxit,f,w)
%WJACOBI Solves a system of equations by the weighted Jacobi method
% 
%   INPUT:
%   A, f: Matrix of the linear system, and right hand side
%   v: initial guess
%   u:
%   w: damping parameter
%   maxit: max number of smoothing steps

[m,n]=size(A);

vold=v;
vnew=zeros(size(v));
vtemp=zeros(size(v));
d=diag(A);

for k=1: maxit
    vtemp(1)=(f(1)-A(1,2:n)*vold(2:n))/d(1);  
    vnew(1)=vtemp(1);
    
    for i=2:m
        vtemp(i)=(f(i)-A(i,1:i-1)*vold(1:i-1)-A(i,i+1:n)*vold(i+1:n))/d(i);
        vnew(i)=vold(i)+w*(vtemp(i)-vold(i));
    end
    
    vold=vnew;
    
end

v=vnew;

end