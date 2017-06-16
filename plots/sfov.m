function [fovA,eigvA,m,M]=sfov(A,AH,v0max,N,k)
%% SFOV 
%  [m,M]= sfov(A,AH, k) computes the field of values of the sparse
%  matrix A and approximations to the inner numerical radius and 
%  spectral radius using the method  of Johnson 
%
%  INPUT
%  A:  function handle that computes  A*x
%  AH: function handle that computes  A'*x (Hermitian transpose)
%  N = size(A)
%  k: number of angles

%  OUTPUT
%  m: Inner numerical radius (distance of the fov to zero)
%  M: Spectral radius
%  fv: Field of values

% Author: Luis Garcia Ramos, TU Berlin
%         (Based on a routine of Nick Higham (Matrix Toolbox))
          % version 0.1 - Apr 2017
         
%%%
if nargin == 1, k = 32; end

theta = linspace(0,2*pi,k); %range of angles
fovA  = zeros(k,1);         %boundary points
eigvA = 0;
for j=1:k
    %We rotate the matrix A to obtain At=exp(i*theta(j))*A 
    %and compute the max eigenvalue and unit eigenvector of 
    %Ht = Hermitian part of At
     j
     et = exp(1i*theta(j));  
     Ht = @(x) 0.5*(et*feval(A,x) + et'*feval(AH,x));
    
     opts.isreal = 0;
     opts.v0     = v0max;
     opts.p = 60;
     [vmaxHt,~,flag] = eigs(Ht,N,1,'LR',opts);
     
     %if flag ~=0
      %  fovA(j) = fovA(j-1);
     %else
        %The boundary point is the rotated max eigenvalue
        v  = vmaxHt/norm(vmaxHt);
        fovA(j) = v'*feval(A,v);
        v0max=v;
     %end
  
end

%The approximate fov is the convex hull of the computed boundary points
%clf
%plot(fv,'LineWidth',3);
%hold on

% To compute the inner numerical radius
% we first check if 0 is in fov(A)
xbd = real(fovA);
ybd = imag(fovA);
in  = inpolygon(0,0,xbd,ybd);

 if in==1
     m=0;
 else
     m = min(abs(fovA));
 end
 
M = max(abs(fovA));


end

