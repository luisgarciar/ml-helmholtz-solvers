function [M] = mass(np,bc)
%% MASS: Constructs a mass matrix for the 1D Helmholtz model problem.
%  Constructs the mass matrx corresponding to the discretization of a 
%  Helmholtz problem  with P1 finite elements in (0,1) with mixed or 
%  Sommerfeld boundary conditions. The model problem is given by
%
%  Mixed bc:
%  - u''- k^2*u-i*eps*u = f in (0,1)
%  u(0)        = 0   
%  u'(1)-ku(1) = 0
%
%  Sommerfeld bc:
%  - u''- k^2*u-i*eps*u = f in (0,1)
%    u'(0) + ku(0) = 0   
%    u'(1) - ku(1) = 0
% 
%  - u''- k^2*u-i*eps*u = f in (0,1)
%    u(0)        = 0   
%    u'(1)-ku(1) = 0
% 
% Use: M  = mass(np,bc)
%      np:  number of interior points
%      bc: 'mix' or 'som'
%
% Output:
%   M = mass matrix
%   size(M) = (np+2,np+2) if bc = 'som'  
%   size(M) = (np+1,np+1) if bc = 'mix' 
%
% References:
% F. Ihlenburg and I. Babuska, Finite element solution of the 
% Helmholtz equation with high wave number Part I: The h-version of the FEM
%
% S Langdon , S Chandler-Wilde, Finite element methods for acoustic scattering 
% 
%  Author: Luis Garcia Ramos, 
%          Institut fur Mathematik, TU Berlin
%  Version 0.2 - Jun 2017
%

%%
h = 1/(np+1);   %gridsize
switch bc
    case 'mix' %mixed boundary conditions
        e = ones(np+1,1);
        M = h*spdiags([(1/6)*e (2/3)*e (1/6)*e],[-1 0 1],np+1,np+1);
        M(np+1,np+1)=h/3; %Mass matrix
        
    case 'som' %Sommerfeld boundary conditions
        e = ones(np+2,1);
        M = h*spdiags([(1/6)*e (2/3)*e (1/6)*e],[-1 0 1],np+2,np+2);
        M(1,1) = h/3;
        M(np+2,np+2) = h/3; %Mass matrix    
    
end

end



