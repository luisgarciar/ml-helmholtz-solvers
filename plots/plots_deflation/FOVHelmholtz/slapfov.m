function [fovSL,min]=slapfov(A,Aeps)
%% SLAPFOV
%  [fovSL,min]= SLapfov(A,Aeps) computes the field of values of the
%  matrix S=A*inv(Aeps), where A is the discrete Helmholtz operator 
%  and Aeps is the shifted Laplacian with shift eps
%   
%  INPUT
%  A: Helmholtz sparse matrix
%  Aeps: shifted Laplace matrix

%  OUTPUT
%   
%  m: Inner numerical radius of A*inv(Aeps) (distance of the fov to zero)
%  fv, eigvA: Field of values and eigenvalues of A
% 
% Author: Luis Garcia Ramos, TU Berlin
%         Based on a routine of Nick Higham (Matrix Toolbox)
          % version 0.1 - Apr 2017
         
%%%