function pde = helmholtz2Dconstantwndata(k,factoreps,poweps)
%% HELMHOLTZ2CONSTANTWNDATA 
%  Data for Helmholtz/shifted Laplace problem with 
%  constant wavenumber
%
% Returns a struct with data related to the
% the Helmholtz/shifted Laplace problem
%   
%  -div(grad u)-(k^2 + i*eps)u = 0 in Omega= (0,1)x(0,1)
%  grad(u) dot n - i*ku = g in bd(Omega) 
%
%
%  Usage:  
%  Input: wavenumber k, factoreps, poweps 
%       (for shifted Laplacian eps= factoreps*k^poweps)
%
%   Output: 
%       pde:    struct containing the following data:
%       All function handles to be applied to input of size (N,2)
%            'k':  wavenumber
%               
%
% Created by Jie Zhou.
%
% Modified by Luis Garcia Ramos, Aug 2017.
%
% Copyright (C)  Long Chen. See COPYRIGHT.txt for details.


kk  = k; 
pde = struct('f',@f,'exactu',@exactu,'k',kk,...,
             'g',@g,'gradu', @gradu);

    %load data (right hand side function)
    function rhs =  f(p)
        %x = p(:,1); %y = p(:,2);
        rhs = zeros(length(p),1);
    end

    % exact solution
    function u =  exactu(p)
        c = cos(t); s = sin(t);
        x = p(:,1); y = p(:,2);
        arg = k*x*c + k*y*s;
        u = cos(arg) + 1i*sin(arg);
    end

    %gradient of the exact solution
    function grad_u =  gradu(p) 
        c = cos(t); s = sin(t);
        x = p(:,1); y = p(:,2);
        arg = k*x*c + k*y*s;
        grad_u(:,1) = -c*k*sin(arg) + 1i*c*k*cos(arg); %du/dx
        grad_u(:,2) = -s*k*sin(arg) + 1i*s*k*cos(arg); %du/dy
    end

    function gB = g(p)
        x = p(:,1); y=p(:,2);
        gB = zeros(size(x));
        
        %boundaries
        south = (y<=x)&(y<=(1-x)); 
        west  = (y<=(1-x))&(y>x);  
        north = (y>(1-x))&(y>=x);
        east  = (y<x) & (y>(1-x));
        
        %compute g using boundary functions
        gB(south) = g_south(x(south),k,t);
        gB(west)  = g_west(y(west),k,t);
        gB(north) = g_north(x(north),k,t);
        gB(east)  = g_east(y(east),k,t);
    end

end