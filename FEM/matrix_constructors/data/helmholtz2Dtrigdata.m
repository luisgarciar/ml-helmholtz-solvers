function pde = helmholtz2Dtrigdata(k)
%% HELMHOLTZ2DTRIGDATA trigonometric data for Helmholtz equation
%
% Returns the struct pde with data related to the
% the Helmholtz problem
%   
%     -div(grad u)-k^2 u = f;
%    with absorbing boundary condition
%    grad(u) dot n - i*ku = g; 
%
%    with exact solution
%    u = cos(pi*x)*cos(pi*y);
%    f = -div(grad u)- k^2 u;
%
%   Usage:  
%   Input: wavenumber k
%
%   Output: 
%       pde: struct containing the following data:
%            All function handles to be applied to input of size (N,2) 
%           'f'    : function handle for right hand side                     
%           exactu': function handle for exact solution
%           'gradu': function handle for the gradient of exact solution u
%
%           'k2': squared wavenumber
%            'g': function handle for boundary data
%
% Created by Jie Zhou.
%
% Modified by Luis Garcia Ramos, Aug 2017.
%
% Copyright (C)  Long Chen. See COPYRIGHT.txt for details.


k2  = k^2; 
pde = struct('f',@f,'exactu',@exactu,'k2',k2,...,
             'g',@g,'gradu', @gradu);

    %load data (right hand side function)
    function rhs =  f(p)
        x = p(:,1); y = p(:,2);
        rhs = (-k^2+2*pi^2)*(cos(pi*x).*cos(pi*y)); 
    end

    % exact solution
    function u =  exactu(p)
        x = p(:,1); y = p(:,2);
        u = cos(pi*x).*cos(pi*y);
    end

    % gradient of the exact solution
    function grad_u =  gradu(p)
        x = p(:,1); y = p(:,2);
        grad_u(:,1) = -pi*sin(pi*x).*cos(pi*y);
        grad_u(:,2) = -pi*sin(pi*y).*cos(pi*x);
    end

    function g_B = g(p)
        x = p(:,1); y = p(:,2);
        g_B = 1i*k*(cos(pi*x).*cos(pi*y));
    end

end