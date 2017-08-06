function pde = helmholtz2Dfemdata
%% HELMHOLTZ2DFEMDATA trigonometric data for Helmholtz equation
%
% Returns the struct pde with data related to the
% the Helmholtz problem
%   
%     -div(grad u)-k^2 u = f;
%    with absorbing boundary condition
%
%    grad(u) dot n - i*ku = g; 
%    u = sin(k*pi*x)^2*sin(k*pi*y)^2;
%    f = -div(grad u)- k^2 u;
%
%   Usage:  Set the wavenumber k as a global variable in 
%           the main .m file
%   Output: 
%       pde: struct containing the following data:
%         All function handles to be applied to input of size (N,2) 
%           'f'     : function handle for right hand side                     
%           exactu': function handle for exact solution
%           'gradu':
%           'k2': squared wavenumber
%            'g': function handle for boundary data
%
% Created by Jie Zhou.
%
% Modified by Luis Garcia Ramos, Aug 2017.
%
% Copyright (C)  Long Chen. See COPYRIGHT.txt for details.

pde = struct('f',@f,'exactu',@exactu,'k2',@k2,...,
             'g',@g,'gradu', @gradu);

    %load data (right hand side function)
    function rhs =  f(p)
        global k
        x = p(:,1); y = p(:,2);
        rhs = -2*pi^2*k^2*(cos(2*k*pi*x).*sin(k*pi*y).^2+cos(2*k*pi*y).*sin(k*pi*x).^2)-k^2*sin(k*pi*x).^2.*sin(k*pi*y).^2;
        rhs = rhs;
    end

    % exact solution
    function u =  exactu(p)
        global k
        x = p(:,1); y = p(:,2);
        u = sin(k*pi*x).^2.*sin(k*pi*y).^2;
    end

    % gradient of the exact solution
    function grad_u =  gradu(p)
        global k
        x = p(:,1); y = p(:,2);
        grad_u(:,1) = k*pi*sin(2*k*pi*x).*sin(k*pi*y).^2;
        grad_u(:,2) = k*pi*sin(2*k*pi*y).*sin(k*pi*x).^2;
    end

    function wavenumber = k2(~)
        global k
        wavenumber = k.^2;
    end
    function g_B = g(p)
        global k
        g_B = p; %?? fix this
    end

end