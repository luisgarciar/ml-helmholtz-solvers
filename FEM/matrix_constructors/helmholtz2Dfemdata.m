function pde = helmholtz2Dfemdata(k)
%% HELMHOLTZ2DFEMDATA trigonometric data for Helmholtz equation
%
%    -div(grad u)-k^2 u = f;
%
%    with absorbing boundary condition
%
%    grad(u) dot n - i*ku = 0; 
%    u = sin(k*pi*x)^2*sin(k*pi*y)^2;
%    f = -div(grad u)- k^2 u;
%
% Created by Jie Zhou.
%
% Copyright (C)  Long Chen. See COPYRIGHT.txt for details.

pde = struct('f',@f,'exactu',@exactu,'k2',@k2,'Du',@Du);

    C = 10; % amplification factor
    % load data (right hand side function)
    function rhs =  f(p)
        global k
        x = p(:,1); y = p(:,2);
        rhs = -2*pi^2*k^2*(cos(2*k*pi*x).*sin(k*pi*y).^2+cos(2*k*pi*y).*sin(k*pi*x).^2)-k^2*sin(k*pi*x).^2.*sin(k*pi*y).^2;
        rhs = C*rhs;
    end
    % exact solution
    function u =  exactu(p)
        global k
        x = p(:,1); y = p(:,2);
        u = sin(k*pi*x).^2.*sin(k*pi*y).^2;
        u = C*u;
    end
    % Derivative of the exact solution
    function uprime =  Du(p)
        global k
        x = p(:,1); y = p(:,2);
        uprime(:,1) = k*pi*sin(2*k*pi*x).*sin(k*pi*y).^2;
        uprime(:,2) = k*pi*sin(2*k*pi*y).*sin(k*pi*x).^2;
        uprime = C*uprime;
    end
    function wavenumber = k2(~)
        global k
        wavenumber = k.^2;
    end
end