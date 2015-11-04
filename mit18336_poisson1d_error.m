function err2= mit18336_poisson1d_error
%MIT18336_POISSON1D_ERROR
%    Sets up and solves sequence of 1d Poisson problems on
%    domain [0 1], with homogeneous Dirichlet b.c.
%    Plots the convergence of the error in log-log scale.

% 02/2009 by Benjamin Seibold
%            http://www-math.mit.edu/~seibold/
% Feel free to modify for teaching and learning.

n = ceil(10.^linspace(1,4,20));                % sequence of cell numbers
h = 1./(n+1);                            % corresponding grid resolutions
err2 = n*0;
clf
for i = 1:length(n)
    x = (h(i):h(i):1-h(i))';       % regular grid without boundary points
                                                      % 1D Poisson matrix
    K1D = spdiags(ones(n(i),1)*[-1 2 -1],-1:1,n(i),n(i))/h(i).^2;
    f = f_corr(x);                                      % right hand side
    u = K1D\f;                                      % solve linear system
    u_diff = u-u_corr(x);                          % calculate difference
    u_diff = u_diff/norm(u_diff,inf);
    subplot(1,2,1), plot(x,u_corr(x),'b-',x,u,'r-',x,u_diff,'k-')
    legend('true solution','numerical approximation','scaled difference')
    title('1d poisson equation')
    subplot(1,2,2)
    err2(i) = norm(u-u_corr(x),inf);  % calculate error in max norm
    figure(111)
    loglog(10.^[-4 0 -4],10.^[-4 0 -8],'k-',h,err2,'b.-')
    title('error convergence')
    drawnow
end

%========================================================================

function y = u_corr(x)                                 % correct solution
c = 9*pi;
y = sin(c*x.^2);

function y = f_corr(x)                                  % right hand side
c = 9*pi;
y = sin(c*x.^2).*(2*c*x).^2-cos(c*x.^2)*(2*c);