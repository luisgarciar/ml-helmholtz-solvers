
clear all
close all

%% Parameters for two grid cycle:
 sigma =10;
 n  = 2^6;
 x  = [0:1/n:1]'; %1-D grid
 f  = -2*pi*exp(x).*cos(pi*x)+(pi^2+sigma-1)*exp(x).*sin(pi*x); %rhs function
 nu = 3;
 mu = 4;

 %% Call Function
[ v, sol, error ] = twogridcycle( f, sigma, nu, mu ) 
 
%% Plot Approximation
plot(x,sol);
hold on
plot(x,v,'r*');


%% 
% xc=[0:2/n:1];%coarse grid  
% fc=fwrestriction(f);
% 
% plot(x,f,'go');
% hold on
% plot(xc,fc,'r*');
% ff=lininterpol(fc);
% plot(x,ff,'b');