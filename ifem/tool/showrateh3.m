function showrateh3(h1,err1,k1,opt1,str1,h2,err2,k2,opt2,str2,h3,err3,k3,opt3,str3)
%% SHOWRATEH3 rate of two error sequences
%
% showrate3(N1,err1,k1,opt1,str1,N2,err2,k2,opt2,str2,N3,err3,k3,opt3,str3)
% plots the err1 vs N1, err2 vs N2, and err3 vs N3 in the loglog scale.
% Additional input
%
%   - k1, k2, k3: specify the starting indices; see showrate
%   - opt1, opt2, opt3: the line color and style 
%   - str1, str2, str3: strings used in legend
%
% Example
%
% showrate2(N,energyErr,1,'r-+','||u-u_h||_A',...
%           N,L2Err,1,'b-+','||u-u_h||');
%
% See also showrate, showresult, showmesh, showsolution
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.

N1 = 1./h1; N2 = 1./h2; N3 =1./h3;
if (nargin<=2) 
    k1 = 1; opt1 = '-*';
end
r1 = showrate(N1,err1,k1,opt1);
hold on
r2 = showrate(N2,err2,k2,opt2);
r3 = showrate(N3,err3,k3,opt3);
title(['Rate of convergence is Ch^{' num2str(-round(r3)) '}'],'FontSize', 14);
h_legend = legend(str1,['C_1h^{' num2str(-round(r1,2)) '}'],...
                  str2,['C_2h^{' num2str(-round(r2,2)) '}'],...
                  str3,['C_3h^{' num2str(-round(r3,2)) '}'],'LOCATION','Best');
set(h_legend,'FontSize',12);
xlabel('log(1/h)');
