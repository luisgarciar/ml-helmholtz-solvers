%% AMG TEST V: DIFFERENT INTERPOLATION OPERATORS
% 
% We consider the effect of using different interpolation operators to
% interpolate fine grid values by coarse grid values. It is tested through
% the following choices in option.interpolation
%
% * 's' standard interpolation. Use the matrix A_fc as a weighted average
% of all connected coarse nodes.
% * 't' two-points interpolation. Use at most two connected coarse
% nodes.
% * 'a' aggegration (one-point) interpolation. Use the strongest connected
% coarse node.

%%
clear all; close all;
%% Unstructured mesh in 2-D
load lakemesh
showmesh(node,elem);

%% Classical coarsening and two-points interpolation
option.coarsen = 'c';
option.interpolation = 't';
[N,itStep,time,err] = amgtest(node,elem,[],option);
%% 
display('Table: Classical AMG')
colname = {'Size','Step','Time','Error'}; 
disptable(colname, N,[],itStep,[],time,'%4.2f',err,'%0.5e');

%% Smoothed aggregation
option.coarsen = 'a';
option.interpolation = 'sa';
[N,itStep,time,err] = amgtest(node,elem,[],option);
%% 
display('Table: Smoothed-aggregation AMG')
colname = {'Size','Step','Time','Error'}; 
disptable(colname, N,[],itStep,[],time,'%4.2f',err,'%0.5e');

%% Unstructured mesh in 3-D
clear all; close all
load bunny;
showboundary3(node,elem);
view([-179 74]);

%% Classical coarsening and two-points interpolation
option.coarsen = 'c';
option.interpolation = 't';
[N,itStep,time,err] = amgtest3(node,elem,1,option);
%% 
display('Table: Classical AMG')
colname = {'Size','Step','Time','Error'}; 
disptable(colname, N,[],itStep,[],time,'%4.2f',err,'%0.5e');

%% Smoothed aggregation
option.coarsen = 'a';
option.interpolation = 'sa';
[N,itStep,time,err] = amgtest3(node,elem,[],option);
%% 
display('Table: Smoothed-aggregation AMG')
colname = {'Size','Step','Time','Error'}; 
disptable(colname, N,[],itStep,[],time,'%4.2f',err,'%0.5e');
