%% AMG TEST I: DIFFERENT MESHES
% 
% We consider linear finite element discretization of the Poisson equation
% with homongenous Dirichlet boundary condition on different meshes. 

%%
clear all; close all;
option.coarsen = 'a';
option.interpolation = 'sa';

%% Uniform mesh
[node,elem] = squaremesh([0,1,0,1],0.1);
[node,elem] = uniformrefine(node,elem);
[node,elem] = uniformrefine(node,elem);
showmesh(node,elem);
snapnow;
[N,itStep,time,err,errHist] = amgtest(node,elem,[],option);
%% 
colname = {'Size','Step','Time','Error'}; 
disptable(colname, N,[],itStep,[],time,'%4.2f',err,'%0.5e');
%%
r = showrate(N,time,2);
xlabel('N'); ylabel('Time');
title(['Complexity is N^{' num2str(r) '}'] ,'Fontsize', 14);

%% Circle mesh
close all;
[node,elem] = circlemesh(0,0,1,0.2);
[node,elem] = uniformrefine(node,elem);
[node,elem] = uniformrefine(node,elem);
showmesh(node,elem);
snapnow;
[N,itStep,time,err] = amgtest(node,elem,[],option);
%% 
colname = {'Size','Step','Time','Error'}; 
disptable(colname, N,[],itStep,[],time,'%4.2f',err,'%0.5e');
%%
r = showrate(N,time,2);
xlabel('N'); ylabel('Time');
title(['Complexity is N^{' num2str(r) '}'] ,'Fontsize', 14);

%% Unstructured mesh
close all;
load lakemesh
showmesh(node,elem);
snapnow;
[N,itStep,time,err] = amgtest(node,elem,[],option);
%% 
colname = {'Size','Step','Time','Error'}; 
disptable(colname, N,[],itStep,[],time,'%4.2f',err,'%0.5e');
%%
r = showrate(N,time,2);
xlabel('N'); ylabel('Time');
title(['Complexity is N^{' num2str(r) '}'],'Fontsize', 14);