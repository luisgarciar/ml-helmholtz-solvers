close all;
clear all;

cube = [-1, 1, -1, 1, -1, 1];
h = 0.025;
surface =spheresurface;

% cube = 2*[-1, 1, -1, 1, -1, 1];
% h = 0.1;
% surface = heartsurface;


% cube = 2*[-1, 1, -1, 1, -1, 1];
% h = 0.05;
% surface = quarticssurface;

% cube = 6*[-1, 1, -1, 1, -1, 1];
% h = 0.1;
% surface = torussurface;


[node,elem, interfaceData] = interfacemesh3(cube, surface.phi, h);

showmesh(node, interfaceData.interface);

flag = checkinterfacemesh3(node,elem, interfaceData)


