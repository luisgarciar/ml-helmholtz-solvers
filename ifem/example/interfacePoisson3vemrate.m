clear all; close all;
cube = [-2, 2, -2, 2, -2, 2];
pde = sphereinterfacedata(1, 1000, 1);
option.solver = 'amg';
option.h0 = 0.2;
option.maxIt = 4;

vemInterfacePoisson3(cube, pde, option)
