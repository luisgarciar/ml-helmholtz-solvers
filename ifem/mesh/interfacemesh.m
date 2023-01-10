function [node,tElem,sElem,bdEdge, interfaceNodeIdx] = interfacemesh(phi,box,hmax,varargin)
%% INTERFACEMESH generated interface fitted mesh by Delaunay algorithm on Cartesian grid
%
%
%  Huayi Wei <weihuayi@xtu.edu.cn>
%

pbox = [-0.25 -0.25 5.25 5.25];
ff = 'eps';

% construct the initial structure mesh
[node,elem,T] = squarequadmesh(box,hmax);
edge = T.edge;
edge2elem = T.edge2elem;
bdEdge = T.bdEdge;
clear T;
N = size(node,1);
NT = size(elem,1);

figure('Color',[1 1 1]);
showmesh(node,elem);
hold on;
plotcircle(0.5);
set(gcf, 'PaperPosition', pbox); 
set(gcf, 'PaperSize', [5 5]); 
saveas(gcf, 'circle1', ff)

% compute the phi value
phiValue = phi(node,varargin{:});
phiValue(abs(phiValue) < hmax^2) = 0;
vSign = msign(phiValue);
% find the edge  intersected  by interface
isCutEdge = phiValue(edge(:,1)).* phiValue(edge(:,2))<0;

A = node(edge(isCutEdge,1),:);
B = node(edge(isCutEdge,2),:);

% construct the intersected points
cutNode = findintersectbisect(phi,A,B);
NC = size(cutNode, 1);
vSign(N+1:N+NC) = 0;

% find interface elem 
isInterfaceElem = false(NT,1);  
isInterfaceElem(edge2elem(isCutEdge,[1,2])) = true;
isInterfaceElem(sum(abs(vSign(elem)), 2) < 3) = true;

% find aux points
isSpecialElem = sum(vSign(elem), 2) == 0 & sum(abs(vSign(elem)), 2) == 2;
selem = elem(isSpecialElem,:);
isSpecialElem = vSign(selem(:, 1)).*vSign(selem(:, 3)) == -1 | vSign(selem(:, 2)).*vSign(selem(:, 4)) == -1;
selem = selem(isSpecialElem,:);
auxPoint = (node(selem(:,1),:) + node(selem(:, 3),:))/2.0;
NA = size(auxPoint,1);
vSign(N+NC+1:N+NC+NA) = 0;

figure('Color',[1 1 1]);
showmesh(node,elem);
findnode(cutNode,'all', 'noindex', 'color','r','MarkerSize',18);
findnode(auxPoint, 'all', 'noindex', 'color','m', 'MarkerSize', 18);


% find interface node
isInterfaceNode = false(N,1);
isInterfaceNode(elem(isInterfaceElem,:)) = true;

findnode(node, isInterfaceNode, 'noindex','color','k','MarkerSize',18);
set(gcf, 'PaperPosition', pbox); 
set(gcf, 'PaperSize', [5 5]); 
saveas(gcf, 'circle2', ff)

% add intersected points to node array and interfaceNode array, respectively 
node = [node;cutNode;auxPoint];
interfaceNode = [node(isInterfaceNode,:);cutNode;auxPoint];

% construct index mapping
interfaceNodeIdx = [find(isInterfaceNode);N + (1:NC)'; N + NC + (1:NA)'];

% construct the Delaunay triangulation of interfaceNode
dt = DelaunayTri(interfaceNode(:,1),interfaceNode(:,2));
tElem = dt.Triangulation;
tElem = fixorder(interfaceNode,tElem); % correct the order of elem node indexs

figure('Color',[1 1 1]);
showmesh(interfaceNode,tElem);
set(gcf, 'PaperPosition', pbox); 
set(gcf, 'PaperSize', [5 5]); 
saveas(gcf, 'circle3', ff)

% get rid of the unnecessary triangle elem
NI = sum(isInterfaceNode);
isUnnecessaryElem = sum(tElem<=NI,2) == 3;
tElem = tElem(~isUnnecessaryElem,:);

figure('Color',[1 1 1]);
showmesh(interfaceNode,tElem);
set(gcf, 'PaperPosition', pbox); 
set(gcf, 'PaperSize', [5 5]); 
saveas(gcf, 'circle4', ff)

% map interfaceNode index to node index
tElem = interfaceNodeIdx(tElem);

% get the remainder quad elems
sElem = elem(~isInterfaceElem,:);

figure('Color',[1 1 1]);
showmesh(node,tElem);
hold on 
showmesh(node, sElem)
findnode(node, vSign == 0, 'noindex','color','r','MarkerSize',18)
set(gcf, 'PaperPosition', pbox); 
set(gcf, 'PaperSize', [5 5]); 
saveas(gcf, 'circle5', ff)

end
