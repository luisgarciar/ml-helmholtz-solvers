function findedge(node,edge,range,varargin)
%% FINDEDGE highlights edges
%
%    FINDEDGE(node,edge,range) finds all elements whose indices are in the
%    range array by displaying these elements in yellow.
%
%    FINDEDGE(node,edge) finds all elements.
%   
%    FINDEDGE(node,edge,'noindex') skip the display of indices.
%
%    FINDEDGE(node,edge,range,'param','value','param','value'...) allows
%    additional patch param/value pairs to highlight the elements.
%    
% Example:
%     [node,elem] = squaremesh([0,1,0,1],1/2);
%     subplot(1,2,1);
%     showmesh(node,elem);
%     T = auxstructure(elem);
%     findedge(node,T.edge,5,'index','color','r','MarkerSize',24);
%     subplot(1,2,2);
%     showmesh(node,elem);
%     findedge(node,T.edge);
%
%   See also findelem3, findnode3.
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.

hold on
% set up range
if (nargin==2) || isempty(range) || (ischar(range) && strcmp(range,'all'))
    range = (1:size(edge,1))'; 
end
if islogical(range)
    range = find(range); 
end
if size(range,2)>size(range,1)
    range = range'; 
end
% draw edges in red
ndim = size(node,2);
midEdge = (node(edge(range,1),:)+node(edge(range,2),:))/2;
if (nargin > 3) && strcmp(varargin{1},'vec') % plot edge vector
    edgeVec = node(edge(range,2),:) - node(edge(range,1),:);    
    if ndim == 2
        h = quiver(midEdge(:,1),midEdge(:,2),edgeVec(:,1),edgeVec(:,2));
    elseif ndim ==3
        h = quiver3(midEdge(:,1),midEdge(:,2),midEdge(:,3), ...
                    edgeVec(:,1),edgeVec(:,2),edgeVec(:,3));
    end
    set(h,'Linewidth',3)
end
if ndim == 2
    h = plot(midEdge(:,1),midEdge(:,2),'r.','MarkerSize', 24);
elseif ndim == 3
    h = plot3(midEdge(:,1),midEdge(:,2),midEdge(:,3),'r.','MarkerSize', 24);
end
if nargin > 3
    if (strcmp(varargin{1},'noindex') || strcmp(varargin{1},'index') ...
            || strcmp(varargin{1},'vec'))
        startidx = 2;
    else
        startidx = 1;
    end
    if  length(varargin) >= startidx
        set(h,varargin{startidx:end});
    end
end
if (nargin <=3) || ~(strcmp(varargin{1},'noindex'))
    if ndim == 2
        text(midEdge(:,1)+0.025,midEdge(:,2)-0.025,int2str(range), ...
             'FontSize',12,'FontWeight','bold','Color','k');
    elseif ndim == 3
        text(midEdge(:,1)+0.025,midEdge(:,2)-0.025,midEdge(:,3),int2str(range), ...
         'FontSize',12,'FontWeight','bold','Color','k');
    end
end