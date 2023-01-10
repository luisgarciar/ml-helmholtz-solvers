%Field of values of Grcar Matrix

G = gallery('grcar',30);
fov(G);
axis equal
set(gca,'FontSize',16);

matlab2tikz('filename','grcarfov.tex','standalone',true, 'extraaxisoptions',['xlabel style={font=\LARGE},'...
                       'ylabel style={font=\LARGE},','ticklabel style={font=\Huge}']);

%set(gca,'Ytick',[-1 -0.5  0 0.5 1],'FontSize',14);


%title('fv(grcar(10))')