%% test of V-cycle
sigma=0;
lmax=10;
itmax=5;
convhist=zeros(itmax,lmax-4);

for l=5:lmax
        nl=2?l; nt=nl+1; x=0:nl; x=x?/nl;
        
        v=sin(6*pi*x)+sin(17*pi*x);
        f=zeros(nt,1);
        convhist(it,l-4)=norm(vt,inf)/norm(v,inf);
        v=vt;

for it = 1:itmax
        vt=twogrid(v,f,2,2);
end

end
%
% print the data from the table
%
convhist