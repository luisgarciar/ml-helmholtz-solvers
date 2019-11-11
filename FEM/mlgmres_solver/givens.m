function [c,s,r] = givens(f,g)
%Computes the angles for a Givens rotation of a 2x1 complex vector
%
% Input:  (a,b)
% Output: [c,s]

%Reference: E. Anderson, Discontinuous Plane Rotations and the Symmetric
%                        Eigenvalue Problem
%%  

% if f and g are both real
if (isreal(f) && isreal(g))
    if g==0 && f~=0
        c = sign(f);
        s = 0;
        r = abs(f);
    elseif f==0
        c = 0;
        s = sign(g);
    elseif (abs(f)<abs(g))
        t = g/f;
        u = sign(f)*sqrt(1+t^2);
        c = 1/u;
        s = t*c;
        r = f*u;
    else
        t=f/g;
        u=sign(g)*sqrt(1+t^2);
        s=1/u;
        c=t*s;
        r=g*u;
    end
    
else
    if g == 0 && f~=0
        c = sign(f);
        s = 0;
        r = f*c;
    elseif f == 0
        c = 0;
        s = sign(g');
        r = abs(g);
    else
        f1 = abs(real(f)) + abs(imag(f));
        g1 = abs(real(g)) + abs(imag(g));
        
        if (f1 >= g1)
            fs = f/f1;
            f2 = real(fs)^2+imag(fs)^2;
            gs = g/f1;
            g2 = real(gs)^2 + imag(gs)^2;
            u  = sign(real(f))*sqrt(1+g2/f2);
            c = 1/u;
            s = gs'*fs*(c/f2);
            r = f*u;
        else
            fs = f/g1;
            f2 = real(fs)^2+imag(fs)^2;
            gs = g/g1;
            g2 = real(gs)^2+imag(fs)^2;
            u  = sign(real(f))*g1*sqrt(f2+g2); 
            f1 = abs(f);
            fs = f/f1;
            c  = f1/u;
            s  = fs*(g'/u);
            r =  fs*u;
        end
    end    
end
    
end


 

