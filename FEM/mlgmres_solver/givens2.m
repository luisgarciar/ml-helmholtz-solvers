function [c,s,r] = givens2(a,b)
%Computes the angles for a Givens rotation of a 2x1 complex vector
%
% Input:  (a,b)
%
% Output: [c,s] such that [c  -s][a] = [r]
%                         [s'  c][b]   [0]
%
%
%References: Golub, Van Loan, Matrix Computations, 4th Ed.
%            E. Anderson, Discontinuous Plane Rotations and the Symmetric
%                        Eigenvalue Problem
%            
%% 
% if f and g are both real

if (isreal(a) && isreal(b))
    if b==0 
        c = 1;
        s = 0;
        r = abs(b);
    elseif (abs(b)>abs(a))
        t = -a/b;
        u = sign(a)*sqrt(1+t^2);
        c = 1/u;
        s = t*c;
        r = a*u;
    else
        t=a/b;
        u=sign(b)*sqrt(1+t^2);
        s=1/u;
        c=t*s;
        r=b*u;
    end
else
    if abs(a)==0
        c = 0;
        s = 1;
        r = b;
    else
        t = abs(a)+abs(b);
        n = t*sqrt(abs(a)/t^2 + abs(b)/t^2);
        alpha = sign(a);
        c = abs(a)/n; s = alpha*b'/n;
        r =  sign(a)*sqrt(abs(a)^2+abs(b)^2);
    end
    %See:
    %Matrix Computations, Golub, Van Loan 4th Ed., Chap5. P244
    %u1 = real(a); u2 = imag(a);
    %v1 = real(b); v2 = imag(b);
    
    %[ca,sa,ru] = givens2(u1,u2);
    %[cb,sb,rv] = givens2(v1,v2);
    %[ct,st,r]  = givens2(ru,rv);
    %t = (ca*cb+sa*sb)+1i*(ca*sb-sb*sa); 
    %c = ct; s = st*t;  
    
end
    
end


 

