clc;
clear all;

A = zeros(16,16);
%s: total number of points in the domain
s = length(A); assert(sqrt(s)==round(sqrt(s))); 
%n: number of points in one dimension
n = sqrt(s);
I = eye(s);

%We construct the permutation matrix for the red-black ordering
if(mod(n,2)==1);
     k=round((s-1)/2);       
     P = sparse(s,s);
     P(:,1:2:(s+1)) = I(:,1:1:(k+1));
     P(:,2:2:(s-1)) = I(:,k+2:1:s);     
else   
    k=round(n/2);
    for i=1:k        
        v1 = ((2*(i-1)*n+1):2:2*(i-1)*n+n)';
        v2 = ((2*i-1)*n+2:2:(2*i-1)*n+n)';
        %red variables
        perm(n*(i-1)+1:n*i,1) = [v1;v2];       
        %black variables
        perm(2*k^2+n*(i-1)+1:2*k^2+n*i,1) = [v1+1;v2-1];       
    end%for    
    P = I(perm,:);      
end%if
     
v=(1:s)';
P*v


