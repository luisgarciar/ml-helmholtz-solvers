function p=polyadd(a,b)
%polyadd Addieren von Polynomen


if nargin<2
    error('Not enough Input Argmuents')
end

a=reshape(a,1,[]); %make sure inputs are polynomial row vectors
b=b(:).';           %this makes a row as well

na=length(a);
nb=length(b);

p=[zeros(1,nb-na) a]+[zeros(1,na-nb) b]; %Auffuellen mit Nullen