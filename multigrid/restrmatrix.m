n=2^3;
I=eye(n-1);
E=zeros(n/2-1,n-1);
EH=zeros(n-1,n/2-1);
Id=eye(n/2-1);

for i=1:length(I)
    E(:,i)  =  fwrestriction(I(:,i));
end

for j=1:length(Id)
    EH(:,j) = lininterpol(Id(:,j));
end