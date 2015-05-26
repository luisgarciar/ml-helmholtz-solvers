function g = rhs(x,y)

h = zeros(size(x)); 
g = zeros(size(y)); 

h(find(x==0.5)) = 1;
g(find(y==0.5)) = 1;
g =  h.*g;

end
