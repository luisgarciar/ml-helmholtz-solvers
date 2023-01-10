function [Q,R] = qrgivens(A)
  [m,n] = size(A);
  Q = eye(m);
  R = A;

  for j = 1:n
    for i = m:-1:(j+1)
      G = speye(m);
      [c,s,r] = givens(R(i-1,j),R(i,j));
      G([i-1, i],[i-1,i]) = [c  s; -s' c];
      R = G*R; R(i-1,j) = r; R(i,j) = 0;
      Q = Q*G';
    end
  end

end