%test of real and complex givens rotations
%
v1 = rand(2,1);
[c1,s1,r1] = givens(v1(1),v1(2));

G1 = [c1, s1; -s1, c1];
G1*v1;
r1;

v2 = rand(2,1) + 1i*rand(2,1);
[c2,s2,r2] = givens(v2(1),v2(2));
G2 = [c2, s2; -s2', c2];

A     = rand(5,5); %+ 1i*rand(5,5);
[Q,R] = qrgivens(A);