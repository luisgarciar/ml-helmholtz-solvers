function cubeIdx = getCubeIdx(p, box, h, ni, nj)

x = p(:, 1) - box(1);
y = p(:, 2) - box(3);
z = p(:, 3) - box(5);

i = floor(x/h)+1;
j = floor(y/h)+1;
k = floor(z/h)+1;

nij = (ni-1)*(nj-1);

cubeIdx = (j-1)*nij + (i-1)*(ni-1) + k;


