diag = [0 5 1 -2 2;1 7 -3 2 1;2 7 4 0 1;-2 6 1 -2 1;1 4 1 -1 1;3 -5 1 -2 0;0 -7 1 3 0;-1 2 0 1 0;1 -3 2 0 0;2 -2 0 0 0];
b = zeros(10,1);
b = [3;-5;15;-3;5;4;-1;-1;-6;4];
offset = [-1 0 1 2 4];

d(1,1) = 5;
d(1,2) = 1;
d(1,3) = -2;
d(1,5) = 2;
for i=2:10
    for j=1:5
        d(i,i+offset(j)) = diag(i,j);
    end
end

size = 10;
A = zeros(size,size);

for i=1:10
    for j=1:10
        A(i,j)=d(i,j);
    end
end


% omegaOpt = creatingIterationMatrix(A, size);

GS(A, b, size);
jacobi(A, b, size);
SORiterations(A, b, size, 1.05);
