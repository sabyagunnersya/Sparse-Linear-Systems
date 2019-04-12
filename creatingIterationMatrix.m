function omegaOpt = creatingIterationMatrix(A, npoints)

D = zeros(npoints,npoints);
E = zeros(npoints,npoints);
F = zeros(npoints,npoints);

for i=1:npoints
    D(i,i) = A(i,i); % forming the diagonal matrix
end

for i=1:npoints
    for j=1:npoints
        if i>j
            E(i,j) = -A(i,j); % lower triangular matrix
        end
        if i<j
            F(i,j) = -A(i,j); % upper triangular matrix
        end
    end
end

GJacobi = D\(E+F);
GGS = (D-E)\F;
omegaOpt = spectralRadius(GJacobi, GGS);
disp('Iteration matrix done');

end
    
    