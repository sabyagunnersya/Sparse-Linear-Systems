% Domain is 0<=x<=1 and 0<=y<=1 ; This has already been normalized
% Left and bottom boundary conditions are T=0; Right and Top boundary conditions are T=1
% Perform the analysis for delX = delY = 0.1 and then plug in 0.05 and 0.01 for the other cases

delX = 0.05;
delY = 0.05;
nx = 1/delX + 1; % number of grid points along x dimension
ny = nx; % we are using a square grid
npoints = nx*ny;

A = zeros(npoints,npoints); % coefficient matrix
b = zeros(npoints,1); % contains the values on RHS

% filling in the matrix for the lower boundary points
for i=1:nx % does not include T_11 for case 1: delX=delY=0.
    b(i,1)=0;
    A(i,i)=1;
end

% filling in the matrix for the left boundary points
for i=1:nx:((ny-1)*nx) % does not include T_111 for case 1: delX=delY=0.
    b(i,1)=0;
    A(i,i)=1;
end

% filling in the matrix for the top boundary points
for i=(nx*ny):-1:(nx*ny)-nx % includes T_111 for case 1: delX=delY=0.1
    b(i,1)=1;
    A(i,i)=1;
end

% filling in the matrix for the right boundary points
for i=nx:nx:(nx*ny) % includes T_11 for case 1: delX=delY=0.1
    b(i,1)=1;
    A(i,i)=1;
end

index = 0;
% filling in the interior points
for i=2:(nx-1)
    for j=2:(ny-1)
        index = (j-1)*nx + i;
        A(index,index)=-4;
        A(index,index+1)=1;
        A(index,index-1)=1;
        A(index,index-nx)=1;
        A(index,index+nx)=1;
    end
end

% The co-efficient matrix is complete now
save ('coeff.dat','A','-ASCII')
disp('Coeff matrix done');
omegaOpt = creatingIterationMatrix(A, npoints);
SORiterations(A, b, npoints, omegaOpt);
GS(A, b, npoints);
jacobi(A, b, npoints);

SD(A, b, npoints, nx, ny);
% MR(A, b, npoints, nx, ny);
CG(A, b, npoints, nx, ny);

