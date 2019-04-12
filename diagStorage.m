% Domain is 0<=x<=1 and 0<=y<=1 ; This has already been normalized
% Left and bottom boundary conditions are T=0; Right and Top boundary conditions are T=1
% Perform the analysis for delX = delY = 0.1 and then plug in 0.05 and 0.01 for the other cases

delX = 0.05;
delY = 0.05;
nx = 1/delX + 1; % number of grid points along x dimension
ny = nx; % we are using a square grid
npoints = nx*ny;

rhs = zeros(npoints,1); % contains the values on RHS

% filling in the matrix for the lower boundary points
for i=1:nx % does not include T_11 for case 1: delX=delY=0.
    rhs(i,1)=0;
end

% filling in the matrix for the left boundary points
for i=1:nx:((ny-1)*nx) % does not include T_111 for case 1: delX=delY=0.
    rhs(i,1)=0;
end

% filling in the matrix for the top boundary points
for i=(nx*ny):-1:(nx*ny)-nx % includes T_111 for case 1: delX=delY=0.1
    rhs(i,1)=1;
end

% filling in the matrix for the right boundary points
for i=nx:nx:(nx*ny) % includes T_11 for case 1: delX=delY=0.1
    rhs(i,1)=1;
end

c=1;

size = (nx-2)*(ny-2);
A = zeros(size,size);
b = zeros(size,1);

for i=2:(nx-1)
    for j=2:(ny-1)
        
        index = (i-1)*ny + j;
        A(c,c)=4;
        
        if j~=(ny-1)
            A(c,c+1)=-1;
        else
            b(c,1) = b(c,1) + rhs(index+1,1);
        end
        
        if j~=2
            A(c,c-1)=-1;
        else
            b(c,1) = b(c,1) + rhs(index-1,1);
        end
        
        if i~=2
            A(c,c-(ny-2))=-1;
        else
            b(c,1) = b(c,1) + rhs(index-ny,1);
        end
        
        if i~=(nx-1)
            A(c,c+(ny-2))=-1;
        else
            b(c,1) = b(c,1) + rhs(index+ny,1);
        end
        
        c = c+1;
        
    end
end

% The co-efficient matrix is complete now
save ('coeffDiag.dat','A','-ASCII')
disp('Coeff matrix done');


SD(A, b, size, (nx-2), (ny-2));
MR(A, b, size, (nx-2), (ny-2));
CG(A, b, size, (nx-2), (ny-2));

