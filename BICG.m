function BICG(A,b,npoints,nx,ny)

T_Old = zeros(npoints,1); % Storing values from previous iteration to check if the values have converged or not
T = zeros(npoints,1); % Initialize unknown vector to zeros
prod0 = zeros(npoints,1); % prod0 represents Ax0

prod = zeros(npoints,1); % prod represents Ap
prodTranspose = zeros(npoints,1);

counter = 0; % Counts number of iterations

% Calculating Ax0
for i=1:npoints
    sum = 0;
    for j=1:npoints
        sum = sum + A(i,j)*T(j,1);
    end
    prod0(i,1) = sum; 
end

r = b(:,1) - prod0(:,1);
rStar = r; % rand(npoints,1); % r0 times r0_Star must not be 0
p = r;
pStar = rStar;

for i=1:npoints
    sum = 0;
    for j=1:npoints
        sum = sum + A(i,j)*p(j,1);
    end
    prod(i,1) = sum; % Calculating Ap
end



while true

%     Points 1, 4 and 34 are the left boundary points with displacement = 0
%     So, set T(2n-1,1) = 0 and T(2n,1) = 0 for each of the points above
%     T(1,1) = 0;
%     T(2,1) = 0;
%     T(7,1) = 0;
%     T(8,1) = 0;
%     T(67,1) = 0;
%     T(68,1) = 0;        
    
    counter = counter +1;
    
    num=0; denom=0;
    for i=1:npoints
        num = num + r(i,1)*rStar(i,1);
        denom = denom + prod(i,1)*pStar(i,1);
    end    
    alpha = num/denom;
    denom = num; % storing for beta calculation
    
    T(:,1) = T(:,1) + alpha*p(:,1);
    
    r(:,1) = r(:,1) - alpha*prod(:,1);
    
    B = transpose(A);
    
    for i=1:npoints
        sum1 = 0;
        for j=1:npoints
            sum1 = sum1 + B(i,j)*pStar(j,1);
        end
        prodTranspose(i,1) = sum1; % Calculating A.T dot pStar
    end

    rStar(:,1) = rStar(:,1) - alpha*prodTranspose(:,1);

    num=0; 
    for i=1:npoints
        num = num + r(i,1)*rStar(i,1); % updated res
    end   
    beta = num/denom;
    
    p(:,1) = r(:,1) + beta*p(:,1);
    pStar(:,1) = rStar(:,1) + beta*pStar(:,1);
    
    % Computing prod for next iteration
    for i=1:npoints
        sum = 0;
        for j=1:npoints
            sum = sum + A(i,j)*p(j,1);
        end
        prod(i,1) = sum;
    end
    
    
    diff = abs(T(:,1) - T_Old(:,1));
    T_Old = T;
    
    if max(diff)<1e-7
        break;
    end

end

fprintf('No. of iterations in BICG: %d \n', counter);
disp('BICG done');

% % Printing the final temperature field
% for i=1:ny
%     for j=1:nx
%         index = (i-1)*nx + j;
%         fprintf('%f ',T(index,1));
%     end
%     fprintf('\n');
% end

% Printing the final displacement field
for i=1:npoints
    fprintf('%f ', T(i,1));
end
fprintf('\n');

save ('finalDisplacementFieldBICG.dat','T','-ASCII')

end

    