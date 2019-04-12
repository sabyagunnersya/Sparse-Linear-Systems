function BICGSTAB(A,b,npoints,nx,ny)

T = NaN(npoints,1);
for i=1:npoints
    T(i,1) = 0.01;
end

T_Old = T; % Storing values from previous iteration to check if the values have converged or not
prod0 = zeros(npoints,1); % prod0 represents Ax0
prod = zeros(npoints,1); % prod represents Ap
prod_S = zeros(npoints,1);

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

for i=1:npoints
    sum = 0;
    for j=1:npoints
        sum = sum + A(i,j)*p(j,1);
    end
    prod(i,1) = sum; % Calculating Ap
end



while true
    
    counter = counter +1;
    
    num=0; denom=0;
    for i=1:npoints
        num = num + r(i,1)*rStar(i,1);
        denom = denom + prod(i,1)*rStar(i,1);
    end    
    alpha = num/denom; % Calculating alpha
    denom = num; % storing for beta calculation
    
    s(:,1) = r(:,1) - alpha*prod(:,1); % Calculating s
    
    for i=1:npoints
        sum1 = 0;
        for j=1:npoints
            sum1 = sum1 + A(i,j)*s(j,1);
        end
        prod_S(i,1) = sum1; % Calculating As
    end
    
    omegaNum=0; omegaDenom=0;
    for i=1:npoints
        omegaNum = omegaNum + prod_S(i,1)*s(i,1);
        omegaDenom = omegaDenom + prod_S(i,1)*prod_S(i,1);
    end    
    omega = omegaNum/omegaDenom; % Calculating omega
    
    T(:,1) = T(:,1) + alpha*p(:,1) + omega*s(:,1);
    
    r(:,1) = s(:,1) - omega*prod_S(:,1);

    num=0; 
    for i=1:npoints
        num = num + r(i,1)*rStar(i,1); % updated r
    end   
    beta = (num*alpha)/(denom*omega);
    
    p(:,1) = r(:,1) + beta*(p(:,1)-omega*prod(:,1));
    
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

fprintf('No. of iterations in BICGSTAB: %d \n', counter);
disp('BICGSTAB done');

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

save ('finalDisplacementFieldBICGSTAB.dat','T','-ASCII')

end

    