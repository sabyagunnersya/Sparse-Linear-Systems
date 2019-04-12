function CG(A,b,npoints,nx,ny)

T_Old = zeros(npoints,1); % Initialize unknown vector to zeros
T = zeros(npoints,1);
prodInit = zeros(npoints,1);
prod = zeros(npoints,1);
temp = zeros(npoints,1);
counter = 0; % Counts number of iterations

for i=1:npoints
    sum = 0;
    for j=1:npoints
        sum = sum + A(i,j)*T(j,1);
    end
    prodInit(i,1) = sum;
end

res = b(:,1) - prodInit(:,1);
p = res;

for i=1:npoints
    sum = 0;
    for j=1:npoints
        sum = sum + A(i,j)*p(j,1);
    end
    prod(i,1) = sum;
end



while true
    
    counter = counter +1;
    
    num=0; denom=0;
    for i=1:npoints
        num = num + res(i,1)*res(i,1);
    end
    for i=1:npoints
        denom = denom + prod(i,1)*res(i,1);
    end    
    alpha = num/denom;
    denom = num; % storing for beta calculation
    
    T(:,1) = T(:,1) + alpha*p(:,1);
    
    res(:,1) = res(:,1) - alpha*prod(:,1);

    num=0; 
    for i=1:npoints
        num = num + res(i,1)*res(i,1); % updated res
    end   
    beta = num/denom;
    
    temp(:,1) = res(:,1) + beta*p(:,1);
    p=temp;
    
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

fprintf('No. of iterations in CG: %d \n', counter);
disp('CG done');

% Printing the final temperature field
% for i=1:ny
%     for j=1:nx
%         index = (i-1)*nx + j;
%         fprintf('%f ',T(index,1));
%     end
%     fprintf('\n');
% end

end

    