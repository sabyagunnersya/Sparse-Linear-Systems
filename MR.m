function MR(A,b,npoints,nx,ny)

T_Old = zeros(npoints,1); % Initialize unknown vector to zeros
T = zeros(npoints,1);
prodInit = zeros(npoints,1);
p = zeros(npoints,1);
counter = 0; % Counts number of iterations

for i=1:npoints
    sum = 0;
    for j=1:npoints
        sum = sum + A(i,j)*T(j,1);
    end
    prodInit(i,1) = sum;
end

res = b(:,1) - prodInit(:,1);

for i=1:npoints
    sum = 0;
    for j=1:npoints
        sum = sum + A(i,j)*res(j,1);
    end
    p(i,1) = sum;
end


while true
    
    counter = counter +1;
    
    num=0; denom=0;
    for i=1:npoints
        num = num + p(i,1)*res(i,1);
    end
    for i=1:npoints
        denom = denom + p(i,1)*p(i,1);
    end
    
    alpha = num/denom;
    
    T(:,1) = T(:,1) + alpha*res(:,1);
    
    res(:,1) = res(:,1) - alpha*p(:,1);

    % Computing p for next iteration
    for i=1:npoints
        sum = 0;
        for j=1:npoints
            sum = sum + A(i,j)*res(j,1);
        end
        p(i,1) = sum;
    end
    
    
    diff = abs(T(:,1) - T_Old(:,1));
    T_Old = T;
    
    if max(diff)<1e-7
        break;
    end

end

fprintf('No. of iterations in MR: %d \n', counter);
disp('MR done');

% Printing the final temperature field
% for i=1:ny
%     for j=1:nx
%         index = (i-1)*nx + j;
%         fprintf('%f ',T(index,1));
%     end
%     fprintf('\n');
% end

end

