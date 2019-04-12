function GS(A,b,npoints)

T_Old = zeros(npoints,1); % Initialize unknown vector to zeros
T = zeros(npoints,1);
counter = 0; % Counts number of iterations

% Don't use vectorized code - we need to mimic C
while true
    
    counter = counter+1;
    
    for i=1:npoints
        sum1=0;
        sum2=0;
        for j=1:npoints
            if i~=j && i<j
                sum1 = sum1 + A(i,j)*T_Old(j,1);
            end
            if i~=j && i>j
                sum2 = sum2 + A(i,j)*T(j,1);
            end
        end
        T(i,1) = (b(i,1) - sum1 - sum2)/A(i,i);
    end
    
    diff = abs(T(:,1) - T_Old(:,1));
    T_Old = T;
    
    if max(diff)<1e-7
        break;
    end
    
end

fprintf('No. of iterations in GS: %d \n', counter);
disp('GS done');

end