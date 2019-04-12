function SORiterations(A, b, npoints, w)

T_Old = zeros(npoints,1); % Initialize unknown vector to zeros
T = zeros(npoints,1);
TStar = zeros(npoints,1);

counter = 0; % Counts number of iterations

% Don't use vectorized code - we need to mimic C
while true
    
%     Points 1, 4 and 34 are the left boundary points with displacement = 0
%     So, set T(2n-1,1) = 0 and T(2n,1) = 0 for each of the points above
%     T(1,1) = 0;
%     T(2,1) = 0;
%     T(7,1) = 0;
%     T(8,1) = 0;
%     T(67,1) = 0;
%     T(68,1) = 0;        
    
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
        TStar(i,1) = (b(i,1) - sum1 - sum2)/A(i,i);
        T(i,1) = T_Old(i,1) + w*(TStar(i,1)-T_Old(i,1));
    end
    
    diff = abs(T(:,1) - T_Old(:,1));
    T_Old = T;
    
    if max(diff)<1e-7
        break;
    end
    
end

fprintf('No. of iterations in SOR: %d \n', counter);
disp('SOR iterations done');

end