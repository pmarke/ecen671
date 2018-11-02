function [P,L,U] = myLU(A, index, pivoting)

% Get the dimensions of the matrix
[r,c] = size(A);

% Find the pivot element
if pivoting == 1
    p_i = index;            % pivot index
    p_e = A(index,index);   % pivot element
    for i = (index+1):r
        if abs(A(i,index)) > abs(p_e)
            p_e = A(i,index);
            p_i = i;
        end
    end

    % Construct the permute matrix
    P_i = eye(r);
    P_i([p_i, index], :) = P_i([index, p_i], :);

else
    P_i = eye(r);
end

% Update U
U_i = P_i*A;

% Construct the elementary matrix
E_i = eye(r);
for i = (index+1):r    
    E_i(i,index) = -U_i(i,index)/U_i(index,index);    
end

% Update U
U_i = E_i*U_i;

% If you are not done traversing the matrix, call this alorithm recursively
if index+1 ~= r
    [P,L,U] = myLU(U_i,index+1,pivoting);
else
    P = eye(r);
    L = eye(r);
    U = U_i;
end

P = P*P_i;
L = inv(P_i)*inv(E_i)*L;

if index ==1
    L = P*L;
end

end