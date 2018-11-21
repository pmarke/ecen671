%% 22
%%
%   a) Find the condition numbers of A and A^TA
A = [10000 10001;
     10001 10002;
     10002 10003;
     10003 10004;
     10004 10005];
 
cond(A)
cond(A'*A)

%%
% b) Compute the least squares solution using the function x_hat =
% (A'A)^-1A^Tb
b = [20001,20003,20005,20007,20009]';
x_hat = inv(A'*A)*A'*b

%%
% c) Compute the solution using QR decomposition
v1 = A(:,1) + sign(A(1,1))*norm(A(:,1))*[1,0,0,0,0]'
Pv1 = v1*v1'/(v1'*v1)
Q = eye(5) - 2*Pv1
R = H1*A
y = Q \ b
x = R \ y

%%
% d) Compute the solution using the Cholesky factorization
