function x = luSolve(P,L,U,b)

% We are solving the equation
% Ax = b   where PA = LU s.t.
% LUx = P^Tb = c where P^tb = c
% Ly = P^tb  where y = Ux


[r,col] = size(L);

% solve for y
c = P*b;
y = zeros(col,1);

for j=1:col
   
    y(j) = c(j) - L(j,:)*y;
    
end

% solve for x
x = zeros(col,1);

for i=col:-1:1
   
    x(i) = (1/U(i,i))*(y(i) - U(i,:)*x);
    
end





end