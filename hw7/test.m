A = rand(3);
x = rand(3,1)
b = A*x;

% A = [2 4 -5; 6 12.001 1; 4 -8 -3]

[P,L,U] = myLU(A,1,1);


luSolve(P,L,U,b)


%% HW 5.1-4

% without pivoting

A = [2, 4, -5; 
    6, 12.001, 1; 
    4, -8, -3];
b = [-5; -33.002; -21];

[P,L,U] = myLU(A,1,0);

x1 = luSolve(P,L,U,b)

% With pivoting

[L,U,P] = lu(A);

x2 = luSolve(P,L,U,b)