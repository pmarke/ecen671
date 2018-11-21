%% Problem 3

% Generate 10 random square complex matrices

num_matrices = 10;;
matrix_dim = 3
M = cell(num_matrices,1);
for k  = 1:num_matrices
cMatrix = rand(matrix_dim,matrix_dim) + i*rand(matrix_dim,matrix_dim);
M{k} = cMatrix;
end


%% Compute QR factorization using my method and matlabs method
myQ = cell(num_matrices,1);
myR = cell(num_matrices,1);
bifQ = cell(num_matrices,1);
bifR = cell(num_matrices,1);


for k = 1:num_matrices

[Q,R] = myqr(M{k},1)
myQ{k} = Q;
myR{k} = R;

[Q,R] = qr(M{k})
bifQ{k} = Q;
bifR{k} = R;

end

%% Check that Q is unitary, R is an upper triangle and QR = A
tests = zeros(3*num_matrices,1);
for k = 1:num_matrices
    
   % Unitary test
   tests(k) =  double(round(trace(myQ{k}'*myQ{k}),8) == matrix_dim);
   
   % Upper triangle test
   tests(k+num_matrices) = double(istriu(round(myR{k},8))) ;
   
   % QR = A test
   A = myQ{k}*myR{k};
   tests(k+2*num_matrices) = double(round(norm(M{k}-A),8) == 0);
    
    
end

% This number shold be 3*num_matrices if it pasts the previous tests.
sum(tests)

%% Compare myqr to the Matlab function qr. Do you get the same results? Why or why not? is QR factorization unique?
%
% By inspection I can see that my myqr implementation yields different 
% results than Matlab's built in function. This supports the fact that 
% QR factorization is not unique. This is understandable since you are
% essentially describing 1 mxn matrix by 2 matrices that are mxm and mxn.
% other words, since QR has more parameters than A, you have an onto
% mapping from QR to A. This is easily proven by looking at the first element
% of A. $a_{11} = q_{11}*r_{11}$. Thus the first element of A is a combination
% of the first elements of q and r. 


