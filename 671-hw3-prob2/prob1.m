%% Problem 1

%% a:
%     Write a script that numerically generates the first 5 Legendre polynomials by using Gram-Schmidt
%     orthogonalization of the set of vectors {1,t,t^2,t^3,t^4} over the
%     interval [-1,1]. Plot each of these polynomials.

num = 10000;
t = linspace(-0.99, .99, num);
p = zeros(num, 5);
for kk = 1:5
 p(:, kk) = t.^(kk-1);
end 

e = zeros(num,5);
q = zeros(num,5);


for i = 1:5
   e(:,i) = p(:,i);
   for j = 1:i-1
       e(:,i) = e(:,i) - dot(p(:,i),q(:,j))*q(:,j);
   end
   q(:,i) = e(:,i)/dot(e(:,i),e(:,i));
end

figure(1),clf;
plot(q);
legend('1','t','t^2','t^3','t^4')

%% b.
%  Perform a linear least-squares approximation of the function f(t)=e^{-t} on the interval [-1,1]
%  using the Legendre polynomials derived in part 1. Plot f(t) and your
%  approximation (or expansion) of f(t). Find the norm of the error vector.
%  That is, find the norm of the difference between f(t) and your
%  approximation of f(t).

f = exp(-t)';

% Perform linear least squares
x = q\f;

f_est = (q*x);

norm_err = sum((f-f_est).^2);

figure(2);
plot(t,f,'b.',t,f_est','r');
legend('f','f\_{est}');
title(["f vs f_{est} with norm error",num2str(norm_err)])

%% c. 
%     Modify your results form (a) to generate the first 5 Chebyshev polynomials on the interval [-1,1].
%     refer to example 2.15.1 on p.120 of the book, and note that
%     generation of the Chebyschev polynomials will be identical to the
%     Legendre polynomials but with the use of a weighted inner product.

w = 1./sqrt(1-t.^2)';

e_c = zeros(num,5);
q_c = zeros(num,5);

for i = 1:5
   e_c(:,i) = p(:,i);
   for j = 1:i-1
       e_c(:,i) = e_c(:,i) - dot(w.*p(:,i),q_c(:,j))*q_c(:,j);
   end
   q_c(:,i) = e_c(:,i)/dot(e_c(:,i),e_c(:,i));
end

figure(3),clf;
plot(q_c);

%% d. 
%  Repeat part (b) using the Chebyshev polynomials instead of the Legendre polynomials.

f_c = exp(-t)';

% Perform linear least squares
x_c = q_c\f_c;

f_c_est = (q_c*x_c);

norm_err_c = sum((f_c-f_c_est).^2);

figure(4);
plot(t,f_c,'b.',t,f_c_est','r');
legend('f','f\_{est}');
title(["f vs f_{est} with norm error",num2str(norm_err_c)])

%% e. Discuss results
% Since both span the same space, they can approximate it equally as well. 