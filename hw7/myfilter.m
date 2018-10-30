function [e_list,h_est] = myfilter(num_coefficients)

h = [5,4,3,2,1]';

num_coeff = num_coefficients;
h_est = zeros(num_coeff,1);
delta = 1;
P = eye(num_coeff)/delta;
P_dl = P;
h_dl = h_est;

q = zeros(1,5)';
q_est = zeros(1,num_coeff)';


e_list = [];

for i=1:100
    
    % Get new input value
    input = 10*(rand-0.5);
    
    % update q and q_est
    q(1:end-1) = q(2:end);
    q(end) = input;
    
    q_est(1:end-1) = q_est(2:end);
    q_est(end) = input;
    
    % Compute new output
    d = q'*h;
    
   
    k = P*q_est/(1+q_est'*P_dl*q_est);
    e = d -q_est'*h_dl;
    h_est = h_dl + k*e;
    P = P_dl - k*q_est'*P_dl;
    
    % Update past values
    h_dl = h_est;
    P_dl = P;    
    
    e_list = [e_list,e];
    
end

end