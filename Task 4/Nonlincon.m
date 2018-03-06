function [c, c_eq] = Nonlincon(z)
alpha = 0.2;
beta = 20;
lambda_t = (2*pi)/3;
q1 = 1; 
q2 = 1; 
N=40;
c = zeros(N,1);
c = alpha*exp(-beta*(z(1:6:N*6)-lambda_t).^2)-z(5:6:N*6);

c_eq = [];
end 

