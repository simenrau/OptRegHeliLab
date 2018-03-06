function [c, c_eq] = Nonlincon1(z)
N = 100;

c = [];
c_eq = 0.15*sin(4*z(1:6:N*6)) - z(5:6:N*6);
end 