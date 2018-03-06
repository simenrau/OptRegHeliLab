function [c, c_eq] = Nonlincon1(z)
N = 100;

<<<<<<< HEAD
c = [];
c_eq = 0.15*sin(4*z(1:6:N*6)) - z(5:6:N*6);

=======


c = [];
c_eq = 0.2*sin(4*z(1:6:N*6)) - z(5:6:N*6);







% a = -11.5;
% b = 36;
% ce = -28.2;
% 
% c = zeros(N,1);
% c = -a*(z(1:6:N*6)).^2 - b*z(1:6:N*6) - ce + z(5:6:N*6);
% c_eq = [];
>>>>>>> 29b455799cc63d52fe3a4b40e62298412dcf4d14
end 