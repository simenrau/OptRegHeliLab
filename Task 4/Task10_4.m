%remember to run 'init.m' prior to this
close all;
clc;

%% Continuous time state space
A_c = [0 1  0         0         0        0
       0 0 -K_2       0         0        0
       0 0  0         1         0        0
       0 0 -K_1*K_pp -K_1*K_pd  0        0
       0 0  0         0         0        1
       0 0  0         0         -K_3*K_ep -K_3*K_ed];
  
B_c = [0        0
       0        0
       0        0
       K_1*K_pp 0
       0        0
       0        K_3*K_ep];
  
%% Discretized using forward Euler
Ts = 0.25;                      %sampling time
A_d = eye(6) + Ts*A_c;          
B_d = Ts*B_c; 

%% Constants 
alpha = 0.2;
beta = 20;
lambda_t = 2*pi/3;
q1 = 1; 
q2 = 1; 

%% QP-problem
x0 = [pi;0;0;0;0;0];
u0 = [0;0];

N = 40;
M = N;

R  = 2*diag([q1,q2]);                        %weighting of q
Q1 = 2*diag([1 0 0 0 0 0]);
Q  = gen_q(Q1,R,N,M);

z = zeros(N*6+2*M,1);
z_0 = z;


A1 = A_d;
B1 = B_d;
mx = size(A1,2);
mu = size(B1,2);

A_eq = gen_aeq(A1,B1,N,mx,mu);
B_eq = zeros(size(A_eq,1),1);
B_eq(1:mx) = A1*x0;
        
% Bounds
xl(1:mx,1) = -Inf*(ones(mx,1)); 
xu(1:mx,1) = Inf*(ones(mx,1));
xl(3) = -30*pi/180;
xu(3) = 30*pi/180;



ul = [-30*pi/180;-inf];
uu = [30*pi/180;inf];

% Constraints
[vlb,vub] = gen_constraints(N,M,xl,xu,ul,uu);
vlb(320) = 0;
vub(320) = 0;

%% Solve QP problem with linear model
tic
fun = @(z) 0.5*z'*Q*z; 
options = optimoptions('fmincon','Algorithm', 'active-set', 'MaxFunEvals',600000);

[Z, ZVAL, EXITFLAG] = fmincon(fun,z_0,[],[],A_eq,B_eq,vlb,vub,@Nonlincon,options); % hint: quadprog. Type 'doc quadprog' for more info 
toc;

% Calculate objective value
phi1 = 0.0;
PhiOut = zeros(N*mx+M*mu,1);
for i=1:N*mx+M*mu
  phi1=phi1+Q(i,i)*Z(i)*Z(i);
  PhiOut(i) = phi1;
end

%% Extract control inputs and states
u_star1  = [Z(N*mx+1:mu:320);Z(320-1)]; % Control input from solution
u_star2  = [Z(N*mx+2:mu:320);Z(320)];

x1 = [x0(1);Z(1:mx:N*mx)];              % State x1 from solution
x2 = [x0(2);Z(2:mx:N*mx)];              % State x2 from solution
x3 = [x0(3);Z(3:mx:N*mx)];              % State x3 from solution
x4 = [x0(4);Z(4:mx:N*mx)];              % State x4 from solution
x5 = [x0(5);Z(5:mx:N*mx)];              % State x5 from solution
x6 = [x0(6);Z(6:mx:N*mx)];              % State x6 from solution
num_variables = 5/Ts;
zero_padding = zeros(num_variables,1);
unit_padding  = ones(num_variables,1);

u_star1  = [zero_padding; u_star1; zero_padding];
u_star2  = [zero_padding; u_star2; zero_padding];
x1  = [pi*unit_padding; x1; zero_padding];
x2  = [zero_padding; x2; zero_padding];
x3  = [zero_padding; x3; zero_padding];
x4  = [zero_padding; x4; zero_padding];
x5  = [zero_padding; x5; zero_padding];
x6  = [zero_padding; x6; zero_padding];


[K,S,e] = dlqr(A_d,B_d,Q1,R,[]);

K_trans = K';

x_star = [x1 x2 x3 x4 x5 x6]';

t = 0:Ts:Ts*(length(u_star1)-1);
ws_x_star =[[t]',x_star'];
ws = [t',u_star1,u_star2];


%% Plotting
simulation = simout_3;
u_est = simulation(:,1);
u_est2 = simulation(:,2);
travel = simulation(:,3);
travel_rate = simulation(:,4);
pitch = simulation(:,5);
pitch_rate = simulation(:,6);
elevation = simulation(:,7);
elevation_rate = simulation(:,8);
t2 = 0:20/(length(pitch)-1):20;

figure
subplot(811)
<<<<<<< HEAD
stairs(t,u_star1*(180/pi),'k'); ; hold on; plot(t2,u_est*(180/pi)); hold on; plot(t,zeros(length(t),1)+30,'--r');
hold on; plot(t,zeros(length(t),1)-30,'--b'); grid
ylabel('u_1^*'); legend('u_1','u_1^*');
subplot(812)
stairs(t,u_star2*(180/pi),'k'); ; hold on; plot(t2,u_est2*(180/pi)); hold off; grid
ylabel('u_2^*'); legend('u_2','u_2^*');
=======
stairs(t,u_star1*(180/pi),'r'); ; hold on; plot(t2,u_est*(180/pi)); hold off; grid
ylabel('u_1^*')
subplot(812)
stairs(t,u_star2*(180/pi),'r'); ; hold on; plot(t2,u_est2*(180/pi)); hold off; grid
ylabel('u_2^*')
>>>>>>> 29b455799cc63d52fe3a4b40e62298412dcf4d14
subplot(813)
plot(t,x1*(180/pi),'ko'); hold on; plot(t2,travel); hold off; grid
ylabel('travel'); legend('Estimated travel','Measured travel','Location','SouthWest');
subplot(814)
plot(t,x2*(180/pi),'ko'); hold on; plot(t2,travel_rate); hold off; grid
ylabel('travel rate'); legend('Estimated travel rate','Measured travel rate','Location','Northeast');
subplot(815)
plot(t,x3*(180/pi),'ko'); hold on; plot(t2,pitch); hold off; grid
ylabel('pitch'); legend('Estimated pitch','Measured pitch');
subplot(816)
plot(t,x4*(180/pi),'ko'); hold on; plot(t2,pitch_rate); hold off; grid
ylabel('pitch rate'); legend('Estimated pitch rate','Measured pitch rate');
subplot(817)
<<<<<<< HEAD
plot(t,x5*(180/pi),'ko'); hold on; 
=======
plot(t,x5*(180/pi),'mo'); hold on; 
>>>>>>> 29b455799cc63d52fe3a4b40e62298412dcf4d14
plot(t2,elevation); hold off; grid
ylabel('elevation'); legend('Estimated elevation','Measured elevation');
subplot(818)
plot(t,x6*(180/pi),'ko'); hold on; plot(t2,elevation_rate); hold off; grid
ylabel('elevation rate'); legend('Estimated elevation rate','Measured elevation rate');
xlabel('time [s]')



