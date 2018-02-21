%remember to run 'init.m' prior to this
close all;
clc;

%% Continuous time state space
A_c = [0 1  0         0
       0 0 -K_2       0
       0 0  0         1
       0 0 -K_1*K_pp -K_1*K_pd];
  
B_c = [0;0;0;K_1*K_pp];
  
%% Discretized using forward Euler
Ts = 0.25;                      %sampling time
A_d = eye(4) + Ts*A_c;          
B_d = Ts*B_c; 
   
%% QP-problem
x0 = [pi;0;0;0];

N = 100;
M = N;

R = 2*1;                        %weighting of q

Q1 = 2*diag([1 0 0 0]);
Q = gen_q(Q1,R,N,M);

f = zeros(500,1);
z = f;

A1 = A_d;
B1 = B_d;
mx = size(A1,2);
mu = size(B1,2);

A_eq = gen_aeq(A1,B1,N,mx,mu);
B_eq = zeros(size(A_eq,1),1);
B_eq(1:mx) = A1*x0;
        
% Bounds
xl = -Inf*(ones(mx,1)); 
xu = Inf*(ones(mx,1));
xl(3) = -30*pi/180;
xu(3) = 30*pi/180;
ul = -30*pi/180;
uu = 30*pi/180;

% Constraints
[vlb,vub] = gen_constraints(N,M,xl,xu,ul,uu);
vlb(500) = 0;
vub(500) = 0;

%% Solve QP problem with linear model
tic
[z,lambda] = quadprog(Q,f,[],[],A_eq,B_eq,vlb,vub,x0); % hint: quadprog. Type 'doc quadprog' for more info 
t1=toc;

% Calculate objective value
phi1 = 0.0;
PhiOut = zeros(N*mx+M*mu,1);
for i=1:N*mx+M*mu
  phi1=phi1+Q(i,i)*z(i)*z(i);
  PhiOut(i) = phi1;
end

%% Extract control inputs and states
u_star  = [z(N*mx+1:N*mx+M*mu);z(N*mx+M*mu)]; % Control input from solution

x1 = [x0(1);z(1:mx:N*mx)];              % State x1 from solution
x2 = [x0(2);z(2:mx:N*mx)];              % State x2 from solution
x3 = [x0(3);z(3:mx:N*mx)];              % State x3 from solution
x4 = [x0(4);z(4:mx:N*mx)];              % State x4 from solution

num_variables = 5/Ts;
zero_padding = zeros(num_variables,1);
unit_padding  = ones(num_variables,1);

u_star   = [zero_padding; u_star; zero_padding];
x1  = [pi*unit_padding; x1; zero_padding];
x2  = [zero_padding; x2; zero_padding];
x3  = [zero_padding; x3; zero_padding];
x4  = [zero_padding; x4; zero_padding];


[K,S,e] = dlqr(A_d,B_d,Q1,R,[]);

K_trans = K';

x_star = [x1 x2 x3 x4]';

t = 0:Ts:Ts*(length(u_star)-1);
ws_x_star =[[t]',x_star'];
ws = [t',u_star];


%% Plotting
simulation = simout_3;
u_est = simulation(:,1);
travel = simulation(:,2);
travel_rate = simulation(:,3);
pitch = simulation(:,4);
pitch_rate = simulation(:,5);
t2 = 0:35/(length(pitch)-1):35;

figure
subplot(511)
stairs(t,u*(180/pi)); ; hold on; plot(t2,u_est*(180/pi)); hold off; grid
ylabel('u')
subplot(512)
plot(t,x1*(180/pi),'mo'); hold on; plot(t2,travel); hold off; grid
ylabel('travel'); legend('Estimated travel','Measured travel');
subplot(513)
plot(t,x2*(180/pi),'mo'); hold on; plot(t2,travel_rate); hold off; grid
ylabel('travel rate'); legend('Estimated travel rate','Measured travel rate','Location','Southeast');
subplot(514)
plot(t,x3*(180/pi),'mo'); hold on; plot(t2,pitch); hold off; grid
ylabel('pitch'); legend('Estimated pitch','Measured pitch');
subplot(515)
plot(t,x4*(180/pi),'mo'); hold on; plot(t2,pitch_rate); hold off; grid
ylabel('pitch rate'); legend('Estimated pitch rate','Measured pitch rate');
xlabel('time [s]')



