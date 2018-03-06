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
t2 = 0:35/(length(pitch)-1):35;

figure(2)
subplot(811)
stairs(t,u_star1*(180/pi),'k','LineWidth',2); hold on; plot(t2,u_est*(180/pi),'Color',[0 114/255 189/255],'LineWidth',2); hold on; 
plot(t,zeros(length(t),1)+30, '--r'); hold on; plot(t,zeros(length(t),1)-30,'--b'); hold off; grid
ylabel('u_1^*')
subplot(812)
stairs(t,u_star2*(180/pi),'k','LineWidth',2); hold on; plot(t2,u_est2*(180/pi),'Color',[0 114/255 189/255],'LineWidth',2); hold off; grid
ylabel('u_2^*')
subplot(813)
plot(t,x1*(180/pi),'ko'); hold on; plot(t2,travel,'r','LineWidth',2); hold off; grid
ylabel('travel'); legend('Estimated travel','Measured travel');
subplot(814)
plot(t,x2*(180/pi),'ko'); hold on; plot(t2,travel_rate,'r','LineWidth',2); hold off; grid
ylabel('travel rate'); legend('Estimated travel rate','Measured travel rate','Location','Northeast');
subplot(815)
plot(t,x3*(180/pi),'ko'); hold on; plot(t2,pitch,'r','LineWidth',2); hold off; grid
ylabel('pitch'); legend('Estimated pitch','Measured pitch');
subplot(816)
plot(t,x4*(180/pi),'ko'); hold on; plot(t2,pitch_rate,'r','LineWidth',2); hold off; grid
ylabel('pitch rate'); legend('Estimated pitch rate','Measured pitch rate');
subplot(817)
plot(t,x5*(180/pi),'ko'); hold on; plot(t2,elevation,'r','LineWidth',2); hold off; grid
ylabel('elevation'); legend('Estimated elevation','Measured elevation');
axis([0 35 -20 20]);
subplot(818)
plot(t,x6*(180/pi),'ko'); hold on; plot(t2,elevation_rate,'r','LineWidth',2); hold off; grid
ylabel('elevation rate'); legend('Estimated elevation rate','Measured elevation rate');
xlabel('time [s]')