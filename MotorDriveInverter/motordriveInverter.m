% Parameters of the motor drive inverter system
R=1;
L=0.1;
K=0.01;
J=0.05;
B=0.1;
Ke=0.01;
T_load=0.1;
s=tf('s');
P_motor = K/((J*s+B)*(L*s+R)+K^2);

A=[0 1 0;0 -B/J K/J;0 -Ke/L R/L];
B =[0; 0; 1/L];
C = eye(3);
D = zeros(3,1);
motor_ss = ss(A,B,C,D);
sys = ss(A,B,C,D);
t =0:0.01:10;
u = ones(size(t));
[y,t,x] = lsim(sys,u,t);

subplot(2,1,1);
plot(t,y(:,1),'b',t,y(:,2),'r',t,y(:,3),'g');
xlabel('Time(s)');
ylabel('Motor variables');
legend({'Position','Velocity','Current'});
title('Motor drive inverter system response');

