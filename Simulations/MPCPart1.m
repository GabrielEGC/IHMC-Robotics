%clc, clear all, close all
% PROBLEM SETUP, DO NOT CHANGE
kc=39.478;
bc=6.283;
T=10^-2;
vmax=200;
kp=300;kd=20;
Km=-[kp kd];
P =[   0.004485757051635   0.000294369065356
   0.000294369065356   0.000019317722835];
P=[0.243839640250665   0.008321491867835
   0.008321491867835   0.000345913223515];
Apre=[0 1 0;-kc -bc 1; 0 0 0];
Bpre=[0;0;1];
A_sys = eye(3)+T*Apre;
B_sys = T*Bpre;
%N = 100;
%x0 = [-1;0;0];%0.2*randn(3,1);
% QP SETUP HERE
Qm=[P, [0;0];zeros(1,3)]/max(max(P))+[-Km 1]'*[-Km 1]/max(max([-Km 1]'*[-Km 1]));
H = 2*[kron(eye(N+1),Qm), zeros(3*N+3,N);zeros(3*N+3,N)', zeros(N)];
%%
f = zeros(4*N+3,1);
A = [zeros(N,3*N+3) eye(N);zeros(N,3*N+3) -eye(N)];
b = [vmax*ones(N,1);vmax*ones(N,1)];

B1=[kron(eye(N),A_sys) zeros(3*N,3) kron(eye(N),B_sys)];
BMPC1=[eye(3*N+3) zeros(3*N+3,N)]-[zeros(3,4*N+3);B1];
%BMPC2=[zeros(4,4*N) eye(4) zeros(4,2*N)];
%c1 = [x0;zeros(3*N,1)];
%c2 = zeros(4,1);

%
B = [BMPC1];
%c = [c1];

% SOLVE QP
options = optimset('Display','Off');
%     [z,fval,exitflag] = quadprog(H,f,A,b,B,c,[],[],[],options);
% u=z(3*N+4:end);
% u=reshape(u,[1 N]);
% %%
% x(:,1)=x0;
% for i=1:N
%     x(:,i+1)=A_sys*x(:,i)+B_sys*u(:,i);
% end
% plot(T*(1:size(x,2)),x')
% legend('x','dx','F')
% figure(2)
% plot(T*(1:size(x,2)),x'*[-Km 1]')
% legend('F-Km*x')
% figure(3)
% plot(T*(1:size(u,2)),u)
% legend('Fdot')