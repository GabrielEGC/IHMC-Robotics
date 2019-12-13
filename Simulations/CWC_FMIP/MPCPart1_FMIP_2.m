clc, clear all, close all
% PROBLEM SETUP, DO NOT CHANGE

Xfo=0.07;
Yfo=0.1;


T=10^-2;
g=9.81;gv=[0;0;-g];
m=150;

M03=zeros(3,3);
Mgcr=cross(repmat(gv,1,3),eye(3));
Apre=[M03 m*Mgcr M03;M03 M03 eye(3);M03 M03 M03];

Bpre=[M03 eye(3); M03 M03; 1/m*eye(3) M03];
nX=size(Apre,2);
nU=size(Bpre,2);
%
A_sys = eye(nX)+T*Apre;
B_sys = T*Bpre;
N = 300;

L0 = [0;0;0];
r0 = [-.4;-.6;1.1];
dr0= [2;2;-3];

L20 = L0-m*cross(r0,dr0);


rfd = [0;0.05;1.1];

X0 = [L20;r0-rfd;dr0];


mu=0.5;
mup=0.5/sqrt(2);
Wcwc=[];
Wred=[1 0 -mup 0 0 0
      -1 0 -mup 0 0 0
      0 1 -mup 0 0 0
      0 -1 -mup 0 0 0
      0 0 -1 0 0 0
      0 0 -Yfo 1 0 0 
      0 0 -Yfo -1 0 0
      0 0 -Xfo 0 1 0
      0 0 -Xfo 0 -1 0
      Yfo Xfo -mup*(Xfo+Yfo) mup mup 1
      Yfo -Xfo -mup*(Xfo+Yfo) +mup -mup 1
      -Yfo Xfo -mup*(Xfo+Yfo) mup -mup 1
      -Yfo -Xfo -mup*(Xfo+Yfo) -mup -mup 1  
      Yfo Xfo -mup*(Xfo+Yfo) -mup -mup -1
      Yfo -Xfo -mup*(Xfo+Yfo) -mup +mup -1
      -Yfo Xfo -mup*(Xfo+Yfo) mup -mup -1
      -Yfo -Xfo -mup*(Xfo+Yfo) mup mup -1];
nWred=size(Wred,1);
bcwcred=zeros(size(Wred,1),1);
b3cwcred=bcwcred+Wred*([m*gv;m*Mgcr*rfd]);
%
% QP SETUP HERE
%Qm=[P, [0;0];zeros(1,3)]/max(max(P))+[-Km 1]'*[-Km 1]/max(max([-Km 1]'*[-Km 1]));
Qm=diag([0.001 0.001 0.001 1 1 1000 0.01 0.01 0.01]);%diag(kron(ones(1,3),[1 1000 1]));%eye(nX);%
Rm=0.00001*eye(nU);

H = 2*[kron(eye(N+1),Qm), zeros(nX*(N+1),nU*N);zeros(nX*(N+1),nU*N)', kron(eye(N),Rm)];
%%
f = zeros(nX*(N+1)+nU*N,1);

Aqp = [zeros(nWred*N,nX*(N+1)) kron(eye(N),Wred)];

bqp = kron(ones(N,1),b3cwcred);
%%


B1=[kron(eye(N),A_sys) zeros(nX*N,nX) kron(eye(N),B_sys)];
BMPC1=[eye(nX*(N+1)) zeros(nX*(N+1),nU*N)]-[zeros(nX,nX*(N+1)+nU*N);B1];
%BMPC2=[zeros(nX,nX*N) eye(nX) zeros(nX,nU*N)];
c1 = [X0;zeros(nX*N,1)];
%c2 = zeros(9,1);
B = [BMPC1];%;BMPC2];
c = [c1];%;c2];
% SOLVE QP
options = optimset('Display','Off');
    [z,fval,exitflag] = quadprog(H,f,Aqp,bqp,B,c,[],[],[],options);
u=z(nX*(N+1)+1:end);
u=reshape(u,[nU N]);
%%
Wgl=kron(eye(N),Wred);
Hol=Wgl*u(:)<bqp;
find(Hol==0)

x(:,1)=X0;
for i=1:N
    x(:,i+1)=A_sys*x(:,i)+B_sys*u(:,i);
end
%plot(x')
figure(1)
subplot(311)
plot(T*(1:size(x,2)),x(1:3,:)')
legend('L2x','L2y','L2z')
subplot(312)
plot(T*(1:size(x,2)),x(4:6,:)')
legend('x','y','z')
subplot(313)
plot(T*(1:size(x,2)),x(7:9,:)')
legend('$\dot{x}$','$\dot{y}$','$\dot{z}$')
figure(2)
subplot(211)
plot(T*(1:size(u,2)),u(1:3,:)')
legend('$f_x$','$f_y$','$f_z$')
subplot(212)
plot(T*(1:size(u,2)),u(4:6,:)')
legend('$\tau_x$','$\tau_y$','$\tau_z$')
figure(3)
xwo=x+[zeros(3,1);rfd;zeros(3,1)]*ones(1,size(x,2));
plot3(xwo(4,:),xwo(5,:),xwo(6,:)); hold on, grid on
plot3(xwo(4,1),xwo(5,1),xwo(6,1),'xr')
plot3(xwo(4,end),xwo(5,end),xwo(6,end),'ok')
plot3([Xfo Xfo -Xfo -Xfo Xfo],[Yfo -Yfo -Yfo Yfo Yfo],zeros(5,1),'r')
xlabel('x')
ylabel('y')
zlabel('z')