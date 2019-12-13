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
N = 400;

L0 = [0;0;0];
r0 = [-.4;-.6;1.1];
dr0= [2;2;-1];
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
Qm=diag([0.001 0.001 0.001 1 1 10000 0.01 0.01 0.01]);%diag(kron(ones(1,3),[1 1000 1]));%eye(nX);%
Rm=0.00001*eye(nU);

H = 2*[kron(eye(N+1),Qm), zeros(nX*(N+1),nU*N);zeros(nX*(N+1),nU*N)', kron(eye(N),Rm)];
%%
f = zeros(nX*(N+1)+nU*N,1);

Aqp = [zeros(nWred*N,nX*(N+1)) kron(eye(N),Wred)];

bqp = kron(ones(N,1),b3cwcred);
%%


B1=[kron(eye(N),A_sys) zeros(nX*N,nX) kron(eye(N),B_sys)];
BMPC1=[eye(nX*(N+1)) zeros(nX*(N+1),nU*N)]-[zeros(nX,nX*(N+1)+nU*N);B1];
BMPC2=[zeros(nX,nX*N) eye(nX) zeros(nX,nU*N)];
c1 = [X0;zeros(nX*N,1)];
c2 = zeros(9,1);
B = [BMPC1];%;BMPC2];
c = [c1];%;c2];
% SOLVE QP
options = optimset('Display','Off');
    [z,fval,exitflag] = quadprog(H,f,Aqp,bqp,B,c,[],[],[],options);
u=z(nX*(N+1)+1:end);
u=reshape(u,[nU N]);
xtraj=z(1:nX*(N+1));
xtraj=reshape(xtraj,[nX N+1]);
%%
zeta=u-[m*gv;m*Mgcr*rfd]*ones(1,N);
nv=[0;0;1];
kP=(sum(zeta(1:3,:).*zeta(4:6,:)))./(nv'*zeta(1:3,:));
tau2P=zeta(4:6,:)-nv*kP;
for i=1:N
    lambP(i)=-(cross(zeta(1:3,i),tau2P(1:3,i))/norm(zeta(1:3,i))^2)'*nv/(nv'*zeta(1:3,i));
    rP(:,i)=cross(zeta(1:3,i),tau2P(1:3,i))/norm(zeta(1:3,i))^2+lambP(i)*zeta(1:3,i);
end
%%
%
Wgl=kron(eye(N),Wred);
Hol=Wgl*zeta(:)<zeros(size(bqp));
if isempty(find(Hol==0))
    disp('Inside CWC - Succed!')
else
    disp('Not inside CWC')
end

[~,Sm]=lqr(Apre,Bpre,Qm,Rm);
K=lqr(Apre,Bpre,Sm,0.001*diag([1 1 1 1 1 1]));

x(:,1)=[L0;r0;dr0];
x(:,1)=x(:,1)+0.03*x(:,1).*(2*rand(9,1)-1);
x(:,1)=[
         0
         0
         0
   -0.4049
   -0.5956
    1.1282
    2.0349
    2.0750
   -1.0025];

for i=1:N
    xm=[x(1:3,i)-m*cross(x(4:6,i),x(7:9,i));x(4:6,i);x(7:9,i)]-[zeros(3,1);rfd;zeros(3,1)];
    v4lamfeed=(Wred*u(:,i)-b3cwcred)./(Wred*K*(xm-xtraj(:,i)));
    lambfeed(i)=min(min(v4lamfeed(find(v4lamfeed>0))),1)*1;%1;%
    uc(:,i)=u(:,i)-lambfeed(i)*K*(xm-xtraj(:,i));%lambfeed(i)*
    zetac(:,i)=uc(:,i)-[m*gv;m*Mgcr*rfd];
    %zetac(:,i)=zeta(:,i);
    x(:,i+1)=x(:,i)+T*[zetac(4:6,i)+cross(x(4:6,i),zetac(1:3,i));x(7:9,i);1/m*zetac(1:3,i)+gv];%A_sys*x(:,i)+B_sys*u(:,i);
    ddx4p(:,i)=1/m*zetac(1:3,i)+gv;
    if max(Wred*zetac(:,i)>0)==1
        disp('error')
        %break
    end
end

rCMP=x(4:6,1:N)-(ones(3,1)*((nv'*x(4:6,1:N))./(nv'*zeta(1:3,:)))).*zeta(1:3,:);

%plot(x')
figure(1)
subplot(311)
plot(T*(1:size(x,2)),x(1:3,:)')
legend('Lx','Ly','Lz')
subplot(312)
plot(T*(1:size(x,2)),x(4:6,:)')
legend('x','y','z')
subplot(313)
plot(T*(1:size(x,2)),x(7:9,:)')
legend('$\dot{x}$','$\dot{y}$','$\dot{z}$')
%
figure(2)
subplot(131)
title('Forces'),hold on, grid on
plot(T*(1:size(u,2)),zeta(1:3,:)')
legend('$f_x$','$f_y$','$f_z$')
subplot(132)
title('Torques'),hold on, grid on
plot(T*(1:size(u,2)),zeta(4:6,:)')
legend('$\tau_x$','$\tau_y$','$\tau_z$')
subplot(133)
title('CoM Accelerations'),hold on, grid on
plot(T*(1:size(u,2)),ddx4p(1:3,:)')
legend('$\ddot{x}$','$\ddot{y}$','$\ddot{z}$')
%
figure(3)
%
plot3(x(4,:),x(5,:),x(6,:)); hold on, grid on
plot3(x(4,1),x(5,1),x(6,1),'xr')
plot3(x(4,end),x(5,end),x(6,end),'ok')
plot3([Xfo Xfo -Xfo -Xfo Xfo],[Yfo -Yfo -Yfo Yfo Yfo],zeros(5,1),'r')
axis([-0.6 0.6 -0.6 0.6 -0.1 1.3])
pause
for i=1:N
plot3(x(4,i),x(5,i),x(6,i),'bx');
plot3(-rP(1,i),-rP(2,i),-rP(3,i),'bx');
plot3(rCMP(1,i),rCMP(2,i),rCMP(3,i),'ro');
%pause
end
plot3(rCMP(1,:),rCMP(2,:),rCMP(3,:),'r');
xlabel('x')
ylabel('y')
zlabel('z')
%
figure(4)
subplot(311)
title('CoM'),hold on, grid on
plot(T*(1:size(x,2)),x(4:6,:)')
legend('x','y','z')
subplot(312)
title('CoP'),hold on, grid on
plot(T*(1:size(rP,2)),-rP(1:3,:)')
legend('rCoPx','rCoPy','rCoPz')
subplot(313)
title('CMP'),hold on, grid on
plot(T*(1:size(rP,2)),rCMP(1:3,:)')
legend('rCMPx','rCMPy','rCMPz')
figure(5)
title('$\lambda$ Feedback'),hold on, grid on
plot(T*(1:size(lambfeed,2)),lambfeed)
legend('$\lambda$')