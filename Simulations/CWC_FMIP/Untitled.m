clc, clear all, close all
xd=1;
k=39.478;
b=6.283;
wn=sqrt(39.478);
xi=b/(2*wn);
vmax=200;
kp=300;kd=20;
Km=-[kp kd];
%   
A=[0 1;-(k+k*kp) -(b+k*kd)];
nbar=(Km*A)';
Q=1*diag([1;1]);%%At will
%Q=[6.995903381204523   0.459335946720598
%   0.459335946720598   0.030159010858415];
P=lyap(A',Q);
areacoef=det(P*(nbar'*P^-1*nbar));

c0=vmax^2/2*(nbar'*P^-1*nbar)^-1;
%%
tplot=linspace(0,2*pi,100);cir=sqrt(2*c0)*[cos(tplot); sin(tplot)];
ellip=inv(chol(P))*cir;
%figure(1)
plot(ellip(1,:),ellip(2,:)), hold on
syms xp yp;
xpl=[-4:0.01:4]*10^-1;ypl1=(vmax-nbar(1)*xpl)/nbar(2);
ypl2=(-vmax-nbar(1)*xpl)/nbar(2);
plot(xpl,ypl1);plot(xpl,ypl2)