clear all, close all, clc
xev=linspace(0,4,20);%%%-4*eps^(1/1.5),4*eps^(1/1.5),20
yev=linspace(-.1,.1,5);
% xev=linspace(0,6,5);%%%-4*eps^(1/1.5),4*eps^(1/1.5),20
% yev=linspace(0,.9,5);
[Xev,Yev]=meshgrid(xev,yev);Zev=zeros(size(Xev));
for i=1:length(yev)
    for j=1:length(xev)
        [v,~,z]=computevcceq2([Xev(i,j);0.1875*(Xev(i,j)-0.2)+Yev(i,j)]);
        ZOev(i,j)=v;
        ZCev(i,j)=norm(z);
    end
end
figure;hold on; grid on
subplot(211);surf(Xev,Yev,ZOev)
subplot(212);surf(Xev,Yev,ZCev)
figure;hold on; grid on
subplot(211);mesh(Xev,Yev,ZOev)
subplot(212);mesh(Xev,Yev,ZCev)
%axis([xev(1) xev(end) yev(1) yev(end) 0 0.005])