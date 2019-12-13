%clear all, close all, clc
ksys=39.478;
bsys=6.283;
vmax=200;
dt=10^-5;
tmax=1;
kp=300;kd=20;
Km=[kp kd];
t=dt*[1:tmax/dt];
e0=[0;0.08];
xa=[e0;-Km*e0];
Aclsys=[0 1;-ksys -bsys];Aa=[Aclsys,[0;ksys];zeros(1,3)];Ba=[0;0;1];
for i=1:length(t)
    v=-vmax*sign(xa(3)+Km*xa(1:2));
    xa=xa+(Aa*xa+Ba*v)*dt;
    Xa(:,i)=xa;
end
D(:,i)=xa;
%%
% figure(1)
% plot(t,Xa'), hold on,
figure(2)
plot3(Xa(1,:),Xa(2,:),Xa(3,:)), hold on, grid on
plot3(Xa(1,1),Xa(2,1),Xa(3,1),'x')
[xx,yy]=ndgrid(linspace(min(Xa(1,:)),max(Xa(1,:)),10),linspace(min(Xa(2,:)),max(Xa(2,:)),10));
zz=(-Km(1)*xx-Km(2)*yy);
surf(xx,yy,zz)