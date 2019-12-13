clear all, close all, clc
ksys=39.478;
bsys=6.283;
vmax=200;
dt=10^-6;
tmax=0.1;
kp=300;kd=20;
Km=[kp kd];
t=dt*[1:tmax/dt];
e1=0;
e2=0.15;
Aclsys=[0 1;-ksys -bsys];Aa=[Aclsys,[0;ksys];zeros(1,3)];Ba=[0;0;1];B=[0;ksys];
nbar=Km*(Aclsys-B*Km);
for i1=1:length(e1)
    i1
    for i2=1:length(e2)
        e0=[e1(i1);e2(i2)];
        xa=[e0;-Km*e0];
        for i=1:length(t)
            v=-vmax*sign(xa(3)+Km*xa(1:2));
            xa=xa+(Aa*xa+Ba*v)*dt;
            Xa(:,i)=xa;
        end
        D(i1,i2)=(norm(xa(1:2))<10^-6);
    end 
end
%%
% figure(1)
% plot(t,Xa'), hold on,
figure(2)
plot3(Xa(1,:),Xa(2,:),Xa(3,:)), hold on, grid on
plot3(Xa(1,1),Xa(2,1),Xa(3,1),'x')
[xx,yy]=ndgrid(linspace(min(Xa(1,:)),max(Xa(1,:)),10),linspace(min(Xa(2,:)),max(Xa(2,:)),10));
zz=(-Km(1)*xx-Km(2)*yy);

xpl=[min(Xa(1,:)) max(Xa(1,:))];ypl1=(vmax-nbar(1)*xpl)/nbar(2);
ypl2=(-vmax-nbar(1)*xpl)/nbar(2);
zpl1=(-Km(1)*xpl-Km(2)*ypl1);zpl2=(-Km(1)*xpl-Km(2)*ypl2);
plot3(xpl,ypl1,zpl1);plot3(xpl,ypl2,zpl2)


surf(xx,yy,zz)