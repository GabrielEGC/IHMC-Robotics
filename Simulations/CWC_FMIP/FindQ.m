clc, clear all, close all
k=39.478;
b=6.283;
kp=300;kd=20;
Km=-[kp kd];
A=[0 1;-(k+kp) -(b+kd)];
nbar=(Km*A)'/norm((Km*A));
areacoefmin= 2.5180e-04;
Qmin=[5.649875618621051   0.092331774758496
   0.092331774758496   0.001508913251155];%eye(2);
for j=0:4
for jc=1:5
for i=1:10000
    m=2*rand-1;
    Q=Qmin+10^-j*[2*rand-1 m; m 2*rand-1];
    if det(Q)>0
    P=lyap(A',Q);
    if det(P*(nbar'*P^-1*nbar))<areacoefmin
        areacoefmin=det(P*(nbar'*P^-1*nbar))
        imin=i;
        Qmin=Q;
    end 
    end
end
end
end