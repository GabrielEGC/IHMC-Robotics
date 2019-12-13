function dXdt = VHIPVarCoP(t,X,zf,rPd,cl)
g=9.81;
gv=[0;0;-g];
u=0;
rP=rPd;
%Koolen
xro=X(2)-rP(2);
dxro=X(5);
zro=X(3)-rP(3);
dzro=X(6);

if cl==1
    %Koolen
Ts=-xro/dxro;
Zc=Ts*dzro+zro-g/2*Ts^2;
a=-1/Ts;
b=(Zc+g/2*Ts^2)/Ts;
u=-7*a^2+(3*zf*a^3-g*a)/b-10*a^3*b/g;
else
    %Garcia
umax=10000;
Ts=-xro/dxro*sqrt(umax);
Zc=zro*umax/g-xro/dxro*dzro*sqrt(umax)*sqrt(umax)/g-1/2*(xro/dxro*sqrt(umax))^2;
k2=(sqrt(zf*umax/g)+1)/2;
k=1+1./(Ts-1);
u=min(max(((1+k.*Ts.*(1-Ts))*k2+(k.*Ts+1).*(Zc-1/2))./(Ts.^4.*(1/2+sqrt(zf*umax/g)./(2*Ts.^2))),0),1);
u=u*umax;
end
u=max(u,0);
dXdt = [X(4:6);u*(X(1:3)-rP)+gv];