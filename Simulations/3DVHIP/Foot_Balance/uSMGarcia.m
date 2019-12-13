function u=uSMGarcia(Tg,Zcg,zf,g)
umax=1000;
Ts=Tg*sqrt(umax);
Zc=Zcg*umax/g;%(zro-xro/dxro*dzro-g/2*(xro/dxro*)^2)*umax/g;
k2=(sqrt(zf*umax/g)+1)/2;
k=1+1./(Ts-1);
u=min(max(((1+k.*Ts.*(1-Ts))*k2+(k.*Ts+1).*(Zc-1/2))./(Ts.^4.*(1/2+sqrt(zf*umax/g)./(2*Ts.^2))),0),1);
u=u*umax;
end
