function dXdt = VHIPVarCoPAug(t,X,zf,rPd,nf,rfo,pRf,cl)

Xfo=0.07;
Yfo=0.1;

g=9.81;
gv=[0;0;-g];
u=0;
% xro=X(2)-rP(2);
% dxro=X(5);
% zro=X(3)-rP(3);
% dzro=X(6);
nfn=nf/nf(3);
kpP=0.5;

Tg=X(7);
Ts=Tg;
Zcg=nfn'*(X(1:3)-rfo)+Ts*nfn'*X(4:6)-g/2*Ts^2;
XiXY=X(1:2)+Ts*X(4:5);
rP=XiXY+kpP*(XiXY-rPd(1:2)) ;
rP(3)=nfn'*rfo-nfn(1:2)'*rP(1:2);

rPproj=pRf^-1*(rP-rfo);rPproj=rPproj(1:2);
lambV=[[Xfo;-Xfo]/rPproj(1);[Yfo;-Yfo]/rPproj(2)];

if ~isempty(find(lambV<1&lambV>0,1))
lambdV=min(lambV(find(lambV<1&lambV>0)));
rP=pRf*[lambdV*rPproj;0]+rfo;
end

if cl==1%Koolen
a=-1/Tg;
b=(Zcg+g/2*Tg^2)/Tg;
u=-7*a^2+(3*zf*a^3-g*a)/b-10*a^3*b/g;
else%Garcia
u=uSMGarcia(Tg,Zcg,zf,g);
end
u=max(u,0);

dXdt = [X(4:6);u*(X(1:3)-rP)+gv;-1+u*X(7)^2];