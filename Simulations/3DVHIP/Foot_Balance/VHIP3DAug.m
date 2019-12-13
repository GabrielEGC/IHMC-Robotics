clc, clear all, close all
Xfo=0.07;
Yfo=0.1;
rfo=[0;-0.5;0];
r0 = [0.3;-1.2;1.2];
dr0= 0.8*[-0.7;2;-0.2];
r0 = [0.0114
   -0.5105
    1.0919];
dr0= [-0.0171
    0.4000
   -0.2000];   
%%
g=9.81;
gv=[0;0;-g];
taub=0:0.05:1;
xb0=r0*ones(size(taub))+dr0*taub+gv*taub.^2/2;
xICC=r0*ones(size(taub))+dr0*taub+gv*taub.^2;
pRf=rotx(15)*roty(15)*rotz(15);
nf=pRf*[0;0;1];   
%%
pr0=r0;
pdr0=dr0;
pr01=r0+dr0;
pr0(3)=(nf'*rfo-nf'*[pr0(1:2);0])/nf(3);
pr01(3)=(nf'*rfo-nf'*[pr01(1:2);0])/nf(3);

pr0=pRf^-1*(pr0-rfo);pr0=pr0(1:2);
pr01=pRf^-1*(pr01-rfo);pr01=pr01(1:2);
pdr0=pr01-pr0;
lambV=[([Xfo;-Xfo]-pr0(1))/pdr0(1);([Yfo;-Yfo]-pr0(2))/pdr0(2)];
lambV=sort(lambV)
tausome=(7*lambV(2)+lambV(3))/8;
zro=nf'*(r0-rfo)/nf(3);
dzro=nf'*dr0/nf(3);
tauICP=(dzro+sqrt(dzro^2+4*g*zro))/(2*g)
tauzcrit=(dzro+sqrt(dzro^2+2*g*zro))/g;
taumin=max(lambV(2),0);
taumax=min(lambV(3),tauzcrit);

rICC=r0+dr0*tauICP+gv*tauICP.^2;

rCBPprev=r0+dr0*tauICP+gv*tauICP.^2/2;
rIBT=r0+dr0*tauzcrit+gv*tauzcrit.^2/2;

tausome=min(taumin+0.95*(taumax-taumin),max(taumin+0.05*(taumax-taumin),tauICP));
rCBPnew=r0+dr0*tausome+gv*tausome.^2/2;

x(:,1)=[r0;dr0;tausome];
%%
xb0some=r0+dr0*tausome+gv*tausome.^2/2;

zsome=(nf'*rfo-nf'*[xb0some(1:2);0])/nf(3);
xb0some(3)=zsome;
xb0some=r0+dr0*tausome+gv*tausome.^2/2;
zf=1.1;
rPd=xb0some;
rPd=rfo;

%%

[t1,X1] = ode45(@(t,X) VHIPVarCoPAug(t,X,zf,rPd,nf,rfo,pRf,1),[0 3],x);
% xn=X1(end,:)';xn(5)=0.4;xn(6)=-0.2;
% 
% ra0=xn(1:3);
% dra0=xn(4:6);
% 
% pr0=ra0;
% pdr0=dra0;
% pr01=pr0+pdr0;
% pr0(3)=(nf'*rfo-nf'*[pr0(1:2);0])/nf(3);
% pr01(3)=(nf'*rfo-nf'*[pr01(1:2);0])/nf(3);
% 
% pr0=pRf^-1*(pr0-rfo);pr0=pr0(1:2);
% pr01=pRf^-1*(pr01-rfo);pr01=pr01(1:2);
% pdr0=pr01-pr0;
% lambV=[([Xfo;-Xfo]-pr0(1))/pdr0(1);([Yfo;-Yfo]-pr0(2))/pdr0(2)];
% lambV=sort(lambV)
% zro=nf'*(ra0-rfo)/nf(3);
% dzro=nf'*dra0/nf(3);
% tauICP=(dzro+sqrt(dzro^2+4*g*zro))/(2*g)
% tauzcrit=(dzro+sqrt(dzro^2+2*g*zro))/g
% taumin=max(lambV(2),0);
% taumax=min(lambV(3),tauzcrit);
% 
% tausome=min(taumin+0.95*(taumax-taumin),max(taumin+0.05*(taumax-taumin),tauICP));
% xn(7)=tausome;
% 
% zf=1.1;

[t2,X2] = ode45(@(t,X) VHIPVarCoPAug(t,X,zf,rPd,nf,rfo,pRf,2),[0 3],x);
rsim1=X1(:,1:3)';
rsim2=X2(:,1:3)';
%%
kpP=0.5;
nfn=nf/nf(3);
Tg=X1(:,7)';
Zcg=nfn'*(X1(:,1:3)'-rfo*ones(1,size(X1,1)))+nfn'*(X1(:,4:6)'.*(ones(3,1)*Tg))-g/2*Tg.^2;
XiXY=X1(:,1:2)'+X1(:,4:5)'.*(ones(2,1)*Tg);
rP=XiXY+kpP*(XiXY-rfo(1:2)*ones(1,size(XiXY,2)));
rP(3,:)=nfn'*rfo-nfn(1:2)'*rP(1:2,:);
%%
figure(1)
XiXY3=X1(:,1:3)'+X1(:,4:6)'.*(ones(3,1)*Tg);
plot3(XiXY3(1,:),XiXY3(2,:),XiXY3(3,:)); 

%%
figure(2)
kpP=0.5;
nfn=nf/nf(3);
Tg=X1(:,7)';
Zcg=nfn'*(X1(:,1:3)'-rfo*ones(1,size(X1,1)))+nfn'*(X1(:,4:6)'.*(ones(3,1)*Tg))-g/2*Tg.^2;
XiXY=X1(:,1:2)'+X1(:,4:5)'.*(ones(2,1)*Tg);
rP=XiXY+kpP*(XiXY-rfo(1:2)*ones(1,size(XiXY,2)));
rP(3,:)=nfn'*rfo-nfn(1:2)'*rP(1:2,:);

rPproj=pRf^-1*(rP-rfo*ones(1,size(rP,2)));rPproj=rPproj(1:2,:);
for irk=1:size(rP,2)
    rPproji=rPproj(:,irk);
    lambV=[[Xfo;-Xfo]/rPproji(1);[Yfo;-Yfo]/rPproji(2)];
    if ~isempty(find(lambV<1&lambV>0,1))
    lambdV=min(lambV(find(lambV<1&lambV>0)));
    rP(:,irk)=pRf*[lambdV*rPproji;0]+rfo;
    end
end

plot3(rP(1,:),rP(2,:),rP(3,:)); 
xlabel('x')
ylabel('y')
zlabel('z')
%%
rCBP=X1(:,1:3)+(X1(:,7)*ones(1,3)).*X1(:,4:6)+X1(:,7).^2/2*gv';

figure(3)
plot3(x(1,1),x(2,1),x(3,1),'xr'); hold on, grid on
plot3(rsim1(1,:),rsim1(2,:),rsim1(3,:)); 
plot3(rsim2(1,:),rsim2(2,:),rsim2(3,:));
plot3(xb0(1,:),xb0(2,:),xb0(3,:));
plot3(xICC(1,:),xICC(2,:),xICC(3,:));
plot3(rP(1,:),rP(2,:),rP(3,:)); 
plot3(rsim1(1,end),rsim1(2,end),rsim1(3,end),'oc')
%plot3(rsim2(1,end),rsim2(2,end),rsim2(3,end),'xc')
plot3(rICC(1,:),rICC(2,:),rICC(3,:),'bo');
plot3(rCBPprev(1,:),rCBPprev(2,:),rCBPprev(3,:),'ro'); 
plot3(rCBPnew(1,:),rCBPnew(2,:),rCBPnew(3,:),'go'); 
%plot3(rCBP(:,1),rCBP(:,2),rCBP(:,3),'g');
%plot3(xb0some(1,:),xb0some(2,:),xb0some(3,:),'bx'); 
%plot3(rIBT(1,:),rIBT(2,:),rIBT(3,:),'bo');

R4pfo=pRf*[[Xfo Xfo -Xfo -Xfo Xfo];[Yfo -Yfo -Yfo Yfo Yfo];zeros(1,5)]+rfo*ones(1,5);
plot3(R4pfo(1,:),R4pfo(2,:),R4pfo(3,:),'r'), hold on
legend('Initial CoM position','Orbital Energy','Sliding Mode','IBT','ICC','CoP','Final Position','ICP','Previous CBP','New CBP')
axis([-0.5 0.5 -0.8 -0.2 -0.2 1.5])

xlabel('x')
ylabel('y')
zlabel('z')