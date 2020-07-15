function PlotSimHybr_FixH_BS_enh_NoMirr_enh(x,alr,simp,n,spP)%,x,xParam
    close all%,clc
    InitSetUp2;
    [Fz0,~,~] = LiftOffEventsFcn(0,xSt0,@Fextfun,Lrv,Aznom,kpLe,kpAz,tau0,l0,I_bod,kdtau,kdLe,kdAz);
    if (norm(xSt0)>1e+12||Fz0(1)<1e-10)%||abs(det(MFFeval))<0.0001
        warning('Invalid initial conditions');
    end
rr=cell(n+1,1);
ESt=cell(n+1,1);
tFl=cell(n,1);
xFl=cell(n,1);
rFlH=cell(n,1);
rBrCoMFl=cell(n,1);
GFLStrp=cell(n,1);
rrFlL=cell(n,1);
rrFlR=cell(n,1);
EFl=cell(n,1);
tSt=cell(n+1,1);
xSt=cell(n+1,1);
rrH=cell(n+1,1);
rrL=cell(n+1,1);
rrR=cell(n+1,1);
rBrCoM=cell(n+1,1);
    
OptStan = odeset('Events', @(t,xF) LiftOffEventsFcn(t,xF,@Fextfun,Lrv,Aznom,kpLe,kpAz,tau0,l0,I_bod,kdtau,kdLe,kdAz));%'AbsTol',1e-10,
OptStanNoMirr = odeset('Events', @(t,xF) LiftOffEventsFcn(t,xF,@Fextfun_NoMirr,[-Lrv(1);Lrv(2:3)],Aznom,kpLe,kpAz,tau0,l0,I_bod,kdtau,kdLe,kdAz));%'AbsTol',1e-10,
[tSt{1},xSt{1}] = ode45(@(t,x) stance2DOut_full(t,x,I_bod,Lrv,Aznom,kpLe,kpAz,tau0,l0,gv,m,kdtau,kdLe,kdAz), [0 tpHM0], xSt0, OptStan);%
%[tSt{1},xSt{1}] = ode45(@(t,x) sysDdeqdt(t,x), linspace(0,tpHM0,10000), xSt0, OptStan);%
xStF=xSt{1}(end,:);
rr{1}=xSt{1}(:,1:3);
i=0;

for jRot = 1:size(xSt{1},1)
    eul2rotjR=eul2rotm1([xSt{i+1}(jRot,4),xSt{i+1}(jRot,5),xSt{i+1}(jRot,6)]);
    rloc=eul2rotjR'*([0;0;0]-xSt{i+1}(jRot,1:3)')-Lrv;
    rrH{1}(jRot,1:3)=xSt{1}(jRot,1:3)'+eul2rotjR*[0;0;Lrv(3)];
    rrL{1}(jRot,1:3)=xSt{1}(jRot,1:3)'+eul2rotjR*[-Lrv(1);Lrv(2:3)];
    rrR{1}(jRot,1:3)=xSt{1}(jRot,1:3)'+eul2rotjR*Lrv;
    rBrCoM{1}(:,jRot)=eul2rotjR'*([0;0;0]-xSt{1}(jRot,1:3)');%rloc%abs(atan2(rloc(1),-rloc(3))+(2*mod(i,2)-1)*Aznom)
    ESt{i+1}(jRot,1)=m/2*(xSt{i+1}(jRot,7:9)*xSt{i+1}(jRot,7:9)')+1/2*(xSt{i+1}(jRot,10:12)*eul2rotjR*(I_bod^-1)*eul2rotjR'*xSt{i+1}(jRot,10:12)')+m*(-gv(3))*xSt{i+1}(jRot,3)+1/2*kpLe*(norm(rloc)-l0)^2+1/2*kpAz*(atan2(rloc(1),-rloc(3))+(2*mod(i,2)-1)*Aznom)^2;
end
EstV=eqEnergy(xSt{i+1}(:,10),xSt{i+1}(:,11),xSt{i+1}(:,12),xSt{i+1}(:,7),xSt{i+1}(:,8),xSt{i+1}(:,9),xSt{i+1}(:,6),xSt{i+1}(:,4),xSt{i+1}(:,1),xSt{i+1}(:,2),xSt{i+1}(:,3),xSt{i+1}(:,5));
norm(EstV-ESt{i+1}(:,1))
for i=1:n
    i
    sid=(2*mod(i,2)-1);
    phasphiStfin=castangphi(xStF',[sid*Lrv(1);Lrv(2:3)]);
    if i>1
        xF0=xStF+[GFLStrp{i-1}(end,1:3),zeros(1,9)];
    else
        xF0=xStF;
    end
    OptFl = odeset('Events', @(t,xF) TouchDoEventsFcnStride(t,xF,-sid*Aznom,l0,phasphi,[-sid*Lrv(1);Lrv(2:3)],alr,phasphiStfin));%'AbsTol',1e-10,%,Lr,l0,q5j,dq5j,av
    [tFl{i},xFl{i}] = ode45(@(t,x) flight2DOut_full(t,x,gv,I_bod), [0 tpHM0], xF0, OptFl);%linspace(0,tpHM0,10000)
    if (tFl{i}(end)==tpHM0)
        i
        error('no TD in tpHM0');
    end
    
    for jRot = 1:size(xFl{i},1)
        eul2rotjR=eul2rotm1([xFl{i}(jRot,4),xFl{i}(jRot,5),xFl{i}(jRot,6)]);
        rFlH{i}(jRot,1:3)=xFl{i}(jRot,1:3)+(eul2rotjR*[0;0;Lrv(3)])';
        phasphiTD=phasphiStfin+alr+phasphi(2)*tFl{i}(jRot); %
        rBrCoMFl{i}(:,jRot)=roty(sid*Aznom)*[0;l0*sin(phasphiTD);-l0*cos(phasphiTD)]+[-sid*Lrv(1);Lrv(2:3)];%OK?
        GFLStrp{i}(jRot,1:3)=xFl{i}(jRot,1:3)+(eul2rotjR*(rBrCoMFl{i}(:,jRot)))';
        rrFlL{i}(jRot,1:3)=xFl{i}(jRot,1:3)+(eul2rotjR*[-Lrv(1);Lrv(2:3)])';
        rrFlR{i}(jRot,1:3)=xFl{i}(jRot,1:3)+(eul2rotjR*Lrv)';


        EFl{i}(jRot,1)=m/2*(xFl{i}(jRot,7:9)*xFl{i}(jRot,7:9)')+1/2*(xFl{i}(jRot,10:12)*eul2rotjR*(I_bod^-1)*eul2rotjR'*xFl{i}(jRot,10:12)')+m*(-gv(3))*xFl{i}(jRot,3);
    end
    rf=xFl{i}(end,1:3)';drf=xFl{i}(end,7:9)';eulf=xFl{i}(end,4:6);Lf=xFl{i}(end,10:12)';
     xFlF=[(rf-GFLStrp{i}(end,1:3)');eulf';drf;Lf];
    
    xSt0=xFlF;
    if mod(i,2)==1
        [tSt{i+1},xSt{i+1}] = ode45(@(t,x) stance2DOut_full_NoMirr(t,x,I_bod,Lrv,Aznom,kpLe,kpAz,tau0,l0,gv,m,kdtau,kdLe,kdAz), [0 tpHM0], xSt0, OptStanNoMirr);%
    else
        [tSt{i+1},xSt{i+1}] = ode45(@(t,x) stance2DOut_full(t,x,I_bod,Lrv,Aznom,kpLe,kpAz,tau0,l0,gv,m,kdtau,kdLe,kdAz), [0 tpHM0], xSt0, OptStan);%
    end
    xStF=xSt{i+1}(end,:);
    rr{i+1}=GFLStrp{i}(end,1:3)+xSt{i+1}(:,1:3);
    for jRot = 1:size(xSt{i+1},1)
        eul2rotjR=eul2rotm1([xSt{i+1}(jRot,4),xSt{i+1}(jRot,5),xSt{i+1}(jRot,6)]);
        rloc=eul2rotjR'*([0;0;0]-xSt{i+1}(jRot,1:3)')-[-sid*Lrv(1);Lrv(2:3)];
        %rrH{i+1}(jRot,1:3)=xSt{i+1}(jRot,1:3)+(eul2rotm1([xSt{i+1}(jRot,4),xSt{i+1}(jRot,5),xSt{i+1}(jRot,6)])*[-sid*Lrv(1);Lrv(2:3)])'+GFLStrp{i}(end,1:3);
        rrH{i+1}(jRot,1:3)=xSt{i+1}(jRot,1:3)+(eul2rotjR*[0;0;Lrv(3)])'+GFLStrp{i}(end,1:3);
        rrL{i+1}(jRot,1:3)=xSt{i+1}(jRot,1:3)+(eul2rotjR*[-Lrv(1);Lrv(2:3)])'+GFLStrp{i}(end,1:3);
        rrR{i+1}(jRot,1:3)=xSt{i+1}(jRot,1:3)+(eul2rotjR*Lrv)'+GFLStrp{i}(end,1:3);
        rBrCoM{i+1}(:,jRot)=eul2rotm1([xSt{i+1}(jRot,4),xSt{i+1}(jRot,5),xSt{i+1}(jRot,6)])'*([0;0;0]-xSt{i+1}(jRot,1:3)');
        ESt{i+1}(jRot,1)=m/2*(xSt{i+1}(jRot,7:9)*xSt{i+1}(jRot,7:9)')+1/2*(xSt{i+1}(jRot,10:12)*eul2rotjR*(I_bod^-1)*eul2rotjR'*xSt{i+1}(jRot,10:12)')+m*(-gv(3))*xSt{i+1}(jRot,3)+1/2*kpLe*(norm(rloc)-l0)^2+1/2*kpAz*(atan2(rloc(1),-rloc(3))+sid*Aznom)^2;
    end
end
% xrv=[xr{1};xFl{1}(2:end,1);xr{2}(2:end)];%
% yrv=[yr{1};xFl{1}(2:end,2);yr{2}(2:end)];%
% t=[tSt{1};tSt{1}(end)+tFl{1}(2:end);tSt{1}(end)+tFl{1}(end)+tSt{2}(2:end)];%
% xB=[0*tSt{1};GFLStxp{1};GFLStxp{1}(end)*ones(size(tSt{2}(2:end)))];
% yB=[0*tSt{1};GFLStyp{1};GFLStyp{1}(end)*ones(size(tSt{2}(2:end)))];


%xrv=[xr{1}];yrv=[yr{1}];zrv=[zr{1}];
%xrHv=[xrH{1}];yrHv=[yrH{1}];zrHv=[zrH{1}];

nt=length(tSt{1});
for i=1:n %%%GOOOOO BACKKKKKK
    nt=nt+length(tFl{i})-2+length(tSt{i+1});
end

rrv=rr{1};
rrHv=rrH{1};
rrHLv=rrL{1};
rrHRv=rrR{1};
rrHLRFv=rrHRv;
t=tSt{1};
tbFl=zeros(size(t));
rB=zeros(length(t),3);
rBrCoMv=rBrCoM{1}';
TE=[ESt{1}];
%tswHM=[0;tSt{1}(end)];
for i=1:n %%%GOOOOO BACKKKKKK
    %tswHM=[tswHM;tswHM(end)+tFl{i}(end);tswHM(end)+tFl{i}(end)+tSt{i+1}(end)];
    rrv=[rrv;xFl{i}(2:end,1:3);rr{i+1}(2:end,1:3)];%
    rrHv=[rrHv;rFlH{i}(2:end,1:3);rrH{i+1}(2:end,1:3)];%
    rrHLv=[rrHLv;rrFlL{i}(2:end,1:3);rrL{i+1}(2:end,1:3)];%
    rrHRv=[rrHRv;rrFlR{i}(2:end,1:3);rrR{i+1}(2:end,1:3)];%
    rrHLRFv=[rrHLRFv;mod(i,2)*rrFlL{i}(2:end,1:3)+(1-mod(i,2))*rrFlR{i}(2:end,1:3);mod(i,2)*rrL{i+1}(2:end,1:3)+(1-mod(i,2))*rrR{i+1}(2:end,1:3)];%
    tbFl=[tbFl;ones(size(tFl{i}(2:end)));zeros(length(tSt{i+1})-1,1)];
    t=[t;t(end)+tFl{i}(2:end);t(end)+tFl{i}(end)+tSt{i+1}(2:end)];
    rB=[rB;GFLStrp{i}(2:end,1:3);ones(length(tSt{i+1})-1,1)*GFLStrp{i}(end,1:3)];
    rBrCoMv=[rBrCoMv;rBrCoMFl{i}(:,2:end)';rBrCoM{i+1}(:,2:end)'];%verify
    TE=[TE;EFl{i}(2:end);ESt{i+1}(2:end)];
end
%rrHLRFv=rrHLv.*tbFl+rrHRv.*(1-tbFl);
% disp('hi')
% length(t)
% nt
% length(xB)
% length(xrv)
% length(xrHv)
% size(rBrCoMv,1)
% length(TE)
figure(3)
title('Total Energy','Interpreter','latex');hold on
%subplot(211)
plot(t,TE,'-x')
xlabel('time(s)','Interpreter','latex')
ylabel('Energy(J)','Interpreter','latex')
%plot(t,tbFl);hold on
figure(1)
figure(2)
title('Leg Swing Trajectory  - Body Coordinates','Interpreter','latex'); hold on, grid on
plot3(0,0,0,'xr')
plot3(Lrv(1),Lrv(2),Lrv(3),'ok')
plot3(-Lrv(1),Lrv(2),Lrv(3),'ok')
phasp=linspace(-pi/2,pi/2,100);
nomtraj=[0*phasp;l0*sin(phasp);-l0*cos(phasp)];
nomtrajrn=roty(-Aznom)*nomtraj+Lrv*ones(1,length(phasp));
nomtrajrp=roty(Aznom)*nomtraj+[-Lrv(1);Lrv(2:3)]*ones(1,length(phasp));
plot3(nomtrajrn(1,:),nomtrajrn(2,:),nomtrajrn(3,:),'-b')
plot3(nomtrajrp(1,:),nomtrajrp(2,:),nomtrajrp(3,:),'-b')
axis([-1 1 -1 1 -1 1])
xlabel('$x(m)$','Interpreter','latex');ylabel('$y(m)$','Interpreter','latex');zlabel('$z(m)$','Interpreter','latex');
view(135,50)
drawnow
axvg=[-1.1 .2 -.5 2 0 1.2];
if(simp)
i=1;
tic;
while(i~=length(t)&&t(i)~=t(end))
figure(1)
clf;
%plot(xw,yw)
%subplot(221)
plot3(rrv(:,1),rrv(:,2),rrv(:,3),'-'),hold on, grid on
a=toc;
[~,i]=min(abs(t-spP*a));
%axvg=[rrv(i,1)-1.75 rrv(i,1)+1.75 rrv(i,2)-1.75 rrv(i,2)+.175 0 2.5];
axvg=[rrv(i,1)-.75 rrv(i,1)+.75 rrv(i,2)-.75 rrv(i,2)+.75 0 1];
%axvg=[xrv(i)-.5 xrv(i)+.5 yrv(i)-1 yrv(i)+1 zrv(i)-.7 zrv(i)+.3];
%axvg=[-.3 .2 yrv(i)-1 yrv(i)+1 zrv(i)-.7 zrv(i)+.3];
axis(axvg);%
%axis([xrv(i)-1 xrv(i)+1 yrv(i)-1 yrv(i)+1 zrv(i)-1 zrv(i)+1])
plot3(rB(i,1),rB(i,2),rB(i,3),'ok');
plot3(rrHv(i,1),rrHv(i,2),rrHv(i,3),'ok');
plot3(rrHLv(i,1),rrHLv(i,2),rrHLv(i,3),'ok');
plot3(rrHRv(i,1),rrHRv(i,2),rrHRv(i,3),'ok');
plot3(rrv(i,1),rrv(i,2),rrv(i,3),'xr');%plot(xw(i),yw(i),'xb');
%plot3([rrHLv(i,1) rrHRv(i,1)],[rrHLv(i,2) rrHRv(i,2)],[rrHLv(i,3) rrHRv(i,3)],'r','LineWidth',3),
plot3([rrHLRFv(i,1) rB(i,1)],[rrHLRFv(i,2) rB(i,2)],[rrHLRFv(i,3) rB(i,3)],'LineWidth',6),
plot3([rrv(i,1) rrHRv(i,1)],[rrv(i,2) rrHRv(i,2)],[rrv(i,3) rrHRv(i,3)],'r','LineWidth',3),
plot3([rrHLv(i,1) rrv(i,1)],[rrHLv(i,2) rrv(i,2)],[rrHLv(i,3) rrv(i,3)],'r','LineWidth',3),

plot3([rrv(i,1) rrHv(i,1)],[rrv(i,2) rrHv(i,2)],[rrv(i,3) rrHv(i,3)],'c','LineWidth',3),
%plot3([rrHLv(i,1) rrHv(i,1)],[rrHLv(i,2) rrHv(i,2)],[rrHLv(i,3) rrHv(i,3)],'c','LineWidth',3),
%plot3([rrHRv(i,1) rrHv(i,1)],[rrHRv(i,2) rrHv(i,2)],[rrHRv(i,3) rrHv(i,3)],'c','LineWidth',3),
%plot([xr(i) xw(i)],[yr(i) yw(i)],'LineWidth',3);
xlabel('$x(m)$','Interpreter','latex');ylabel('$y(m)$','Interpreter','latex');zlabel('$z(m)$','Interpreter','latex');title('FastRunner 3D','Interpreter','latex');%title('Isometric View')
%view(135,20)
%view(100,0)
view(120,20)
%view(0,10)
% 
% 
%     subplot(223)
%     plot3(xrv,yrv,zrv,'-'),hold on, grid on
%     plot3(xB(i),yB(i),zB(i),'ok');
%     plot3(xrHv(i),yrHv(i),zrHv(i),'or');
%     plot3(xrv(i),yrv(i),zrv(i),'or');%plot(xw(i),yw(i),'xb');
%     plot3([xrHv(i) xB(i)],[yrHv(i) yB(i)],[zrHv(i) zB(i)],'LineWidth',3),
%     plot3([xrv(i) xrHv(i)],[yrv(i) yrHv(i)],[zrv(i) zrHv(i)],'LineWidth',3),%plot([xr(i) xw(i)],[yr(i) yw(i)],'LineWidth',3);
%     xlabel('$x(m)$');ylabel('$y(m)$');zlabel('$z(m)$');
%     
%     hold off
%     axis(axvg), 
%     view(90,0)
%     title('Lateral View')
%     
%     subplot(222)
%     plot3(xrv,yrv,zrv,'-'),hold on, grid on
%     plot3(xB(i),yB(i),zB(i),'ok');
%     plot3(xrHv(i),yrHv(i),zrHv(i),'or');
%     plot3(xrv(i),yrv(i),zrv(i),'or');%plot(xw(i),yw(i),'xb');
%     plot3([xrHv(i) xB(i)],[yrHv(i) yB(i)],[zrHv(i) zB(i)],'LineWidth',3),
%     plot3([xrv(i) xrHv(i)],[yrv(i) yrHv(i)],[zrv(i) zrHv(i)],'LineWidth',3),%plot([xr(i) xw(i)],[yr(i) yw(i)],'LineWidth',3);
%     xlabel('$x(m)$');ylabel('$y(m)$');zlabel('$z(m)$');
%     hold off
%      axis(axvg), 
%      view(0,0)
%      title('Frontal View')
% 
%     subplot(224)
%     plot3(xrv,yrv,zrv,'-'),hold on, grid on
%     plot3(xB(i),yB(i),zB(i),'ok');
%     plot3(xrHv(i),yrHv(i),zrHv(i),'or');
%     plot3(xrv(i),yrv(i),zrv(i),'or');%plot(xw(i),yw(i),'xb');
%     plot3([xrHv(i) xB(i)],[yrHv(i) yB(i)],[zrHv(i) zB(i)],'LineWidth',3),
%     plot3([xrv(i) xrHv(i)],[yrv(i) yrHv(i)],[zrv(i) zrHv(i)],'LineWidth',3),%plot([xr(i) xw(i)],[yr(i) yw(i)],'LineWidth',3);
%     xlabel('$x(m)$');ylabel('$y(m)$');zlabel('$z(m)$');
%     hold off
%         axis(axvg), 
%         view(0,90)
%         title('Top View')
% 
% 
% figure(2)
% if tbFl(i)
%     plot3(rBrCoMv(i,1),rBrCoMv(i,2),rBrCoMv(i,3),'xr')
% else
%     plot3(rBrCoMv(i,1),rBrCoMv(i,2),rBrCoMv(i,3),'xb')
% 
% rloc=rBrCoMv(i,:)'-[0;Lrv(2:3)];
% drloc=0*rloc;%-rZYXR'*(dr0+cross(omgc,([0;0;0]-r0)));
% azhat=[-rloc(3);0;rloc(1)];azhat=azhat/norm(azhat)^2;
% ttauhat=[-rloc(1)*rloc(2); rloc(1)^2 + rloc(3)^2; -rloc(2)*rloc(3)];ttauhat=ttauhat/(norm(ttauhat)*norm(rloc));%
% 
% dphasest=ttauhat'*drloc;
% dnormrloc=rloc'/norm(rloc)*drloc;
% dazrloc=azhat'*drloc;
% 
% Fl=(-kpLe*(norm(rloc)-l0)-kdLe*dnormrloc)*rloc/norm(rloc);
% FAz=(-kpAz*(atan2(rloc(1),-rloc(3))-Aznom)-kdAz*dazrloc)*azhat;
% Ftau=(-kdtau*dphasest+tau0(1))*ttauhat;%%
% 
% 
% quiver3(rBrCoMv(i,1),rBrCoMv(i,2),rBrCoMv(i,3),Fl(1),Fl(2),Fl(3),0.1,'r');
% %norm(rloc)
% quiver3(rBrCoMv(i,1),rBrCoMv(i,2),rBrCoMv(i,3),FAz(1),FAz(2),FAz(3),0.1,'g');
% quiver3(rBrCoMv(i,1),rBrCoMv(i,2),rBrCoMv(i,3),Ftau(1),Ftau(2),Ftau(3),.1,'b');
% drawnow
%end


end
end
%plot(t,-0.7+0.1*sin(2*t))
%axis([0 5 -1 -.5])
%lg211=legend('$y_w$','$y_{wref}$');set(lg211,'FontSize',12);
%ylabel('$y(m)$');
end
function dydx = sysDdeqdt(x,y)%VERIFY L0P
    r0x=y(1,:);r0y=y(2,:);r0z=y(3,:);
    psi=y(4,:);th=y(5,:);phi=y(6,:);
    dr0x=y(7,:);dr0y=y(8,:);dr0z=y(9,:);
    L0p1=y(10,:);L0p2=y(11,:);L0p3=y(12,:);   
    dydx=eqdSRBMdt(L0p1,L0p2,L0p3,dr0x,dr0y,dr0z,phi,psi,r0x,r0y,r0z,th);
end