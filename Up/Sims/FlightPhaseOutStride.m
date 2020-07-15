function [tFlF,xFlF,phasphiTD] = FlightPhaseOutStride(xStF,Aznom,l0,phasphi,Lrv,gv,I_bod,tpHM0,alr,phasphiStfin)%,varargin
zv=xStF([3 9]);
L0=xStF(10:12);
Aznom=-Aznom;
Opt = odeset('Events', @(t,xF) TouchDoEventsFcnStride_2(t,xF,Aznom,l0,phasphi,Lrv,alr,zv,gv,phasphiStfin),'AbsTol',1e-6,'RelTol',1e-6);%'AbsTol',1e-10,%,Lr,l0,q5j,dq5j,av
[tFl,xFl] = ode45(@(t,x) flight2DOut_full_red(t,x,I_bod,L0), [0 tpHM0], xStF(4:6), Opt);%Opt(t,x,I_bod,L0,g)
%[tFl,xFl] = ode4eve(@(t,x) flight2DOut_full_red(t,x,I_bod,L0), linspace(0,tpHM0,300), xStF(4:6), Opt);%Opt(t,x,I_bod,L0,g)
eulf=xFl(end,1:3);
drf=xStF(7:9)+gv*tFl(end);
nMirr=[-cos(-eulf(1)/2);sin(-eulf(1)/2);0];
RMirr=eye(3)-2*nMirr*nMirr';
xFlF=[-eulf(2);eulf(3);RMirr*drf;-RMirr*L0];%%%
tFlF=tFl(end);
phasphiTD=phasphiStfin+alr+phasphi(2)*tFl(end);%
end
function [position,isterminal,direction] = TouchDoEventsFcnStride_2(t,xF,Aznom,l0,phasphi,Lrv,alr,zv,gv,phasphiStfin)%,varargin
    rZYXR = eul2rotm1(xF);
    phasphiTD=phasphiStfin+alr+phasphi(2)*t;
    zt=zv(1)+zv(2)*t+gv(3)/2*t^2;
    rfooz=zt+rZYXR(3,:)*((roty(-Aznom)*[0;l0*sin(phasphiTD);-l0*cos(phasphiTD)])+[-Lrv(1);Lrv(2:3)]);
    position = [rfooz;(abs(cos(xF(2)))>0.01)-.5];
    isterminal = [1;1];
    direction = [-1;0];
    
end


% 
% 
% Aznom=-Aznom;%%%ANALIZEEEE
% Opt = odeset('Events', @(t,xF) TouchDoEventsFcnStride(t,xF,Aznom,l0,phasphi,Lr,alr,phasphiStfin),'AbsTol',1e-6,'RelTol',1e-6);%,Lr,l0,q5j,dq5j,av
% [tFl,xFl] = ode45(@(t,x) flight2DOut_full(t,x,gv,I_bod), [0 tpHM0], xStF, Opt);%Opt(t,x,I_bod,L0,g)
% %[tFl,xFl] = ode4eve(@(t,x) flight2DOut_full(t,x,gv,I_bod), linspace(0,tpHM0,300), xStF, Opt);%Opt(t,x,I_bod,L0,g)
% rf=xFl(end,1:3)';drf=xFl(end,7:9)';eulf=xFl(end,4:6);Lf=xFl(end,10:12)';
% rRf=eul2rotm1(eulf);
% %omegaf=rRf*I_bod^-1*rRf'*L0;%%%%%%CONSIDERRRRRRRRRRR
% %alphaf=-(q5j+dq5j*tFl(end)+av(1)*(tFl(end)-0.252227016293890)^2+av(2)*(tFl(end)-0.252227016293890)^3);
% RMirr=diag([-1,1,1]);
% RrotM=RMirr*[-rRf(1:3,1),rRf(1:3,2:3)];
% eulfM=rotm2eul(RrotM);
% phasphiTD=phasphiStfin+alr+phasphi(2)*tFl(end);%
% rFfoot=rf+rRf*((roty(-Aznom)*[0;l0*sin(phasphiTD);-l0*cos(phasphiTD)])+[0;0;Lr]);
% xFlF=[RMirr*(rf-rFfoot);eulfM';RMirr*drf;-RMirr*Lf];%%%NOIK TODO
% tFlF=tFl(end);

% 
% %
% Aznom=-Aznom;%%%ANALIZEEEE
% Opt = odeset('Events', @(t,xF) TouchDoEventsFcnStride(t,xF,Aznom,l0,phasphi,Lr,alr,phasphiStfin),'AbsTol',1e-6,'RelTol',1e-6);%,Lr,l0,q5j,dq5j,av
% [tFl,xFl] = ode45(@(t,x) flight2DOut_full(t,x,gv,I_bod), [0 tpHM0], xStF, Opt);%Opt(t,x,I_bod,L0,g)
% %[tFl,xFl] = ode4eve(@(t,x) flight2DOut_full(t,x,gv,I_bod), linspace(0,tpHM0,300), xStF, Opt);%Opt(t,x,I_bod,L0,g)
% 
% drf=xFl(end,7:9)';eulf=xFl(end,4:6);Lf=xFl(end,10:12)';
% nMirr=[-cos(-eulf(1)/2);sin(-eulf(1)/2);0];
% RMirr=eye(3)-2*nMirr*nMirr';
% xFlF=[-eulf(2);eulf(3);RMirr*drf;-RMirr*Lf];%%%
% tFlF=tFl(end);
% 
% phasphiTD=phasphiStfin+alr+phasphi(2)*tFl(end);%
