function [tFlF,xFlF] = FlightPhaseOutNoStride(xStF,Aznom,l0,phasphi,Lrv,gv,I_bod,tpHM0)%,varargin
zv=xStF([3 9]);
L0=xStF(10:12);
Aznom=-Aznom;
Opt = odeset('Events', @(t,xF) TouchDoEventsFcnStride_red(t,xF,Aznom,l0,phasphi,Lrv,zv,gv),'AbsTol',1e-6,'RelTol',1e-6);%'AbsTol',1e-10,%,Lr,l0,q5j,dq5j,av
[tFl,xFl] = ode45(@(t,x) flight2DOut_full_red(t,x,I_bod,L0), [0 tpHM0], xStF(4:6), Opt);%Opt(t,x,I_bod,L0,g)
%[tFl,xFl] = ode4eve(@(t,x) flight2DOut_full_red(t,x,I_bod,L0), linspace(0,tpHM0,300), xStF(4:6), Opt);%Opt(t,x,I_bod,L0,g)
%[tFl,eulFlsol] = ode45(@(t,x) flight2DOut_full_red(t,x,I_bod,sol1.y(10:12,end)), [0 tFlpsol], sol1.y(4:6,end));%Opt(t,x,I_bod,L0,g)%(tanh(parave(2))+1)
eulf=xFl(end,1:3);
drf=xStF(7:9)+gv*tFl(end);
nMirr=[-cos(-eulf(1)/2);sin(-eulf(1)/2);0];
RMirr=eye(3)-2*nMirr*nMirr';
xFlF=[-eulf(2);eulf(3);RMirr*drf;-RMirr*L0];%%%
tFlF=tFl(end);
end
function [position,isterminal,direction] = TouchDoEventsFcnStride_red(t,xF,Aznom,l0,phasphi,Lrv,zv,gv,varargin)
    rZYXR = eul2rotm1(xF);
    phasphiTD=phasphi(1);
    zt=zv(1)+zv(2)*t+gv(3)/2*t^2;
    rfooz=zt+rZYXR(3,:)*((roty(-Aznom)*[0;l0*sin(phasphiTD);-l0*cos(phasphiTD)])+[-Lrv(1);Lrv(2:3)]);
    position = [rfooz;(abs(cos(xF(2)))>0.01)-.5];
    isterminal = [1;1];
    direction = [-1;0];
end
% rRf=eul2rotm1(eulf);
% RMirrloc=diag([-1,1,1]);
% nMirr=[-cos(-.5);sin(-.5);0];
% RMirr=eye(3)-2*nMirr*nMirr';
% RrotM=RMirr*rRf*RMirrloc;%transposed?
% eulfM=rotm2eul(RrotM);
% xFlF=[eulfM';RMirr*drf;-RMirr*Lf];%%%NOIK TODO
% tFlF=tFl(end);
% drf=xFl(end,7:9)';eulf=xFl(end,4:6);Lf=xFl(end,10:12)';
% nMirr=[-cos(-0);sin(-0);0];
% RMirr=eye(3)-2*nMirr*nMirr';
% xFlF=[2*0-eulf(1);-eulf(2);eulf(3);RMirr*drf;-RMirr*Lf];%%%NOIK TODO
% tFlF=tFl(end);



% %%zv=xStF([3 6]);
% L0=xStF(10:12);
% Aznom=-Aznom;
% Opt = odeset('Events', @(t,xF) TouchDoEventsFcnStride_red(t,xF,Aznom,l0,phasphi,Lr,zv,gv),'AbsTol',1e-6,'RelTol',1e-6);%'AbsTol',1e-10,%,Lr,l0,q5j,dq5j,av
% [tFl,xFl] = ode45(@(t,x) flight2DOut_full_red(t,x,I_bod,L0), [0 tpHM0], xStF(4:6), Opt);%Opt(t,x,I_bod,L0,g)
% %[tFl,xFl] = ode4eve(@(t,x) flight2DOut_full(t,x,gv,I_bod), linspace(0,tpHM0,300), xStF, Opt);%Opt(t,x,I_bod,L0,g)
% %[tFl,eulFlsol] = ode45(@(t,x) flight2DOut_full_red(t,x,I_bod,sol1.y(10:12,end)), [0 tFlpsol], sol1.y(4:6,end));%Opt(t,x,I_bod,L0,g)%(tanh(parave(2))+1)
% eulf=xFl(end,1:3);
% drf=xStF(7:9)+gv*tFl(end);
% nMirr=[-cos(-eulf(1)/2);sin(-eulf(1)/2);0];
% RMirr=eye(3)-2*nMirr*nMirr';
% xFlF=[-eulf(2);eulf(3);RMirr*drf;-RMirr*L0];%%%
% tFlF=tFl(end);