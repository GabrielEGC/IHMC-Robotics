function [tFlF,xFlF,phasphiTD] = FlightPhaseOutStride_av(xStF,Aznom,l0,phasphi,Lrv,gv,I_bod,tpHM0,alr,phasphiStfin,av)%,varargin
zv=xStF([3 9]);
L0=xStF(10:12);
Aznom=-Aznom;
Opt = odeset('Events', @(t,xF) TouchDoEventsFcnStride_2(t,xF,Aznom,l0,phasphi,Lrv,alr,zv,gv,phasphiStfin,av),'AbsTol',1e-6,'RelTol',1e-6);%'AbsTol',1e-10,%,Lr,l0,q5j,dq5j,av
[tFl,xFl] = ode45(@(t,x) flight2DOut_full_red(t,x,I_bod,L0), [0 tpHM0], xStF(4:6), Opt);%Opt(t,x,I_bod,L0,g)
%[tFl,xFl] = ode4eve(@(t,x) flight2DOut_full_red(t,x,I_bod,L0), linspace(0,tpHM0,300), xStF(4:6), Opt);%Opt(t,x,I_bod,L0,g)
eulf=xFl(end,1:3);
drf=xStF(7:9)+gv*tFl(end);
nMirr=[-cos(-eulf(1)/2);sin(-eulf(1)/2);0];
RMirr=eye(3)-2*nMirr*nMirr';
xFlF=[-eulf(2);eulf(3);RMirr*drf;-RMirr*L0];%%%
tFlF=tFl(end);
phasphiTD=phasphiStfin+alr+phasphi(2)*tFl(end)+av(1)*(tFl(end)-0.010220029143547)^2+av(2)*(tFl(end)-0.010220029143547)^3;%
end
function [position,isterminal,direction] = TouchDoEventsFcnStride_2(t,xF,Aznom,l0,phasphi,Lrv,alr,zv,gv,phasphiStfin,av)%,varargin
    rZYXR = eul2rotm1(xF);
    phasphiTD=phasphiStfin+alr+phasphi(2)*t+av(1)*(t-0.010220029143547)^2+av(2)*(t-0.010220029143547)^3;
    zt=zv(1)+zv(2)*t+gv(3)/2*t^2;
    rfooz=zt+rZYXR(3,:)*((roty(-Aznom)*[0;l0*sin(phasphiTD);-l0*cos(phasphiTD)])+[-Lrv(1);Lrv(2:3)]);
    position = [rfooz;(abs(cos(xF(2)))>0.01)-.5];
    isterminal = [1;1];
    direction = [-1;0];
    
end