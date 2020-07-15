% RMirr=diag([-1,1,1]);
% rFfoot=rf-RMirr*(ya(1:3)-[0;0;0]);
% RMirr*(rf-rFfoot)=ya(1:3);
% eulfM=ya(4:6);
% RMirr*drf=ya(7:9);
% -RMirr*Lf=ya(10:12);
% 
% xFlF=[RMirr*(rf-rFfoot);eulfM';RMirr*drf;-RMirr*Lf];%%%NOIK TODO
% Opt = odeset('Events', @(t,xF) TouchDoEventsFcn(t,xF));%'AbsTol',1e-10,%,Lr,l0,q5j,dq5j,av
% [tFl,xFl] = ode45(@(t,x) flight2DOut_full(t,x,gv,I_bod), [0 tComFl], xStF, Opt);%Opt(t,x,I_bod,L0,g)
% rf=xFl(end,1:3)';drf=xFl(end,7:9)';eulf=xFl(end,4:6);Lf=xFl(end,10:12)';
% rRf=eul2rotm1(eulf);
% 
% RrotM=RMirr*[-rRf(1:3,1),rRf(1:3,2:3)];
% eulfM=rotm2eul(RrotM);
%%
xStF=yb(1:12);
%Opt = odeset('Events', @(t,xF) TouchDoEventsFcn(t,xF));%'AbsTol',1e-10,%,Lr,l0,q5j,dq5j,av
[tFl,xFl] = ode45(@(t,x) flight2DOut_full(t,x,gv,I_bod), [0 tFlpa], xStF);%Opt(t,x,I_bod,L0,g)
rf=xFl(end,1:3)';drf=xFl(end,7:9)';eulf=xFl(end,4:6);Lf=xFl(end,10:12)';
rRf=eul2rotm1(eulf);
%MIRRORING
RMirr=diag([-1,1,1]);
RrotM=RMirr*[-rRf(1:3,1),rRf(1:3,2:3)];
eulfM=rotm2eul(RrotM);
rFfoot=rf-RMirr*(xSt0(1:3)-[0;0;0]);
%disp('zFfoot: ');disp(rFfoot(3));
xFlF=[RMirr*(rf-rFfoot);eulfM';RMirr*drf;-RMirr*Lf];%%%NOIK TODO
tFlF=tFl(end);