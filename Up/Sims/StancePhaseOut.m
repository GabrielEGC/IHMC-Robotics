function [tStF,xStF] = StancePhaseOut(xSt0,I_bod,Lrv,Aznom,kpLe,kpAz,tau0,l0,gv,m,tpHM0,kdtau,kdLe,kdAz)
Opt = odeset('Events', @(t,xF) LiftOffEventsFcn(t,xF,@Fextfun,Lrv,Aznom,kpLe,kpAz,tau0,l0,I_bod,kdtau,kdLe,kdAz),'AbsTol',1e-6,'RelTol',1e-6);%,'AbsTol',1e-10,'RelTol',1e-10
[tSt,xSt] = ode45(@(t,x) stance2DOut_full(t,x,I_bod,Lrv,Aznom,kpLe,kpAz,tau0,l0,gv,m,kdtau,kdLe,kdAz), [0 tpHM0], xSt0, Opt);%
%[tSt,xSt] = ode4eve(@(t,x) stance2DOut_full(t,x,I_bod,Lr,Aznom,kpLe,kpAz,tau0,l0,gv,m,kdtau,kdLe,kdAz), linspace(0,tpHM0,300), xSt0, Opt);%
xStF=xSt(end,:)';
tStF=tSt(end);
end