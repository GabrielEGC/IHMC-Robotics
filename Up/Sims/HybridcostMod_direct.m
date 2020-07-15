function [z,xSt0]=HybridcostMod_direct(x,wrV)
    I_bod=[];gv=[];Lr=[];m=[];l0=[];Aznom=[];kpLe=[];kpAz=[];kdtau=[];tau0=[];g=[];
    tpHM0=[];
    InitSetUp_direct;
    %xSt0
    [Fz0,~,~] = LiftOffEventsFcn(0,xSt0,@Fextfun,Lr,Aznom,kpLe,kpAz,tau0,l0,I_bod,kdtau);
    if (norm(xSt0)>1e+12||Fz0(1)<=-1e-5)
        z=wrV;disp('ret. wrV1');return;
    end
    %CAREFULL WITH xStFS
    [tStF,xStFS] = StancePhaseOut(xSt0,I_bod,Lr,Aznom,kpLe,kpAz,tau0,l0,gv,m,tpHM0,kdtau);
    [Fz0,~,~]=LiftOffEventsFcn(tStF,xStFS,@Fextfun,Lr,Aznom,kpLe,kpAz,tau0,l0,I_bod,kdtau);
    if ((tStF==tpHM0)||(tStF<0.001)||xStFS(3)<0||abs(Fz0(1))>1e-6)
        z=wrV;disp('ret. wrV2');return;
    end
    %tStF
    %xStFS
    z=sysDres(xSt0,xStFS,0);
    z=z(1:12);%z=(xSt0-xFlF);%disp(norm(z));
   
function res = sysDres(ya,yb,parave) %%%%%%NESTEDDDDDDDD
    xStF=yb(4:6);
    L0=yb(10:12);
    tFlp=(ya(9)-yb(9))/(gv(3));
    %Opt = odeset('Events', @(t,xF) TouchDoEventsFcn(t,xF));%'AbsTol',1e-10,%,Lr,l0,q5j,dq5j,av
    if tFlp==0
        xFl=xStF';
    else
        [~,xFl] = ode45(@(t,x) flight2DOut_full_red(t,x,I_bod,L0), [0 tFlp], xStF);%Opt(t,x,I_bod,L0,g)%(tanh(parave(2))+1)
    end
    RMirr=diag([-1,1,1]);
    xFlF312=[yb(3)+tFlp*yb(9)+gv(3)*(tFlp^2)/2;-xFl(end,1:2)';xFl(end,3);-yb(7);yb(8);-RMirr*L0];
    rZYXRa=rotz(ya(4))*roty(ya(5))*rotx(ya(6));rloca=rZYXRa'*([0;0;0]-ya(1:3))-[0;0;Lr];
    phasvpo=roty(Aznom)*rloca;%TOO SENSIBLEE
    atan2(phasvpo(2),-phasvpo(3))
    res=[rloca'*rloca-l0^2;%ya(9)+0.5;%
        ya(7)-.7;%%
        ya(8)-8;%atan2(phasvpo(2),-phasvpo(3))-0.1;
        ya([3:8 10:12]')-xFlF312;
        eqFextWz(yb(6),yb(4),yb(1),yb(2),yb(3),yb(5))];%%Consider signed initial conditions with parameters
end
    
end
function dydx = sysD(x,y,parave)%VERIFY L0P
    r0x=y(1,:);r0y=y(2,:);r0z=y(3,:);
    psi=y(4,:);th=y(5,:);phi=y(6,:);
    dr0x=y(7,:);dr0y=y(8,:);dr0z=y(9,:);
    L0p1=y(10,:);L0p2=y(11,:);L0p3=y(12,:);   
    dydx=parave(1)*eqdSRBMdt(L0p1,L0p2,L0p3,dr0x,dr0y,dr0z,phi,psi,r0x,r0y,r0z,th);
end