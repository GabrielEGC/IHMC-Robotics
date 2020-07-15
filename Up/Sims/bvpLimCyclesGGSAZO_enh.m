function sol1=bvpLimCyclesGGSAZO(x,solver)
if nargin < 2
   solver = 'bvp4c';
end
bvpsolver = fcnchk(solver);
tpHM0=[];phasphi=[];I_bod=[];gv=[];Lrv=[];m=[];l0=[];Aznom=[];kpLe=[];kpAz=[];kdtau=[];tau0=[];g=[];kdLe=[];kdAz=[];
InitSetParam;
InitSetParamVar;
xindf=0.15;%scaled5
parave=[1];
%cstjacBC={[eye(3);zeros(3)],[zeros(3);eye(3)],zeros(6,3)};
%optsbvp = bvpset('FJacobian',@jacsysD,'BCJacobian',cstjacBC,'Vectorized','on','Stats','on');
optsbvp = bvpset('FJacobian',@jacsysD,'Vectorized','on','RelTol',1e-8,'AbsTol',1e-8,'Stats','on');
solinit1 = bvpinit(linspace(0,xindf,10),@sysDinit1,parave);
tic
sol1 = bvpsolver(@sysD,@sysDres,solinit1,optsbvp);
toc
%xint = linspace(0,xindf);
%solingue1=sysDinit1(xint);
tFlpsol=(sol1.y(9,1)-sol1.y(9,end))/(gv(3));
[tFl,eulFlsol] = ode45(@(t,x) flight2DOut_full_red(t,x,I_bod,sol1.y(10:12,end)), [0 tFlpsol], sol1.y(4:6,end));%Opt(t,x,I_bod,L0,g)%(tanh(parave(2))+1)
eulFlsol=eulFlsol';
rfFlsol=sol1.y(1:3,end)*ones(size(tFl'))+sol1.y(7:9,end)*tFl'+gv*(tFl'.^2)/2;
drfFlsol=sol1.y(7:9,end)*ones(size(tFl'))+gv*tFl';
%tFl=tFl+xindf*sol1.parameters(2);

figure
title('Solution/Angles');
subplot(221); hold on; grid on
plot(sol1.x,sol1.y(1,:),'o-',sol1.x,sol1.y(2,:),'o-',sol1.x,sol1.y(3,:),'o-');
plot(tFl,rfFlsol(1,:),'o-',tFl,rfFlsol(2,:),'o-',tFl,rfFlsol(3,:),'o-')
legend('$x$','$y$','$z$')
xlabel('x');
ylabel('Position');
subplot(222); hold on; grid on
plot(sol1.x,sol1.y(4,:),'o-',sol1.x,sol1.y(5,:),'o-',sol1.x,sol1.y(6,:),'o-');
plot(tFl,eulFlsol(1,:),'o-',tFl,eulFlsol(2,:),'o-',tFl,eulFlsol(3,:),'o-')
legend('$\psi$','$\theta$','$\phi$')
xlabel('x');
ylabel('Euler Angles');
subplot(223); hold on; grid on
plot(sol1.x,sol1.y(7,:),'o-',sol1.x,sol1.y(8,:),'o-',sol1.x,sol1.y(9,:),'o-');
legend('$\dot{x}$','$\dot{y}$','$\dot{z}$')
xlabel('x');
ylabel('Velocity');
subplot(224); hold on; grid on
plot(sol1.x,sol1.y(10,:),'o-',sol1.x,sol1.y(11,:),'o-',sol1.x,sol1.y(12,:),'o-');
legend('$L_x$','$L_y$','$L_z$')
xlabel('x');
ylabel('Angular Momentum');

figure;
solr0x=sol1.y(1,:);solr0y=sol1.y(2,:);solr0z=sol1.y(3,:);
solpsi=sol1.y(4,:);solth=sol1.y(5,:);solphi=sol1.y(6,:);
title('Force and Energy');
subplot(211); hold on; grid on
plot(sol1.x,eqFextWz(solphi,solpsi,solr0x,solr0y,solr0z,solth),'o-');
legend('$F_z$')
xlabel('x');
ylabel('Vertical Force');
subplot(212); hold on; grid on
plot(sol1.x,eqEnergy(sol1.y(10,:),sol1.y(11,:),sol1.y(12,:),sol1.y(7,:),sol1.y(8,:),sol1.y(9,:),solphi,solpsi,solr0x,solr0y,solr0z,solth),'o-');
legend('$E_T$')
xlabel('x');
ylabel('Total Energy');
%fprintf('Angular Momentum required %7.3f.\n',sol1.parameters)
%Sxint1 = deval(sol1,xint);
%plot(xint,Sxint1(1,:),'o-',xint,solingue1(1,:));

function dydx = sysD(x,y,parave)%VERIFY L0P
    r0x=y(1,:);r0y=y(2,:);r0z=y(3,:);
    psi=y(4,:);th=y(5,:);phi=y(6,:);
    dr0x=y(7,:);dr0y=y(8,:);dr0z=y(9,:);
    L0p1=y(10,:);L0p2=y(11,:);L0p3=y(12,:);   
    dydx=parave(1)*eqdSRBMdt(L0p1,L0p2,L0p3,dr0x,dr0y,dr0z,phi,psi,r0x,r0y,r0z,th);
end
function [dfdy,dfdp] = jacsysD(x,y,parave)% p, [dfdy,dfdp] = fjac(x,y,p)
    r0x=y(1,:);r0y=y(2,:);r0z=y(3,:);
    psi=y(4,:);th=y(5,:);phi=y(6,:);
    dr0x=y(7,:);dr0y=y(8,:);dr0z=y(9,:);
    L0p1=y(10,:);L0p2=y(11,:);L0p3=y(12,:);
    dfdy=parave(1)*eqdfdy(L0p1,L0p2,L0p3,phi,psi,r0x,r0y,r0z,th);
    dfdp=[eqdSRBMdt(L0p1,L0p2,L0p3,dr0x,dr0y,dr0z,phi,psi,r0x,r0y,r0z,th)];
end
function res = sysDres(ya,yb,parave) %%%%%%NESTEDDDDDDDD
    xStF=yb(4:6);
    L0=yb(10:12);
    tFlp=(ya(9)-yb(9))/(gv(3));
    %Opt = odeset('Events', @(t,xF) TouchDoEventsFcn(t,xF));%'AbsTol',1e-10,%,Lr,l0,q5j,dq5j,av
    Opt = odeset('RelTol',1e-10,'AbsTol',1e-10);%'AbsTol',1e-10,%,Lr,l0,q5j,dq5j,av
    if tFlp==0
        xFl=xStF';
    else
        [~,xFl] = ode45(@(t,x) flight2DOut_full_red(t,x,I_bod,L0), [0 tFlp], xStF,Opt);%Opt(t,x,I_bod,L0,g)%(tanh(parave(2))+1)
    end
    RMirr=diag([-1,1,1]);
    xFlF312=[yb(3)+tFlp*yb(9)+gv(3)*(tFlp^2)/2;-xFl(end,1:2)';xFl(end,3);-yb(7);yb(8);-RMirr*L0];
    rZYXRa=rotz(ya(4))*roty(ya(5))*rotx(ya(6));rloca=rZYXRa'*([0;0;0]-ya(1:3))-[0;0;Lr];
    phasvpo=roty(Aznom)*rloca;%TOO SENSIBLEE
    atan2(phasvpo(2),-phasvpo(3));
    res=[rloca'*rloca-l0^2;%ya(9)+0.5;%
        ya(7)-.7;%%
        ya(8)-2;%atan2(phasvpo(2),-phasvpo(3))-0.1;
        ya([3:8 10:12]')-xFlF312;
        eqFextWz(yb(6),yb(4),yb(1),yb(2),yb(3),yb(5))];%%Consider signed initial conditions with parameters
end
function yinit = sysDinit1(x)%%%%DOESN'T CONTAIN PARAMETERS
    yinit=[-.05+1*(-xindf/4+x-x^2/xindf);3*(x-xindf/2);0.4-0.6*(x-x^2/xindf);
        zeros(3,1);
        1*(1-2*x/xindf);3;-0.6*(1-2*x/xindf);
        zeros(3,1)];%TEST FOR atan2(phasvpo(2),-phasvpo(3))-.26;TOO SENSIBLEE
end
end
%xindf+.5*(xindf/4)
% function [dBCdya,dBCdyb,dBCdp] = jacBC(ya,yb,p) %%[dbcdya,dbcdyb,dbcdp] = bcjac(ya,yb,p)%%%REPLACE BY CELL CONSTANT JACOBIAN
%     dBCdya=[eye(3);zeros(3)];
%     dBCdyb=[zeros(3);eye(3)];
%     dBCdp=[zeros(6,3)];
% end