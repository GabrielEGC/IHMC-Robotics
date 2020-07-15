clear all, close all, clc
syms I1 I2 I3 L0p1 L0p2 L0p3 dL1 dL2 dL3
syms psi th phi
syms dpsi dth dphi
I_bod=diag([I1;I2;I3]);
L0=[L0p1;L0p2;L0p3];
dL=[dL1;dL2;dL3];
rZYXR = rotz(psi)*roty(th)*rotx(phi);
DmEA=[ (cos(psi)*sin(th))/cos(th), (sin(psi)*sin(th))/cos(th), 1
              -sin(psi),                   cos(psi), 0
       cos(psi)/cos(th),           sin(psi)/cos(th), 0];
eul=[psi;th;phi];
deul=[dpsi;dth;dphi];
deuleq=simplify(DmEA*rZYXR*(I_bod^-1)*transpose(rZYXR)*L0);
ddeul=jacobian(deuleq,eul)*deul+jacobian(deuleq,L0)*dL;
L0=simplify(rZYXR*I_bod*transpose(rZYXR)*(DmEA^-1)*deul);
L0p1=L0(1);
L0p2=L0(2);
L0p3=L0(3);
ddeul=simplify(subs(ddeul));
omgc=DmEA^-1*deul;%deuldt=DmEA*omgc;
Ek=1/2*simplify(transpose(omgc)*rZYXR*(I_bod)*transpose(rZYXR)*omgc);
V=0;
Lg=Ek-V;
dLgddeul=transpose(jacobian(Lg,deul));
ddLgddeuldt=simplify(jacobian(dLgddeul,eul)*deul+jacobian(dLgddeul,deul)*ddeul);
genforc=simplify(ddLgddeuldt-transpose(jacobian(Lg,eul)))

%%
clear all, close all, clc
syms psi th phi
rRf=rotz(psi)*roty(th)*rotx(phi);
RMirr=diag([-1,1,1]);
RrotM=RMirr*rRf*RMirr;
rRf=rotz(-psi)*roty(-th)*rotx(phi);
simplify(RrotM-rRf)

%%
rRf=rotz(xFl(end,1))*roty(xFl(end,2))*rotx(xFl(end,3));
%MIRRORING
RMirr=diag([-1,1,1]);
RrotM=RMirr*rRf*RMirr;%[-rRf(1:3,1),rRf(1:3,2:3)]
eulfM=[xFl(end,1:2),-xFl(end,3)];%rotz(-psi)*roty(-th)*rotx(phi);