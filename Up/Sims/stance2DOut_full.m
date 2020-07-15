function dxdt = stance2DOut_full(t,x,I_bod,Lrv,Aznom,kpLe,kpAz,tau0,l0,gv,m,kdtau,kdLe,kdAz)%kpSe,kdSe,q60
	%g=-gv(3);%q=x(1:6);%dq=x(7:12);
    r0=x(1:3);eul=x(4:6);dr0=x(7:9);L0=x(10:12);
    psi=eul(1);th=eul(2);phi=eul(3);
    rZYXR = eul2rotm1([psi th phi]);
    DmEA=[ (cos(psi)*sin(th))/cos(th), (sin(psi)*sin(th))/cos(th), 1
                  -sin(psi),                   cos(psi), 0
           cos(psi)/cos(th),           sin(psi)/cos(th), 0];
    omgc=rZYXR*(I_bod^-1)*rZYXR'*L0;
    deuldt=DmEA*omgc;
    rloc=rZYXR'*([0;0;0]-r0)-Lrv;
    drloc=-rZYXR'*(dr0+cross(omgc,([0;0;0]-r0)));
    FextNoRot = Fextfun(rloc,Aznom,kpLe,kpAz,tau0,l0,drloc,kdtau,kdLe,kdAz);
    FextW=rZYXR*FextNoRot;
    dxdt=[dr0;deuldt;FextW/m+gv;cross([0;0;0]-r0,FextW)];
    %c[x(3),FextW(3),norm(rloc)]
end