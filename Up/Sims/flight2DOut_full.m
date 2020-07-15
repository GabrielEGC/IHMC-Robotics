function dxdt = flight2DOut_full(t,x,gv,I_bod)%I_bod,L0
    %r0=x(1:3);
    eul=x(4:6);dr0=x(7:9);L0=x(10:12);
    psi=eul(1);th=eul(2);phi=eul(3);
    rZYXR = eul2rotm1([psi th phi]);
    DmEA=[ (cos(psi)*sin(th))/cos(th), (sin(psi)*sin(th))/cos(th), 1
                  -sin(psi),                   cos(psi), 0
           cos(psi)/cos(th),           sin(psi)/cos(th), 0];
    deuldt=DmEA*rZYXR*(I_bod^-1)*rZYXR'*L0;
    dxdt=[dr0;deuldt;gv;[0;0;0]];
end