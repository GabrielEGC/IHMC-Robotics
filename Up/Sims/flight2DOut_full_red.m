function dxdt = flight2DOut_full_red(t,x,I_bod,L0)
    th=x(2);phi=x(3);%psi=eul(1);
    rZYXR = eul2rotm1(x);
    DmEAR=[ 0,           sin(phi)/cos(th),           cos(phi)/cos(th)
             0,                   cos(phi),                  -sin(phi)
         1, (sin(phi)*sin(th))/cos(th), (cos(phi)*sin(th))/cos(th)];
    dxdt=DmEAR*(I_bod^-1)*rZYXR'*L0;
end