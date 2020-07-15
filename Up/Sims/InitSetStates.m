    rZYXR=eul2rotm1([0;x(16:17)]);
    initrvec=[0;0;0]-rZYXR*((roty(-Aznom)*[0;l0*sin(phasphi(1));-l0*cos(phasphi(1))])+Lrv);
    xSt0=[initrvec;0;x(16:end)];