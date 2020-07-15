run('InitSetParam.m')
run('InitSetParamVar.m')
rZYXR=eul2rotm1([0;xIn(2:3)]);
initrvec=[0;0;0]-rZYXR*((roty(-Aznom)*[0;l0*sin(xIn(1));-l0*cos(xIn(1))])+Lrv);
xSt0=[initrvec;0;xIn(2:end)];