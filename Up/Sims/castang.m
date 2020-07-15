ansbvp=bvpLimCyclesGGSAZO('bvp4c');
ya=ansbvp.y(:,1);
run('InitSetParam.m')
% rZYXRa=rotz(ya(4))*roty(ya(5))*rotx(ya(6));
% rloca=rZYXRa'*([0;0;0]-ya(1:3))-[0;0;Lr];
% phasvpo=roty(Aznom)*rloca;%TOO SENSIBLEE
% phasphi=atan2(phasvpo(2),-phasvpo(3))
phasphi=castangphi(ya(1:6),Lr);
xopsw=[phasphi;0;Aznom;ansbvp.y(4:end,1)]
chckinco(xopsw,ya);

function chckinco(x,ya)
InitSetUp;
if norm(xSt0(1:3)-ya(1:3))>1e-3
    xSt0(1:3)-ya(1:3)
    error('Mismatch in casted initial conditions');
end
end
