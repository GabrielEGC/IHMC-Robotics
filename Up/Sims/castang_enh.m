xopsw=[   0.174547472865521
  -2.989793406358777
   0.026837652588568
 -21.528213393064714
   3.514985407896324
  10.077073427530705
  40.532435962135409
   9.905351026048002
   5.600030561417103
   0.149571359363582
   0.668710097692177
   0.712013054788482
   0.251196211298017
   0.460794340023554
  -0.275514285450344
   0.164045834692033
   0.590760461674338
   1.267428077143814
   2.390099225927656
  -0.761383583297375
   0.588464103617993
   0.572493751016180
   0.301142771652694
];
computevcceq(xopsw)
ansbvp=bvpLimCyclesGGSAZO_enh(xopsw,'bvp4c');


%%
ansbvp=bvpLimCyclesGGSAZO_enh('bvp4c');
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
