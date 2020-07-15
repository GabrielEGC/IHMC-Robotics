%clear all, close all, clc
MaxFunEvals=20000;
MaxIter=1000;
TolFro=1e-8;
TolXs=1e-8;
optionsfsol = optimoptions('fsolve','TolFun',TolFro,'TolX',TolXs,'MaxFunEvals',MaxFunEvals,'MaxIter',MaxIter);%FunctionTolerance
%load('xopsw')
xopsw=[  0.257284216037292
                   0
  -0.021696386802502
  -0.084791924720442
  -0.015113518410895
   0.700000000000000
   2.000000000000000
  -1.041424973611196
  -0.031073899191838
  -0.119145459240669
   0.000000000000039];
%xopsw=BestSol.Position;
%xopsw=xopsw2;
for i=1:1
    xopsw = fsolve(@(x) FunCostFull(x,1e+04), xopsw,optionsfsol);%2*rand(12,1)-1
    FunCostFull(xopsw,1e+04)
end
%%
%clear all, close all, clc
MaxIter=100;
TolFro=1e-8;
TolXs=1e-8;
optionsf = optimoptions('fmincon','TolFun',TolFro,'TolX',TolXs,...
'FiniteDifferenceStepSize',eps^(1/3),...
'MaxIter',MaxIter,'Algorithm','sqp');%FunctionTolerance
%
%     
% xopsw=[   0.000244180253760
%    0.036083314241827
%   -0.021672570619034
%   -0.084578868626387
%   -0.015089480711711
%    0.699002558577425
%    2.000025292004221
%   -1.039978509681531
%   -0.031066546121688
%   -0.119029106666637
%    0.000000178453952
xopsw=[0.228049483343706
  -0.000000000000002
  -0.010742352622381
  -0.005508968801729
  -0.005036432134696
   0.320843564525583
   1.604932499049682
  -0.344498334445252
  -0.028934832201468
  -0.023221319389143
   0.000000000027759];
%
%xopsw=BestSol.Position;
%xopsw2 = fmincon(@(x) FunCostFull(x,1e+4), xopsw,[],[]);%optionsf
%xopsw2 = fmincon(@(x) FunCostFloq([x;xopsw(3:end)],1e+10), xopsw(1:2),[],[],[],[],[],[],@(x) HybridcostMod([x;xopsw(3:end)],1e+04),optionsf);%
xopsw2 = fmincon(@(x) FunCostFloq(x,1e+10), xopsw,[],[],[],[],[],[],@(x) conHybridcostMod(x,1e+04),optionsf);%
%xopsw2 = fmincon(@(x) FunCostFull([x;xopsw(3:end)],1e+10), xopsw(1:2),[],[],[],[],[],[],[],optionsf);%
%%
xopsw=[  -0.000023227263413
  -0.036571690879223
  -0.647922684168409
   0.000398020690377
   1.102380689849323
  -0.000437631000645
  -0.002965710950230
   0.115159516373456
  -0.000020475171806
  -2.012481819101979
  -1.205837035558581
   2.924655651762201
   4.098432110239302];