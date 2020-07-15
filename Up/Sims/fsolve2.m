clear all, close all, clc
MaxFunEvals=20000;
MaxIter=1000;
TolFro=1e-8;
TolXs=1e-8;
optionsfsol = optimoptions('fsolve','TolFun',TolFro,'TolX',TolXs,'MaxFunEvals',MaxFunEvals,'MaxIter',MaxIter);%FunctionTolerance
%load('xopsw')
xopsw=[   0.257284216037292
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
% for i=1:1
%     xopsw = fsolve(@(x) HybridcostMod([0;0;x],1e+04), xopsw(3:end),optionsfsol);%2*rand(12,1)-1
%     HybridcostMod([0;0;xopsw],1e+04)
% end
%%
for i=1:1
    xopsw = fsolve(@(x) HybridcostMod(x,1e+04), rand(11,1),optionsfsol);%2*rand(12,1)-1
    HybridcostMod(xopsw,1e+04)
end

%HybridcostMod_direct
%%
%clear all, close all, clc
MaxIter=200;
TolFro=1e-8;
TolXs=1e-8;
optionsf = optimoptions('fmincon','TolFun',TolFro,'TolX',TolXs,...
    'FiniteDifferenceStepSize',eps^(1/2),...
    'MaxIter',MaxIter,'Algorithm','sqp');%FunctionTolerance
%xopsw=[-0.050482374690964  -0.021001156905058  -0.016625156401357   0.874223157956022   0.002198387270519];
%xopsw=2*rand(5,1)-1;%[0 0 0 0 0];
xopsw=[   0.257284216037292
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
xopsw2 = fmincon(@(x) nHybridcostMod(x,1e+4), xopsw,[],[],[],[],[],[],[],optionsf);%optionsf

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