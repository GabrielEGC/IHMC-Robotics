w%clear all, close all, clc
MaxFunEvals=20000;
MaxIter=1000;
TolFro=1e-8;
TolXs=1e-8;
optionsfsol = optimoptions('fsolve','TolFun',TolFro,'TolX',TolXs,'MaxFunEvals',MaxFunEvals,'MaxIter',MaxIter);%FunctionTolerance
%load('xopsw')
xopsw=[   0
                   0
  -0.021696575977682
  -0.084791460158090
  -0.015113675183019
   0.699999999999999
   2.000000000000000
  -1.041435004762093
  -0.031073903450261
  -0.119143642856341
   0.000000000002166];
%xopsw=BestSol.Position;
%xopsw=xopsw2;
for i=1:1
    xopsw = fsolve(@(x) FunCostFull(x,1e+04), xopsw,optionsfsol);%2*rand(12,1)-1
    FunCostFull(xopsw,1e+04)
end
%%
%clear all, close all, clc
MaxIter=20;
TolFro=1e-8;
TolXs=1e-8;
optionsf = optimoptions('fmincon','TolFun',TolFro,'TolX',TolXs,...
    'FiniteDifferenceStepSize',eps^(1/8),...
    'MaxIter',MaxIter);%FunctionTolerance
xopsw=[   0.000244180253760
   0.036083314241827
  -0.021672570619034
  -0.084578868626387
  -0.015089480711711
   0.699002558577425
   2.000025292004221
  -1.039978509681531
  -0.031066546121688
  -0.119029106666637
   0.000000178453952];%ode4
%   -0.021696575977682
%   -0.084791460158090
%   -0.015113675183019
%    0.699999999999999
%    2.000000000000000
%   -1.041435004762093
%   -0.031073903450261
%   -0.119143642856341
%    0.000000000002166];%%%ODE45
%xopsw=BestSol.Position;
%xopsw2 = fmincon(@(x) FunCostFull(x,1e+4), xopsw,[],[]);%optionsf
xopsw2 = fmincon(@(x) FunCostFull(x,1e+4), xopsw,[],[],[],[],[],[],[],optionsf);%optionsf
%xopsw2 = fmincon(@(x) FunCostFull([x;xopsw(3:end)],1e+10), xopsw(1:2),[],[],[],[],[],[],[],optionsf);%
%%   
   
xopsw2=0.036083314241827;
xev=xopsw2+linspace(-4*eps^(1/1.5),4*eps^(1/1.5),20);yev=xev;
for i=1:length(xev)
yev(i)=feval(@(x) FunCostFull([0.000244180253760...
    ;x;xopsw(3:end)],1e+10),xev(i));
end
figure;hold on; grid on
plot(xev,yev)
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