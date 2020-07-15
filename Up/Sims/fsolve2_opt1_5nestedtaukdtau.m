%clear all, close all, clc
MaxIter=1000;
TolFro=1e-8;
TolXs=1e-8;
optionsf = optimoptions('fmincon','TolFun',TolFro,'TolX',TolXs,...
'FiniteDifferenceStepSize',eps^(1/3),...
'ConstraintTolerance',1e-4,...
'MaxIter',MaxIter,'Algorithm','sqp');%FunctionTolerance
%

% xopsw=[0.257283988650583
%                    0
%    0.174532925199433
%   -0.021696561186663
%   -0.084791505694651
%   -0.015113662721195
%    0.700000000000000
%    2.000000000000000
%   -1.041434123794297
%   -0.031073905904270
%   -0.119143808945155
%   -0.000000000000001];
xopsw=[0.190799045177082
  -1.268252036378003
   0.162225732664051
   4.105104487425344
   0.716349462229682
   1.591511967400386
   2.582756515959443
   0.021921239208429
   0.089795071640447
   0.015412500130244
  -0.679383424283145
   3.241831928767350
  -0.099472447202356
  -0.023484920832579
   1.172623446508131
   0.013286108178605];
computevcceq2(xopsw(4:5))
%%
%xopsw=BestSol.Position;
%xopsw2 = fmincon(@(x) FunCostFull(x,1e+4), xopsw,[],[]);%optionsf
%xopsw2 = fmincon(@(x) FunCostFloq([x;xopsw(3:end)],1e+10), xopsw(1:2),[],[],[],[],[],[],@(x) HybridcostMod([x;xopsw(3:end)],1e+04),optionsf);%
%xopsw2 = fmincon(@(x) FunCostFloq(x,wrV), xopsw,[],[],[],[],[],[],@(x) conHybridcostMod(x,1e+04),optionsf);%
%xopsw2 = fmincon(@(x) FunCostFull([x;xopsw(3:end)],1e+10), xopsw(1:2),[],[],[],[],[],[],[],optionsf);%
[xopsw2,f,eflag,outpt]=runobjcons(xopsw(4:5),optionsf);

function [x,f,eflag,outpt] = runobjcons(x0,opts)
if nargin == 1
    opts = [];
end
xLast = [];
myf = [];
myc = [];
myceq = [];
fun = @objfun;
cfun = @constr;
[x,f,eflag,outpt] = fmincon(fun,x0,[],[],[],[],[],[],cfun,opts);
f
    function y = objfun(x)
        if ~isequal(x,xLast)
            [myf,myc,myceq] = computevcceq2(x);
            xLast = x;
        end
        y = myf;
    end
    function [c,ceq] = constr(x)
        if ~isequal(x,xLast)
            [myf,myc,myceq] = computevcceq2(x);
            xLast = x;
        end
        c = myc;
        ceq = myceq;
    end
end
%%
% xopsw=[0.374733655961109
%   -1.592314196854576
%   pi/18
%   0
%   0
%   1
%   1
%    0.015835068853870
%    0.146192185591264
%   -0.006630385844639
%   -0.174039117146840
%    2.490903341255652
%   -0.508557900387606
%   -0.032291828701717
%    0.422540236355706
%   -0.000003968984767];INTO
% xopsw2=[0.190799045177082
%   -1.268252036378003
%    0.162225732664051
%   -0.496314111856229
%   -0.131940856628591
%    1.591511967400386
%    2.582756515959443
%    0.021921239208429
%    0.089795071640447
%    0.015412500130244
%   -0.679383424283145
%    3.241831928767350
%   -0.099472447202356
%   -0.023484920832579
%    1.172623446508131
%    0.013286108178605]