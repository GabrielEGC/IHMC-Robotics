%-0.050449585495789  -0.020894780187015  -0.016590704316786   0.906921314425220   0.000023195841719
clear all, close all, clc
xopsw2=[0.046334284562866
  -5.355010042230490
   0.113524936450297
  -2.035906118674400
   4.954170143877920
   1.550635410871368
   5.165420748579254
   2.435748668787137
   2.674809901488776
   0.149457073807350
   0.846646376468854
   0.743743990645261
   0.054642291958716
   0.355987019879191
  -0.237653168642405
   0.125572488975305
   0.566559485228246
   1.113299684743702
   1.995390842602884
  -0.779869923535340
   0.535150016261177
   0.549135727176047
   0.248226304619958];
    x=xopsw2;
    InitSetUp;
    [x0p,alr,tFlF] = Striderest(xSt0,x);
    xSt0 = xSt0(5:12);
    z=x0p(2:end)-xSt0;zn=norm(z)
    delx=1e-7;
    [v,A]=Floqfunopt([phasphi(1);xSt0],x,delx,x0p,alr)
Q = diag(ones(9,1))*1e+10;
P = dlyap(A',Q);P=P/norm(P);
%%
xSt0=[phasphi(1);xSt0];
%%
% Nopt2=20;
% xgnonor=2*rand(9,Nopt2)-1;V0nonor=zeros(Nopt2,1);
% for i=1:Nopt2
%     V0nonor(i)=xgnonor(:,i)'*P*xgnonor(:,i);
%     xgnonor(:,i)=xgnonor(:,i)/sqrt(V0nonor(i))*.017+xSt0;
% end

Nopt2=20;
xgnonor=2*rand(9,Nopt2)-1;V0nonor=zeros(Nopt2,1);
for i=1:Nopt2
    V0nonor(i)=xgnonor(:,i)'*P*xgnonor(:,i);
    xgnonor(:,i)=xgnonor(:,i)/sqrt(V0nonor(i))*.017+xSt0;
    xgnonor(:,i)=xSt0;xgnonor(2:3,i)=xgnonor(2:3,i)+(2*rand(2,1)-1)*0.04;
end
%%
SP=xgnonor
MaxIter=1000;TolFro=1e-12;TolXs=1e-12;
optionsf = optimoptions('fmincon','TolFun',TolFro,'TolX',TolXs,'MaxIter',MaxIter);%FunctionTolerance
%xopsw=2*rand(5,1)-1;%[0 0 0 0 0];
avopprev=[0;0];%10*rand(2,1);%[0;0];%
avop = fmincon(@(av) FunCostFull2(av,xSt0,xopsw2,P,SP,100,alr,A), avopprev,[],[],[],[],[],[],[],optionsf);%optionsf
%%
SP=xgnonor
MaxIter=1000;TolFro=1e-12;TolXs=1e-12;
optionsf = optimoptions('fminimax','TolFun',TolFro,'TolX',TolXs,'MaxIter',MaxIter);%FunctionTolerance
%xopsw=2*rand(5,1)-1;%[0 0 0 0 0];
avopprev=10*rand(2,1);%[0;0];%
avop = fminimax(@(av) FunCostFull2mm(av,xSt0,xopsw2,P,SP,100,alr,A), avopprev,[],[],[],[],[],[],[],optionsf);%optionsf