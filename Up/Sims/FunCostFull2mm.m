function v=FunCostFull2(av,xSt0,xopsw2,P,SP,wrV,alr,A)
%     thre=3e-6;
%     [z,xSt0]=HybridcostMod(x,wrV);
%     xSt0=xSt0(2:5);
%     v=10*z/thre;delx=0.0001;
%     if z<thre
%         [v,~]=Floqfun(xSt0,x,delx);v=v+5*z/thre
%     end
%PxSt0=Stride_av(av,xSt0, wrV,xopsw2,alr);
Vn=@(x)(x-xSt0)'*P*(x-xSt0);%Vnp=@(x)(x-x0p)'*P*(x-x0p);
%Vn=@(x)(x-(PxSt0+A*(x-xSt0)))'*P*(x-(PxSt0+A*(x-xSt0)));%Vnp=@(x)(x-x0p)'*P*(x-x0p);
Vnp1=@(x)Vn(Stride_av(av,x, wrV,xopsw2,alr));%Stride_av(av,xIn,wrV,x,alr)
vCost=zeros(size(SP,2),1);
for i=1:size(SP,2)
    vCost(i)=Vnp1(SP(:,i));
end
v=(vCost).^(1/4)
end