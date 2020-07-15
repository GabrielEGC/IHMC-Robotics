function v=FunCostFull(x,wrV)
    %x(1)=x(1)*1e-4;
    thre=1e-4;
    [z,xSt0]=HybridcostMod(x,wrV);z=norm(z);
    xSt0=xSt0(4:12);
    v=20*z/thre;
    delx=1e-5;%VERIFY delx
    if z<thre
        [v,~]=Floqfun(xSt0,x,delx);v=v+.05*z/thre
    end
end