function v=FunCostFloq(x,wrV)
    [z,xSt0]=HybridcostMod(x,wrV);
    xSt0=xSt0(4:12);
    thre=1e-4;
    z=norm(z);
    delx=1e-5;%VERIFY delx
    try
        [v,~]=Floqfun(xSt0,x,delx)
    catch
        v=1000+10*z;
    end
    %end
end