function z2=nHybridcostMod(x,wrV)
    [z,~]=HybridcostMod(x,wrV);
    z2=norm(z);
end