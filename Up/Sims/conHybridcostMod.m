function [c,z]=conHybridcostMod(x,wrV)
    [z,~]=HybridcostMod(x,wrV);
    c=[];
end