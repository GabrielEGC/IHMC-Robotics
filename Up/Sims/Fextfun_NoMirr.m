function Fext = Fextfun_NoMirr(rloc,Aznom,kpLe,kpAz,tau0,l0,drloc,kdtau,kdLe,kdAz)
    FextNoRot = Fextfun(diag([-1;1;1])*rloc,Aznom,kpLe,kpAz,tau0,l0,diag([-1;1;1])*drloc,kdtau,kdLe,kdAz);%_NoMirr
    Fext=diag([-1;1;1])*FextNoRot;
end

