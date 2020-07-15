function Fext = Fextfun(rloc,Aznom,kpLe,kpAz,tau0,l0,drloc,kdtau,kdLe,kdAz)
    azhat=[-rloc(3);0;rloc(1)];azhat=azhat/norm(azhat)^2;
    ttauhat=[-rloc(1)*rloc(2); rloc(1)^2 + rloc(3)^2; -rloc(2)*rloc(3)];ttauhat=ttauhat/(norm(ttauhat)*norm(rloc));%
    
    dphasest=ttauhat'*drloc;
    dnormrloc=rloc'/norm(rloc)*drloc;
    dazrloc=azhat'*drloc;
    
    Fl=(-kpLe*(norm(rloc)-l0)-kdLe*dnormrloc)*rloc/norm(rloc);
    FAz=(-kpAz*(atan2(rloc(1),-rloc(3))-Aznom)-kdAz*dazrloc)*azhat;
    Ftau=(-kdtau*dphasest+tau0(1))*ttauhat;%%
    Fext=-(Fl+FAz+Ftau);
end
% Fl=(-kpLe*(norm(rloc)-l0)-kdLe*dnormrloc+tau0(3))*rloc/norm(rloc);
%     FAz=(-kpAz*(atan2(rloc(1),-rloc(3))-Aznom)-kdAz*dazrloc+tau0(2))*azhat;
%     Ftau=(-kdtau*dphasest+tau0(1))*ttauhat;%%