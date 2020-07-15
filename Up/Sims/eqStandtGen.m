InitSetParam;
%omgc4t=rZYXR*(I_bod^-1)*transpose(rZYXR)*[L01;L02;L03];
%rloc=transpose(rZYXR)*([0;0;0]-[xr01;yr02;zr03])-[0;0;Lr];
%drloc=-transpose(rZYXR)*([vx0;vy0;vz0]+cross(omgc4t,([0;0;0]-[xr01;yr02;zr03])));
%dphasest4tau=(rloc(3)*drloc(2)-drloc(3)*rloc(2))/(rloc(3)^2+rloc(2)^2);
%tau0=kdtau*dphasest4tau+x(14)^2;

%initrvec=[0;0;0]-rZYXR*((roty(-Aznom)*[0;l0*sin(phasphi);-l0*cos(phasphi)])+[0;0;Lr]);
%xr01=initrvec(1);yr02=initrvec(2);zr03=initrvec(3);%0.604293541844519
%xSt0=[xr01;yr02;zr03;eua1;eua2;eua3;vx0;vy0;vz0;L01;L02;L03];

%%CONSIDER INITIALIZATION IN SWING TRAJECTORY

syms r0x r0y r0z psi th phi dr0x dr0y dr0z L0p1 L0p2 L0p3
%(t,x,I_bod,Lr,Aznom,kpLe,kpAz,tau0,l0,gv,m,kdtau)
    r0=[r0x;r0y;r0z];
    dr0=[dr0x;dr0y;dr0z];
    L0=[L0p1;L0p2;L0p3];
    rZYXR = rotz(psi)*roty(th)*rotx(phi);
    DmEA=[ (cos(psi)*sin(th))/cos(th), (sin(psi)*sin(th))/cos(th), 1
                  -sin(psi),                   cos(psi), 0
           cos(psi)/cos(th),           sin(psi)/cos(th), 0];
    omgc=rZYXR*(I_bod^-1)*transpose(rZYXR)*L0;
    deuldt=DmEA*omgc;
    
    rloc=transpose(rZYXR)*([0;0;0]-r0)-[0;0;Lr];
    drloc=-transpose(rZYXR)*(dr0+cross(omgc,([0;0;0]-r0)));
    
    azhat=[-rloc(3);0;rloc(1)];azhat=azhat/(norm(azhat)^2);%%FUCKKKKKKKKK
    ttauhat=[0;-rloc(3);rloc(2)];ttauhat=ttauhat/norm(ttauhat);
    
    dphasest=(rloc(3)*drloc(2)-drloc(3)*rloc(2))/(rloc(3)^2+rloc(2)^2);
    Fl=(-kpLe*(norm(rloc)-l0))*rloc/norm(rloc);
    FAz=-kpAz*(atan2(rloc(1),-rloc(3))-Aznom)*azhat;
    Ftau=-(tau0-kdtau*dphasest)/norm(rloc)*roty(-Aznom)*ttauhat;%%-kdLe*dnormrloc
    FextNoRot=-(Fl+FAz+Ftau);
    FextW=simplify(rZYXR*FextNoRot);
    dxdt=[dr0;deuldt;FextW/m+gv;cross([0;0;0]-r0,FextW)];
    dxdt=simplify(dxdt);
    TEn=1/2*m*(transpose(dr0)*dr0)+1/2*(transpose(omgc)*rZYXR*(I_bod)*transpose(rZYXR)*omgc)+m*(-gv(3))*r0(3)+1/2*kpLe*(norm(rloc)-l0)^2+1/2*kpAz*(atan2(rloc(1),-rloc(3))-Aznom)^2;
    matlabFunction(dxdt,'File','eqdSRBMdt');
    matlabFunction(FextW(3),'File','eqFextWz');
    matlabFunction(TEn,'File','eqEnergy');
    %%
    dfdy=jacobian(dxdt,[r0x r0y r0z psi th phi dr0x dr0y dr0z L0p1 L0p2 L0p3]);
    matlabFunction(dfdy,'File','eqdfdy');%%SLOW
    
    %%
    

    deuldt=simplify(DmEA*rZYXR*(I_bod^-1)*transpose(rZYXR)*L0p);
    matlabFunction(deuldt,'File','eqdeuldt');
    dfdyp=simplify(jacobian(deuldt,[psi th phi L0p1 L0p2 L0p3]))
    matlabFunction(dfdyp,'File','eqdfdyp');