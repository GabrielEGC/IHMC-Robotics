function [position,isterminal,direction] = LiftOffEventsFcn(t,xF,Fextfunhand,Lrv,Aznom,kpLe,kpAz,tau0,l0,I_bod,kdtau,kdLe,kdAz)
    r0=xF(1:3);eul=xF(4:6);dr0=xF(7:9);L0=xF(10:12);
    rZYXR = eul2rotm1(eul);
    omgc=rZYXR*(I_bod^-1)*rZYXR'*L0;
    rloc=rZYXR'*([0;0;0]-r0)-Lrv;
    drloc=-rZYXR'*(dr0+cross(omgc,([0;0;0]-r0)));
    FextNoRot = Fextfunhand(rloc,Aznom,kpLe,kpAz,tau0,l0,drloc,kdtau,kdLe,kdAz);
    FextW=rZYXR*FextNoRot;
    position = [FextW(3);(abs(cos(eul(2)))>0.01)-.5;r0(3)-0.001;rloc(3)+0.001;((rloc(1)^2+rloc(3)^2)>0.01)-.5];
    %[((l0-l)>=0||(dl<0))||(t<1e-10);yCoM]; % The value that we want to be zero[(l0-l)>0||(dl<0);yCoM]
    isterminal = [1;1;1;1;1];%[1;1];  % Halt integration 
    direction = [-1;0;0;0;0];%[0;0];   % The zero can be approached from either direction
end