function [position,isterminal,direction] = TouchDoEventsFcnStride(t,xF,Aznom,l0,phasphi,Lrv,alr,phasphiStfin)%,varargin
    r0=xF(1:3);
    eul=xF(4:6);
    rZYXR = eul2rotm1(eul');
    phasphiTD=phasphiStfin+alr+phasphi(2)*t;
    rfoo=r0+rZYXR*((roty(-Aznom)*[0;l0*sin(phasphiTD);-l0*cos(phasphiTD)])+Lrv);
    position = [rfoo(3);(abs(cos(eul(2)))>0.01)-.5]; % pf4(3) The value that we want to be zero
    isterminal = [1;1];  % Halt integration sssssssssssssss
    direction = [-1;0];   % The zero can be approached from either direction
end