function ph=castangphi(xStF,Lrv)
    rZYXR=eul2rotm1(xStF(4:6));
    rloc=rZYXR'*([0;0;0]-xStF(1:3))-Lrv;
    ph=pi/2-acos(dot(rloc,[0;1;0])/norm(rloc));
end