syms psi th phi L0p1 L0p2 L0p3
L0p=[L0p1;L0p2;L0p3];
I_bod=diag([2.5;1.0;0.5]);
    rZYXR = [ cos(psi)*cos(th), cos(psi)*sin(phi)*sin(th) - cos(phi)*sin(psi), sin(phi)*sin(psi) + cos(phi)*cos(psi)*sin(th)
 cos(th)*sin(psi), cos(phi)*cos(psi) + sin(phi)*sin(psi)*sin(th), cos(phi)*sin(psi)*sin(th) - cos(psi)*sin(phi)
         -sin(th),                              cos(th)*sin(phi),                              cos(phi)*cos(th)];
    DmEA=[ (cos(psi)*sin(th))/cos(th), (sin(psi)*sin(th))/cos(th), 1
                  -sin(psi),                   cos(psi), 0
           cos(psi)/cos(th),           sin(psi)/cos(th), 0];
    deuldt=simplify(DmEA*rZYXR*(I_bod^-1)*transpose(rZYXR)*L0p);
    matlabFunction(deuldt,'File','eqdeuldt');
    dfdyp=simplify(jacobian(deuldt,[psi th phi L0p1 L0p2 L0p3]))
    matlabFunction(dfdyp,'File','eqdfdyp');