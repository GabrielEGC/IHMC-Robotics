function R=rotx(psi)
c=cos(psi);
s=sin(psi);
R=[1,0,0;0,c,-s;0,s,c];
end