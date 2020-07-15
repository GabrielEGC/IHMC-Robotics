function R=roty(theta)
theta=-theta;
c=cos(theta);
s=sin(theta);
R=[c,0,-s;0,1,0;s,0,c];
end