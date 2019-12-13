function dxdt=massforclimvpi(t,x)
ksys=39.478;
bsys=6.283;
vmax=200;
A=[0 1;-ksys -bsys];
Aa=[A,[0;ksys];zeros(1,3)];Ba=[0;0;1];
v=vmax;
dxdt=-(Aa*x+Ba*v);
end