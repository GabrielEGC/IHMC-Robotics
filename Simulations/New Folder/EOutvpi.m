function [value,isterminal,direction]=EOutvpi(T,X,emax)
kp=300;kd=20;
value=[((X(3)+kp*X(1)+kd*X(2))-1000*eps),(abs(X(1))-emax)];%double((abs(X(1))<emax)||((X(3)+kp*X(1)+kd*X(2))>0));%
isterminal=[1,1];
direction=[-1,0];