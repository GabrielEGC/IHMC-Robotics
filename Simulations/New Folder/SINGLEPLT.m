Opt=odeset('RelTol',1e-10,'AbsTol',1e-10,'Events',@(t,x)EOutvpi(t,x,emax));
[t,X]=ode45(@(t,x) massforclimi(t,x),[0 0.7],[-0.355040165614578;11.456473750172634;9.780665567173706],Opt);
Xa=X';
figure(2)
plot3(Xa(1,:),Xa(2,:),Xa(3,:)), hold on, grid on
plot3(Xa(1,1),Xa(2,1),Xa(3,1),'or')
plot3(Xa(1,end),Xa(2,end),Xa(3,end),'ok')