syms psiaz phi r
%roty(-psiaz)*[0;0;-1]
%rotx(phi)*[0;0;-1]
%%
%rotx(phi)*roty(-psiaz)*[0;0;-1]
%%
xT=roty(-psiaz)*rotx(phi)*[0;0;-1]*r;
%%
Jrq=jacobian(xT,[psiaz;phi;r]).'
%%
syms x1 x2 x3; assume([x1 x2 x3],'real')
jacobian(atan2(x1,-x3),[x1 x2 x3])
jacobian(sqrt(sum([x1;x2;x3].^2)),[x1 x2 x3])
%%
jacobian(atan2(x2,sqrt(x3^2+x1^2)),[x1 x2 x3])
%%