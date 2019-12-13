XallsolfM=cell2mat(Xallsolf);
v=(Xfi(:,2:end)-Xfi(:,1:end-1))./sqrt(ones(3,1)*sum((Xfi(:,2:end)-Xfi(:,1:end-1)).^2));
vave=sum(v'); vave=vave/norm(vave);

v=(Xfmin(:,2:end)-Xfmin(:,1:end-1))./sqrt(ones(3,1)*sum((Xfmin(:,2:end)-Xfmin(:,1:end-1)).^2));
vave=sum(v,2); vave=vave/norm(vave);
npla=[Km 1]'/norm([Km 1]);
nplac=cross(vave,npla);
%
NPP=[npla nplac];
rpp=XallsolfM'*NPP;
[ABP]=fmincon(@(x)x'*[0 1;1 0]*x,[0;0],[-rpp.^2],-ones(size(rpp,1),1));
tploell=linspace(0,2*pi,20);xplell=1/sqrt(ABP(1))*cos(tploell);yplell=1/sqrt(ABP(2))*sin(tploell)
figure(4)
plot(rpp(:,1),rpp(:,2)), hold on
plot(xplell,yplell);
%%
syms tplocil t2

xplcil=1/sqrt(ABP(1))*cos(tplocil);
yplcil=1/sqrt(ABP(2))*sin(tplocil);
rpp2=NPP*[xplcil;yplcil]+vave*t2;
figure(2)
epmax=emax;
epmax=2;
EZS=ezsurf(rpp2(1),rpp2(2),rpp2(3),[-epmax/vave(1)*0.7,epmax/vave(1)*0.7]);EZS.EdgeColor='none';
%%

%%
RMNF=[npla nplac vave];
P=RMNF*diag([ABP;0])*RMNF';