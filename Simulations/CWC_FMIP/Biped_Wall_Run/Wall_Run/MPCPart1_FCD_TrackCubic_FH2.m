%clc, clear all, close all
% PROBLEM SETUP, DO NOT CHANGE
tic
Xfo=0.12;
Yfo=0.07;
Xha=0.05;
Yha=0.05;
rfo=[0;-1;0];
rha=[0;0;1.2];
rfo1=[-0.7;0;1.2];
rfo2=[-1.5;0;1.3];
rfo3=[-2;0;1.2];
rfo4g=[-2.5;-1.4;0];
rfo5g=[-3;-0.3;0];

%rha=[0;0.5;1.6];
rpoi={rfo,rha,rfo1,rfo2,rfo3,rfo4g,rfo5g};
Xcell={Xfo,Xha,Xfo,Xfo,Xfo,Xfo,Xfo};
Ycell={Yfo,Yha,Yfo,Yfo,Yfo,Yfo,Yfo};
mcont=length(rpoi);
for irp=1:mcont
    Mrcr{irp}=cross(repmat(rpoi{irp},1,3),eye(3));
end

T=3*10^-2;
g=9.81;gv=[0;0;-g];
m=100;

M03=zeros(3);
Mgcr=cross(repmat(gv,1,3),eye(3));
Apre=[M03 -m*Mgcr M03;M03 M03 eye(3);M03 M03 M03];


Bpre=[];
for irp=1:length(rpoi)
    Bprecell{irp}=[Mrcr{irp} eye(3); M03 M03; 1/m*eye(3) M03];
    Bpre=[Bpre,Bprecell{irp}];
end
nX=size(Apre,2);
nU=size(Bpre,2);
%
A_sys = eye(nX)+T*Apre;
B_sys = T*Bpre;
N = 100;

tfs=T*N;

L0 = [0;0;0];
r0 = [0;-0.5;0.6];%[-.4;-.6;1.1];%[0;-1.5;-1.5];
dr0= 0.5*[-2;2;3];
rfd = (rfo4g+rfo5g)/2+[0;0;1];

cpre=[zeros(6,1);gv];%m*[Mgcr*rfd;zeros(6,1)];
c_sys = T*cpre;

Mdr0cr=cross(repmat(dr0,1,3),eye(3));

muf=0.8;
mupf=muf/sqrt(2);
muh=0.4;
muph=muh/sqrt(2);
rpoi={rfo,rha,rfo1,rfo2,rfo3,rfo4g,rfo5g};
mupv={mupf,muph,mupf,mupf,mupf,mupf,mupf};
Xlv={Xfo,Xha,Xfo,Xfo,Xfo,Xfo,Xfo};
Ylv={Yfo,Yha,Yfo,Yfo,Yfo,Yfo,Yfo};

pRf=[eye(3)];Rf=blkdiag(pRf,pRf);
pRh=[1 0 0; 0 0 -1; 0 1 0];Rh=blkdiag(pRh,pRh);
pRf1=pRh;Rf1=blkdiag(pRf1,pRf1);
pRf2=pRh;Rf2=blkdiag(pRf2,pRf2);
pRf3=pRh;Rf3=blkdiag(pRf3,pRf3);
pRf4g=pRf;Rf4g=blkdiag(pRf4g,pRf4g);
pRf5g=pRf;Rf5g=blkdiag(pRf5g,pRf5g);
Rcell={Rf,Rh,Rf1,Rf2,Rf3,Rf4g,Rf5g};
nv{1}=pRf*[0;0;1];
nv{2}=pRh*[0;0;1];
nv{3}=pRf1*[0;0;1];
nv{4}=pRf2*[0;0;1];
nv{5}=pRf3*[0;0;1];
nv{6}=pRf4g*[0;0;1];
nv{7}=pRf5g*[0;0;1];
%pRv

for j=1:mcont
    mup=mupv{j};
    Xl=Xlv{j};
    Yl=Ylv{j};
Wcwcj{j}=[1 0 -mup 0 0 0
      -1 0 -mup 0 0 0
      0 1 -mup 0 0 0
      0 -1 -mup 0 0 0
      0 0 -1 0 0 0
      0 0 1 0 0 0
      0 0 -Yl 1 0 0 
      0 0 -Yl -1 0 0
      0 0 -Xl 0 1 0
      0 0 -Xl 0 -1 0
      Yl Xl -mup*(Xl+Yl) mup mup 1
      Yl -Xl -mup*(Xl+Yl) +mup -mup 1
      -Yl Xl -mup*(Xl+Yl) mup -mup 1
      -Yl -Xl -mup*(Xl+Yl) -mup -mup 1  
      Yl Xl -mup*(Xl+Yl) -mup -mup -1
      Yl -Xl -mup*(Xl+Yl) -mup +mup -1
      -Yl Xl -mup*(Xl+Yl) mup -mup -1
      -Yl -Xl -mup*(Xl+Yl) mup mup -1];
Wredj{j}=Wcwcj{1}*Rcell{j}';
bcwcredj{j}=zeros(size(Wredj{j},1),1);
end
bcwcredj{1}(6)=1.6*m*g;
bcwcredj{2}(6)=0.5*m*g;
bcwcredj{3}(6)=1.6*m*g;
bcwcredj{4}(6)=1.6*m*g;
bcwcredj{5}(6)=1.6*m*g;
bcwcredj{6}(6)=1.6*m*g;
bcwcredj{7}(6)=1.6*m*g;

%%
Wbcwc=[];
bbcwc=[];
for j=1:mcont
Wbcwc=blkdiag(Wbcwc,Wredj{j});
bbcwc=[bbcwc;bcwcredj{j}];
end

load('xreft.mat')
rHf=xreft(1:3,:);
drHf=xreft(4:6,:);

%MrHfcr=cross(repmat(rHf,1,3),eye(3));
%pQLHe=[eye(3) zeros(3) -m*MrHfcr];QLHe=transpose(pQLHe)*pQLHe;
al1=100000;
al2=1000;
al3=1;

%QrHe=diag([0 0 0 1 1 1 0 0 0]);
%QdrHe=diag([0 0 0 0 0 0 1 1 1]);
%frHe=[zeros(3,1);-al2*rHf;-al3*drHf];
%%

%Qmp=diag([0 0 0 0 0 0 0 0 0]);%diag(kron(ones(1,3),[1 1000 1]));%eye(nX);%
%Qm=Qmp+al1*QLHe+al2*QrHe+al3*QdrHe;
%%
Rm=0.0000001*eye(nU);

H = [kron(eye(N+1),zeros(nX)), zeros(nX*(N+1),nU*N);zeros(nX*(N+1),nU*N)', kron(eye(N),Rm)];
fqp = [zeros(nX*(N+1)+nU*N,1)];
H = sparse(H);
fqp = sparse(fqp);
for i=1:(N+1)
    tHe=T*(i-1);
    MrHfcr=cross(repmat(rHf(:,i),1,3),eye(3));
    pQLHe=[eye(3) zeros(3) -m*MrHfcr];QLHe=transpose(pQLHe)*pQLHe;
    QrHe=diag([0 0 0 1 1 1 0 0 0]);
    QdrHe=diag([0 0 0 0 0 0 1 1 1]);
    H(nX*(i-1)+(1:nX),nX*(i-1)+(1:nX)) = al1*QLHe+al2*QrHe+al3*QdrHe;%al1*QLHe+al2*QrHe+al3*QdrHe;
    frHe=[zeros(3,1);-al2*rHf(:,i);-al3*drHf(:,i)];
    fqp(nX*(i-1)+(1:nX),1) = frHe;
end

%%

%fqp = [zeros(nX*(N+1)+nU*N,1)];

% bcwcred=zeros(size(Wred,1),1); bcwcred(6)=1.2*m*g;
% b3cwcred=bcwcred+Wred*([m*gv;-m*Mgcr*rfd]);
%
nWred=size(Wbcwc,1);
Wred=Wbcwc;
b3cwcred=bbcwc;%+Wred*([m*gv;-m*Mgcr*(rfd-rpoi{1});zeros(6,1)]);


Aqppre1 = [zeros(nWred*N,nX*(N+1)) kron(eye(N),Wred)];
bqppre1 = kron(ones(N,1),b3cwcred);
Aqppre2 = [kron(eye(N+1),[zeros(1,5) 1 zeros(1,3)]) zeros(N+1,nU*N)];
bqppre2 = kron(ones(N+1,1),1.7);

Aqppre3 = -[kron(eye(N+1),[zeros(1,5) 1 zeros(1,3)]) zeros(N+1,nU*N)];
bqppre3 = kron(ones(N+1,1),-0.4);

Aqppre4 = [kron(eye(N+1),[zeros(1,4) 1 zeros(1,4)]) zeros(N+1,nU*N)];
bqppre4 = kron(ones(N+1,1),-0.1);

Aqp = [Aqppre1;Aqppre2;Aqppre3;Aqppre4];
bqp = [bqppre1;bqppre2;bqppre3;bqppre4];

B1=[kron(eye(N),A_sys) zeros(nX*N,nX) kron(eye(N),B_sys)];

BMPC1=[eye(nX) zeros(nX,nX*N+nU*N)];
BMPC2=[zeros(nX*N,nX) eye(nX*N) zeros(nX*N,nU*N)]-B1;
BMPC3=[zeros(nX,nX*N) eye(nX) zeros(nX,nU*N)];
%BMPC4=[zeros(1,nX*(N+1)+nU*N) (nv{2})'];
%BMPC3=[zeros(3,nX*N) eye(nX/2) eye(nX/2) zeros(3,nU*N) zeros(nX/2,3)];

c1 = [L0+m*cross(r0,dr0);r0;dr0];
c2 = kron(ones(N,1),c_sys);%zeros(nX*N,1);
c3 = [zeros(3,1);rfd;zeros(3,1)];%zeros(nX/2,1);


f = [zeros(nX*(N+1)+nU*N,1)];
B = [BMPC1;BMPC2;BMPC3];
c = [c1;c2;c3];

options = optimset('Display','Off');
%[z,fval,exitflag,outp] = linprog(f,[Aqp],[bqp],B,c,[],[],[],options);
toc
tic
[z,fval,exitflag] = quadprog(H,fqp,Aqp,bqp,B,c,[],[],[],options);
toc
exitflag
%%
% SOLVE QP
options = optimset('Display','Off');
    %[z,fval,exitflag] = quadprog(H,f,Aqp,bqp,B,c,[],[],[],options);

rVt=[];
rccpI=[];
eFlag=[];
Voutp={};
Z={};   
nVjsplit=10;
% for thvj=linspace(-muph,muph,nVjsplit)%linspace(0,2*pi,nVjsplit)
% for phivj=linspace(-muph,muph,nVjsplit)%linspace(0,pi,nVjsplit)
%     tic
%     vj=[sin(phivj)*cos(thvj);sin(phivj)*sin(thvj);cos(phivj)];
%     %vj=[cos(thvj);0;sin(thvj)];
%     f = [zeros(nX*(N+1)+nU*N+3,1);1];
%     vj=[thvj;-1;phivj];
%     BMPC4=[zeros(1,nX*(N+1)+nU*N) vj' 0];
%     BMPC5=[zeros(3,nX-3) -eye(3) zeros(3,nX*(N-1)) zeros(3,nX-3) eye(3) zeros(3,nU*N) zeros(3) -vj];
%     
%     B = [BMPC1;BMPC2;BMPC3;BMPC4;BMPC5];
%     c = [c1;c2;c3;0;zeros(3,1)];
%     
%     lambbd = [zeros(nX*(N+1)+nU*N+3,1);1];
%     
%     [z,fval,exitflag,outp] = linprog(f,[Aqp],[bqp],B,c,[],[],[],options);
%     z(end)
%     %exitflag
% %     while(exitflag~=1 && exitflag~=-5)
% %         disp('Warning: Exitflag not 1')
% %         disp(exitflag)
% %         phivj=phivj+(pi/(2*nVjsplit))*(2*rand-1);thvj=thvj+(2*pi/(2*nVjsplit))*(2*rand-1);
% %     vj=[sin(phivj)*cos(thvj);sin(phivj)*sin(thvj);cos(phivj)];
% %     %vj=[cos(thvj);0;sin(thvj)];
% %     f = [zeros(nX*(N+1)+nU*N+3,1);1];
% %     
% %     BMPC4=[zeros(3,nX*(N+1)+nU*N) eye(3) -vj];
% %     B = [BMPC1;BMPC2;BMPC3;BMPC4];
% %     c = [c1;c2;c3;zeros(3,1)];
% %     
% %     [z,fval,exitflag,outp] = linprog(f,Aqp,bqp,B,c,[],[],[],options);
% %     end
%     Voutp{end+1}=outp;
%     Z{end+1}=z;
%     eFlag(end+1)=exitflag;
%     xtraj=z(1:nX*(N+1));
%     xtraj=reshape(xtraj,[nX N+1]);
%     pIfe2=xtraj(nX-2:nX,N+1)-xtraj(nX-2:nX,1);
%     %nv{2}'*pIfe2;
%     rccpI(:,end+1)=z((end-3):(end-1));
%     disp(rccpI(:,end)'*pIfe2)
%     rVt(:,end+1)=cross(nv{2},rccpI(:,end))/(nv{2}'*pIfe2);
%     toc
%     figure(1)
%     plot3(vj(1),vj(2),vj(3),'o'),hold on, grid on, drawnow
%     figure(8)
%     plot3(z(end-2),z(end-1),z(end),'o'),hold on, grid on, drawnow
% end
% end
%%
% figure(8)
% rVt1ef=find(eFlag>0);
% Fplot3(rccpI(1,rVt1ef),rccpI(2,rVt1ef),rccpI(3,rVt1ef),'ob'), hold on, grid on
% 
% figure(3);
% rVt1ef=find(eFlag>0);
% plot3(rVt(1,rVt1ef),rVt(2,rVt1ef),rVt(3,rVt1ef),'ob'), hold on, grid on
% rVtn1ef=find(eFlag<0);
% plot3(rVt(1,rVtn1ef),rVt(2,rVtn1ef),rVt(3,rVtn1ef),'or')
% rVtx=rVt(1,:);
% rVty=rVt(2,:);
% rVtz=rVt(3,:);
% axis([-10 10 -10 10 -10 10])
% kk=convhull(rVtx,rVty,rVtz);
% trisurf(kk,rVtx,rVty,rVtz);
%%
bcwcred=bbcwc;
for Vn=1:1%size(Z,2)
%z=Z{Vn};
u=z(nX*(N+1)+1:end);
u=reshape(u,[nU N]);
xtraj=z(1:nX*(N+1));
xtraj=reshape(xtraj,[nX N+1]);
%
zeta=u;%as%-[m*gv;-m*Mgcr*(rfd-rpoi{1});zeros(6,1)]*ones(1,N);


for j=1:mcont
taun=(sum(zeta((1:3)+6*(j-1),:).*zeta((4:6)+6*(j-1),:)))./(nv{j}'*zeta((1:3)+6*(j-1),:));
tau2P=zeta((4:6)+6*(j-1),:)-nv{j}*taun;
for i=1:N
    if norm(zeta((1:3)+6*(j-1),i))>0.001
    lambP(i)=-(cross(zeta((1:3)+6*(j-1),i),tau2P((1:3),i))/norm(zeta((1:3)+6*(j-1),i))^2)'*nv{j}/(nv{j}'*zeta((1:3)+6*(j-1),i));
    rP{j}(:,i)=rpoi{j}+cross(zeta((1:3)+6*(j-1),i),tau2P((1:3),i))/norm(zeta((1:3)+6*(j-1),i))^2+lambP(i)*zeta((1:3)+6*(j-1),i);
    else
    rP{j}(:,i)=rpoi{j};%[0;0;0]; 
    end
end
end

%
% %
Wgl=kron(eye(N),Wred);
Hol=Wgl*zeta(:)<kron(ones(N,1),bcwcred);
if isempty(find(Hol==0))
    disp('Inside CWC - Succed!')
else
    disp('Not inside CWC')
    max(Wgl*zeta(:)-kron(ones(N,1),bcwcred))
end
% 
% [~,Sm]=lqr(Apre,Bpre,Qm,Rm);
% K=lqr(Apre,Bpre,Sm,0.001*diag([1 1 1 1 1 1]));

x(:,1)=[L0;r0;dr0];
x(:,1)=x(:,1);%+0.03*x(:,1).*(2*rand(nX,1)-1);
% x(:,1)=[
%          0
%          0
%          0
%    -0.4049
%    -0.5956
%     1.1282
%     2.0349
%     2.0750
%    -1.0025];
Tottorq=[];
Totf=[];
for i=1:N
%     xm=[x(1:3,i)-m*cross(x(4:6,i),x(7:9,i));x(4:6,i);x(7:9,i)]-[zeros(3,1);rfd;zeros(3,1)];
%     v4lamfeed=(Wred*u(:,i)-b3cwcred)./(Wred*K*(xm-xtraj(:,i)));
%     lambfeed(i)=min(min(v4lamfeed(find(v4lamfeed>0))),1)*1;%1;%
%     uc(:,i)=u(:,i)-lambfeed(i)*K*(xm-xtraj(:,i));%lambfeed(i)*
%     zetac(:,i)=uc(:,i)-[m*gv;-m*Mgcr*(rfd-rpoi{1})];
        
    zetac(:,i)=zeta(:,i);
    totf=0;
    for jm=1:mcont
        totf=totf+zetac((1:3)+6*(jm-1),i);
    end
    ddx4p(:,i)=1/m*(totf)+gv;
    tottorq=0;
    for jm=1:mcont
        tottorq=tottorq+zetac((4:6)+6*(jm-1),i)-cross(x(4:6,i)-rpoi{jm},zetac((1:3)+6*(jm-1),i));
    end
    Tottorq(:,i)=tottorq;
    Totf(:,i)=totf;
    
    x(:,i+1)=x(:,i)+T*[tottorq;x(7:9,i);ddx4p(:,i)];%A_sys*x(:,i)+B_sys*u(:,i);
%     if max(Wred*zetac(:,i)>bcwcred)==1
%         disp('error')
%         break
%     end
end
xreft=x(4:end,:);
save('xreft','xreft')

rCMP=x(4:6,1:N)-(ones(3,1)*((nv{1}'*x(4:6,1:N))./(nv{1}'*zeta(1:3,:)))).*zeta(1:3,:);

%plot(x')
figure(1)
subplot(311)
plot(T*(1:size(x,2)),x(1:3,:)')
legend('Lx','Ly','Lz')
subplot(312)
plot(T*(1:size(x,2)),x(4:6,:)')
legend('x','y','z')
subplot(313)
plot(T*(1:size(x,2)),x(7:9,:)')
legend('$\dot{x}$','$\dot{y}$','$\dot{z}$')
%
figure(2)
subplot(131)
title('Total Force'),hold on, grid on
plot(T*(1:size(u,2)),Totf)
legend('$f_x$','$f_y$','$f_z$')
subplot(132)
title('Total Torque around CoM'),hold on, grid on
plot(T*(1:size(u,2)),Tottorq)
legend('$\tau_x$','$\tau_y$','$\tau_z$')
subplot(133)
title('CoM Accelerations'),hold on, grid on
plot(T*(1:size(u,2)),ddx4p(1:3,:)')
legend('$\ddot{x}$','$\ddot{y}$','$\ddot{z}$')
%%
figure(3)
plot3(x(4,:),x(5,:),x(6,:)); hold on, grid on
plot3(x(4,1),x(5,1),x(6,1),'xr')
plot3(x(4,end),x(5,end),x(6,end),'ok')
for j=1:mcont
R4pha=pRcell{j}*[[Xcell{j} Xcell{j} -Xcell{j} -Xcell{j} Xcell{j}];[Ycell{j} -Ycell{j} -Ycell{j} Ycell{j} Ycell{j}];zeros(1,5)]+rpoi{j}*ones(1,5);
plot3(R4pha(1,:),R4pha(2,:),R4pha(3,:),'r')
end
axis([-3 3 -3 3 0 3])

tHe=[0:T:T*N];
plot3(rHf(1,:),rHf(2,:),rHf(3,:),'k'); hold on, grid on
%%
%pause
for j=1:2
plot3(rP{j}(1,:),rP{j}(2,:),rP{j}(3,:),'g');
end
%plot3(rCMP(1,:),rCMP(2,:),rCMP(3,:),'r');
%pause
% for i=1:N
% plot3(x(4,i),x(5,i),x(6,i),'bx');
% plot3(rP(1,i),rP(2,i),rP(3,i),'bx');
% plot3(rCMP(1,i),rCMP(2,i),rCMP(3,i),'ro');
% %pause
% end
%plot3(rCMP(1,:),rCMP(2,:),rCMP(3,:),'r');
xlabel('x')
ylabel('y')
zlabel('z')
%
figure(4)
subplot(311)
title('CoM'),hold on, grid on
plot(T*(1:size(x,2)),x(4:6,:)')

tHe=[0:T:T*N];
plot(T*(1:size(x,2)),rHf(1:3,:)')
legend('x','y','z','xref','yref','zref')

subplot(312)
title('CoP Fo'),hold on, grid on
plot(T*(1:size(rP{1},2)),rP{1}(1:3,:)'-(rpoi{1}*ones(1,size(rP{1},2)'))');
legend('rCoPx','rCoPy','rCoPz')
subplot(313)
title('CoP Ha'),hold on, grid on
plot(T*(1:size(rP{2},2)),rP{2}(1:3,:)'-(rpoi{2}*ones(1,size(rP{2},2)'))');
legend('rCoPx','rCoPy','rCoPz')
%
figure(5)
subplot(131)
title('Total Force Foot'),hold on, grid on
plot(T*(1:size(u,2)),(zeta(1:3,:))')
legend('$f_x$','$f_y$','$f_z$')
subplot(132)
title('Total Force Hand'),hold on, grid on
plot(T*(1:size(u,2)),(zeta(7:9,:))')
legend('$f_x$','$f_y$','$f_z$')
% title('CMP'),hold on, grid on
% plot(T*(1:size(rP,2)),rCMP(1:3,:)')
% legend('rCMPx','rCMPy','rCMPz')
% figure(5) 
% title('$\lambda$ Feedback'),hold on, grid on
% plot(T*(1:size(lambfeed,2)),lambfeed)
% legend('$\lambda$')
figure(3)
drawnow
end