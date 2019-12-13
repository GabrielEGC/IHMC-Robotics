clc, clear all, close all
% PROBLEM SETUP, DO NOT CHANGE
tic
Xfo=0.12;
Yfo=0.08;
rfo=[0.44;2.06;0];
rfo1=[0;1.4;1];
rfo2=[0;0.6;1.6];
rfo3=[0.9;0;1.75];
rfo4=[1.7;0;1.4];
rfo5g=[2;0.8;0];%[2.2;0.6;0];
rfo6g=[2.4;0.4;0];%[2.2;0.6;0];

%rfo6g=[0;];

%rha=[0;0.5;1.6];
rpoi={rfo,rfo1,rfo2,rfo3,rfo4,rfo5g,rfo6g};
Xcell={Xfo,Xfo,Xfo,Xfo,Xfo,Xfo,Xfo};
Ycell={Yfo,Yfo,Yfo,Yfo,Yfo,Yfo,Yfo};
mcont=length(rpoi);
for irp=1:mcont
    Mrcr{irp}=cross(repmat(rpoi{irp},1,3),eye(3));
end


N = 100;

ActCont={[1];[2];[3];[4];[5];[6];[6;7]};
Nsw=round(linspace(0,N,length(ActCont)+1));

%Nsw = [0    2    11    20   29    38    47    56    65    74    75   100];
%Nsw = [0    15  25  33  38   43     49   55  60  65  75    100];
%Nsw = [0     30    40    45    50    55    60    66    74    82    91   100];
Nsw = [0    30  40 50 60 70  80   100];
T=2.5*10^-2;
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


for delta=1:length(ActCont)
    Bpredelta{delta}=[];
    for irp=1:length(rpoi)
        if isempty(find(ActCont{delta}==irp))
            Bprecell{irp}=0*[Mrcr{irp} eye(3); M03 M03; 1/m*eye(3) M03];
        else
            Bprecell{irp}=[Mrcr{irp} eye(3); M03 M03; 1/m*eye(3) M03];
        end
    Bpredelta{delta}=[Bpredelta{delta},Bprecell{irp}];
    end
end
%%
nX=size(Apre,2);
nU=size(Bpre,2);
%
A_sys = eye(nX)+T*Apre;
B_sys = T*Bpre;


tfs=T*N;

L0 = [0;0;0];
r0 = rfo+[0;0;0.7];%[0.4;2.06;1.2];
dr0 = 1.5*[-0.2;-0.6;0.5];%[0;0;0];%8*[-0.2;-0.6;0];%0.5
rfd = (rfo5g+rfo6g)/2+[0;0;1.2];
    
cpre=[zeros(6,1);gv];%m*[Mgcr*rfd;zeros(6,1)];
c_sys = T*cpre;

Mdr0cr=cross(repmat(dr0,1,3),eye(3));

muf=1;
mupf=muf/sqrt(2);
muh=0.8;
muph=muh/sqrt(2);
rpoi={rfo,rfo1,rfo2,rfo3,rfo4,rfo5g,rfo6g};
mupv={mupf,mupf,mupf,mupf,mupf,mupf,mupf};
Xlv={Xfo,Xfo,Xfo,Xfo,Xfo,Xfo,Xfo};
Ylv={Yfo,Yfo,Yfo,Yfo,Yfo,Yfo,Yfo};

pRf=rotz(-90);Rf=blkdiag(pRf,pRf);
pRf1=rotx(-45)*roty(90)*rotz(-90);Rf1=blkdiag(pRf1,pRf1);
pRf2=roty(90)*rotz(-90);Rf2=blkdiag(pRf2,pRf2);
pRf3=rotx(-90);Rf3=blkdiag(pRf3,pRf3);
pRf4=roty(45)*rotx(-90);Rf4=blkdiag(pRf4,pRf4);
pRf5g=rotz(45);Rf5g=blkdiag(pRf5g,pRf5g);
pRf6g=rotz(55);Rf6g=blkdiag(pRf6g,pRf6g);
Rcell={Rf,Rf1,Rf2,Rf3,Rf4,Rf5g,Rf6g};
pRcell={pRf,pRf1,pRf2,pRf3,pRf4,pRf5g,pRf6g};
nv{1}=pRf*[0;0;1];
nv{2}=pRf1*[0;0;1];
nv{3}=pRf2*[0;0;1];
nv{4}=pRf3*[0;0;1];
nv{5}=pRf4*[0;0;1];
nv{6}=pRf5g*[0;0;1];
nv{7}=pRf6g*[0;0;1];
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
Wredj{j}=Wcwcj{j}*Rcell{j}';
bcwcredj{j}=zeros(size(Wredj{j},1),1);
end
bcwcredj{1}(6)=1.5*m*g;
bcwcredj{2}(6)=1.5*m*g;
bcwcredj{3}(6)=1.5*m*g;
bcwcredj{4}(6)=1.5*m*g;
bcwcredj{5}(6)=1.5*m*g;
bcwcredj{6}(6)=1.5*m*g;
bcwcredj{7}(6)=1.5*m*g;


%%
Wbcwc=[];
bbcwc=[];
for j=1:mcont
Wbcwc=blkdiag(Wbcwc,Wredj{j});
bbcwc=[bbcwc;bcwcredj{j}];
end

%%

syms tHe;
Msolcoef=[tfs^2*eye(3) tfs^3*eye(3);2*tfs*eye(3) 3*tfs^2*eye(3)]\[rfd-r0-dr0*tfs;-dr0];
coea3=Msolcoef(1:3);coea4=Msolcoef(4:6);clear Msolcoef;

%rHf=(r0-rfd)*(1-tHe/(T*N))+rfd;
rHf=r0+dr0*tHe+coea3*tHe^2+coea4*tHe^3;
drHf=diff(rHf,tHe);
MrHfcr=cross(repmat(rHf,1,3),eye(3));
MdrHfcr=cross(repmat(drHf,1,3),eye(3));
pQLHe=[eye(3) m*MdrHfcr -m*MrHfcr];QLHe=transpose(pQLHe)*pQLHe;

vQLHe=-m*cross(rHf,drHf);

qfHe=transpose(pQLHe)*vQLHe;

al1=100000;
al2=1000;
al3=1;
QrHe=diag([0 0 0 1 1 1 0 0 0]);
QdrHe=diag([0 0 0 0 0 0 1 1 1]);
frHe=[zeros(3,1);-al2*rHf;-al3*drHf]-al1*qfHe;
%%

Qmp=diag([0 0 0 0 0 0 0 0 0]);%diag(kron(ones(1,3),[1 1000 1]));%eye(nX);%
Qm=Qmp+al1*QLHe+al2*QrHe+al3*QdrHe;
%%
Rm=0.0000001*eye(nU);
toc
tic
H = [kron(eye(N+1),zeros(size(Qm))), zeros(nX*(N+1),nU*(N+1));zeros(nX*(N+1),nU*(N+1))', kron(eye(N+1),Rm)];
fqp = [zeros(nX*(N+1)+nU*(N+1),1)];
H = sparse(H);
fqp = sparse(fqp);
for i=1:(N+1)
    tHe=T*(i-1);
    H(nX*(i-1)+(1:nX),nX*(i-1)+(1:nX)) = eval(Qm);
    fqp(nX*(i-1)+(1:nX),1) = eval(frHe);
end
toc
tic

nWred=size(Wbcwc,1);
Wred=Wbcwc;
b3cwcred=bbcwc;

Aqppre1 = [zeros(nWred*(N+1),nX*(N+1)) kron(eye(N+1),Wred)];
bqppre1 = kron(ones(N+1,1),b3cwcred);
Aqppre2 = [kron(eye(N+1),[zeros(1,5) 1 zeros(1,3)]) zeros(N+1,nU*(N+1))];
bqppre2 = kron(ones(N+1,1),2.4);

Aqppre3 = -[kron(eye(N+1),[zeros(1,5) 1 zeros(1,3)]) zeros(N+1,nU*(N+1))];
bqppre3 = kron(ones(N+1,1),-0.4);

Aqppre4 = -[kron(eye(N+1),[zeros(1,4) 1 zeros(1,4)]) zeros(N+1,nU*(N+1))];
bqppre4 = kron(ones(N+1,1),-0.1);

Aqppre5 = -[kron(eye(N+1),[zeros(1,3) 1 zeros(1,5)]) zeros(N+1,nU*(N+1))];
bqppre5 = kron(ones(N+1,1),-0.1);

Aqp = [Aqppre1;Aqppre2;Aqppre3;Aqppre4;Aqppre5];
bqp = [bqppre1;bqppre2;bqppre3;bqppre4;bqppre5];

BB1=[];
for delta=1:length(ActCont)
   BB1=blkdiag(BB1,kron(eye(Nsw(delta+1)-Nsw(delta)),T*Bpredelta{delta}));
end
B1=[kron(eye(N),A_sys) zeros(nX*N,nX) BB1 zeros(nX*N,nU)];
B1f=[zeros(nX,nX*N) A_sys zeros(nX,nU*N) T*Bpredelta{end}];

BMPC1=[eye(nX) zeros(nX,nX*N+nU*(N+1))];
BMPC2=[zeros(nX*N,nX) eye(nX*N) zeros(nX*N,nU*(N+1))]-B1;
BMPC3=[zeros(nX,nX*N) eye(nX) zeros(nX,nU*(N+1))]-B1f;
%BMPC4=[zeros(nX,nX*N) eye(nX) zeros(nX,nU*N)];

c1 = [L0+m*cross(r0,dr0);r0;dr0];
c2 = kron(ones(N,1),c_sys);
c3 = kron(ones(1,1),c_sys);

%c4 = [zeros(3,1);rfd;zeros(3,1)];

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
rVt=[];
rccpI=[];
eFlag=[];
Voutp={};
Z={};

bcwcred=bbcwc;
u=z(nX*(N+1)+1:end);
u=reshape(u,[nU N+1]);
xtraj=z(1:nX*(N+1));
xtraj=reshape(xtraj,[nX N+1]);
%
zeta=u;

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

% %
Wgl=kron(eye(N+1),Wred);
Hol=Wgl*zeta(:)<kron(ones(N+1,1),bcwcred);
if isempty(find(Hol==0))
    disp('Inside CWC - Succed!')
else
    disp('Not inside CWC')
    max(Wgl*zeta(:)-kron(ones(N+1,1),bcwcred))
end
% 
% [~,Sm]=lqr(Apre,Bpre,Qm,Rm);
% K=lqr(Apre,Bpre,Sm,0.001*diag([1 1 1 1 1 1]));

x(:,1)=[L0;r0;dr0];
x(:,1)=x(:,1);%+0.03*x(:,1).*(2*rand(nX,1)-1);

Tottorq=[];
Totf=[];
for i=1:N
    deltai=min(find((i<=Nsw)))-1;
    zetac(:,i)=zeta(:,i);
    totf=0;
    for jm=ActCont{deltai}'
        totf=totf+zetac((1:3)+6*(jm-1),i);
    end
    ddx4p(:,i)=1/m*(totf)+gv;
    tottorq=0;
    for jm=ActCont{deltai}'
        tottorq=tottorq+zetac((4:6)+6*(jm-1),i)-cross(x(4:6,i)-rpoi{jm},zetac((1:3)+6*(jm-1),i));
    end
    Tottorq(:,i)=tottorq;
    Totf(:,i)=totf;
    x(:,i+1)=x(:,i)+T*[tottorq;x(7:9,i);ddx4p(:,i)];
end
xreft=x(4:end,:);
save('xreft','xreft')
%rCMP=x(4:6,1:N)-(ones(3,1)*((nv{1}'*x(4:6,1:N))./(nv{1}'*zeta(1:3,:)))).*zeta(1:3,:);
%
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
plot(T*(1:(size(u,2)-1)),Totf)
legend('$f_x$','$f_y$','$f_z$')
subplot(132)
title('Total Torque around CoM'),hold on, grid on
plot(T*(1:(size(u,2)-1)),Tottorq)
legend('$\tau_x$','$\tau_y$','$\tau_z$')
subplot(133)
title('CoM Accelerations'),hold on, grid on
plot(T*(1:(size(u,2)-1)),ddx4p(1:3,:)')
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

axis([0 2.5 0 2.5 0 2.5])
view(150,30)
tHe=[0:T:T*N];
% rHf1=rHf(1);rHf2=rHf(2);rHf3=rHf(3);
% rHfrefp1=eval(rHf1);
% rHfrefp2=eval(rHf2);
% rHfrefp3=eval(rHf3);
% plot3(rHfrefp1,rHfrefp2,rHfrefp3,'k'); hold on, grid on
%%
%pause
for j=1:mcont
plot3(rP{j}(1,:),rP{j}(2,:),rP{j}(3,:),'g'), hold on, grid on 
%quiver3(rP{j}(1,:),rP{j}(2,:),rP{j}(3,:),zeta(1+6*(j-1),:),zeta(2+6*(j-1),:),zeta(3+6*(j-1),:),100)
end
xlabel('x')
ylabel('y')
zlabel('z')
%
figure(4)
subplot(311)
title('CoM'),hold on, grid on
plot(T*(1:size(x,2)),x(4:6,:)')

tHe=[0:T:T*N];
rHfrefp=eval(rHf);
plot(T*(1:size(x,2)),rHfrefp(1:3,:)')
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
figure(3)

% for isimpl=1:N
%     deltais=min(find((i<=Nsw)))-1;
%     plot3(x(4,isimpl),x(5,isimpl),x(6,isimpl),'x'); hold on, grid on
%     for j=ActCont{deltais}
%     %plot3(R4pha(1,:),R4pha(2,:),R4pha(3,:),'r')
%     %plot3(rP{j}(1,:),rP{j}(2,:),rP{j}(3,:),'g'), hold on, grid on 
%     quivhand=quiver3(rP{j}(1,isimpl),rP{j}(2,isimpl),rP{j}(3,isimpl),zeta(1+6*(j-1),isimpl),zeta(2+6*(j-1),isimpl),zeta(3+6*(j-1),isimpl),0.0005);
%     end
%     pause(0.1)    
% end
drawnow


%%Iterations


for itttt=1:40
itttt
pause(0.001)

load('xreft.mat')
rHf=xreft(1:3,:);
drHf=xreft(4:6,:);

al1=100000;
al2=1000;
al3=1;


Rm=0.001*eye(nU);




H = [kron(eye(N+1),zeros(nX)), zeros(nX*(N+1),nU*(N+1));zeros(nX*(N+1),nU*(N+1))', kron(eye(N+1),Rm)];
fqp = [zeros(nX*(N+1)+nU*(N+1),1)];
H = sparse(H);
fqp = sparse(fqp);
for i=1:(N+1)
    MdrHfcr=cross(repmat(drHf(:,i),1,3),eye(3));
    MrHfcr=cross(repmat(rHf(:,i),1,3),eye(3));
    pQLHe=[eye(3) m*MdrHfcr -m*MrHfcr];QLHe=transpose(pQLHe)*pQLHe;
    QrHe=diag([0 0 0 1 1 1 0 0 0]);
    QdrHe=diag([0 0 0 0 0 0 1 1 1]);
    
    vQLHe=-m*cross(rHf(:,i),drHf(:,i));
    qfHe=transpose(pQLHe)*vQLHe;
    
    H(nX*(i-1)+(1:nX),nX*(i-1)+(1:nX)) = al1*QLHe+al2*QrHe+al3*QdrHe;%al1*QLHe+al2*QrHe+al3*QdrHe;
    frHe=[zeros(3,1);-al2*rHf(:,i);-al3*drHf(:,i)]-al1*qfHe;
    fqp(nX*(i-1)+(1:nX),1) = frHe;
end

[z,fval,exitflag] = quadprog(H,fqp,Aqp,bqp,B,c,[],[],[],options);
toc
exitflag
%%

rVt=[];
rccpI=[];
eFlag=[];
Voutp={};
Z={};
%%
bcwcred=bbcwc;

u=z(nX*(N+1)+1:end);
u=reshape(u,[nU N+1]);
xtraj=z(1:nX*(N+1));
xtraj=reshape(xtraj,[nX N+1]);
%
zeta=u;


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

Wgl=kron(eye(N+1),Wred);
Hol=Wgl*zeta(:)<kron(ones(N+1,1),bcwcred);
if isempty(find(Hol==0))
    disp('Inside CWC - Succed!')
else
    disp('Not inside CWC')
    max(Wgl*zeta(:)-kron(ones(N+1,1),bcwcred))
end

x(:,1)=[L0+m*cross(r0,dr0);r0;dr0];
x(:,1)=x(:,1);%+0.03*x(:,1).*(2*rand(nX,1)-1);
Tottorq=[];
Totf=[];
for i=1:N
    deltai=min(find((i<=Nsw)))-1;
    zetac(:,i)=zeta(:,i);
    totf=0;
    for jm=ActCont{deltai}'
        totf=totf+zetac((1:3)+6*(jm-1),i);
    end
    
    ddx4p(:,i)=1/m*(totf)+gv;
    tottorq=m*cross(x(4:6,i),gv);
    for jm=ActCont{deltai}'
        %tottorq=tottorq+zetac((4:6)+6*(jm-1),i)-cross(x(4:6,i)-rpoi{jm},zetac((1:3)+6*(jm-1),i));
        tottorq=tottorq+zetac((4:6)+6*(jm-1),i)+cross(rpoi{jm},zetac((1:3)+6*(jm-1),i));
    end
    Tottorq(:,i)=tottorq;
    Totf(:,i)=totf;
    x(:,i+1)=x(:,i)+T*[tottorq;x(7:9,i);ddx4p(:,i)];
    LCoMest(:,i)=x(1:3,i)-m*cross(x(4:6,i),x(7:9,i));
end
xreft=x(4:end,:);
save('xreft','xreft')

%rCMP=x(4:6,1:N)-(ones(3,1)*((nv{1}'*x(4:6,1:N))./(nv{1}'*zeta(1:3,:)))).*zeta(1:3,:);

%plot(x')
figure(1)
subplot(311)
plot(T*(1:(size(x,2)-1)),LCoMest(1:3,:)')
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
plot(T*(1:(size(u,2)-1)),Totf)
legend('$f_x$','$f_y$','$f_z$')
subplot(132)
title('Total Torque around CoM'),hold on, grid on
plot(T*(1:(size(u,2)-1)),Tottorq)
legend('$\tau_x$','$\tau_y$','$\tau_z$')
subplot(133)
title('CoM Accelerations'),hold on, grid on
plot(T*(1:(size(u,2)-1)),ddx4p(1:3,:)')
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

axis([0 2.5 0 2.5 0 2.5])
view(150,30)
tHe=[0:T:T*N];
plot3(rHf(1,:),rHf(2,:),rHf(3,:),'k'); hold on, grid on
%%
%pause
for j=1:mcont
plot3(rP{j}(1,:),rP{j}(2,:),rP{j}(3,:),'g');
%quiver3(rP{j}(1,:),rP{j}(2,:),rP{j}(3,:),zeta(1+6*(j-1),:),zeta(2+6*(j-1),:),zeta(3+6*(j-1),:),100)
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

% figure(3)
% 
% for isimpl=1:N
%     deltais=min(find((isimpl<=Nsw)))-1;
%     plot3(x(4,isimpl),x(5,isimpl),x(6,isimpl),'x'); hold on, grid on
%     for j=ActCont{deltais}'
%     %plot3(R4pha(1,:),R4pha(2,:),R4pha(3,:),'r')
%     %plot3(rP{j}(1,:),rP{j}(2,:),rP{j}(3,:),'g'), hold on, grid on 
%     quivhand=quiver3(rP{j}(1,isimpl),rP{j}(2,isimpl),rP{j}(3,isimpl),zeta(1+6*(j-1),isimpl),zeta(2+6*(j-1),isimpl),zeta(3+6*(j-1),isimpl),0.0005);
%     end
%     pause(0.1)    
% end
end 