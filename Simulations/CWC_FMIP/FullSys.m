clc, clear all, close all
N=100;
M=10;%update
ksys=39.478;
bsys=6.283;
kp=300;kd=20;
vmax=200;
dt=10^-5;
P=[0.243839640250665   0.008321491867835
   0.008321491867835   0.000345913223515];

Km=-[kp kd];
Aclsys=[0 1;-(ksys+kp) -(bsys+kd)];
nbar=(Km*Aclsys)';
contF=0;
c0=vmax^2/2*(nbar'*P^-1*nbar)^-1;
timQP=-Inf;
t=0;t0=0;
x=-1;dx=0.5;F=1;
    MPCPart1;
	c=[[x;dx;F];zeros(3*N,1)];
    [z,fval,exitflag] = quadprog(H,f,A,b,B,c,[],[],[],options);
    fval;t0p=t0;vpoints=z(3*N+4:end);t0p=t0;xpoints=reshape(z(1:3*N+3),[3 N+1]);
%z=zeros(1,4*N+3);
emax=0.6;
tmax=1;
tv=dt*[1:tmax/dt];
Pcyl=10^2*[   1.189882950697608   0.079320369897852   0.003947265989137
   0.079320369897852   0.005353334057673   0.000505007706091
   0.003947265989137   0.000505007706091   0.000904182047887];
CylTrigg=0;
BBtrigg=0.8;
BBtrigg2=min(1.1*BBtrigg,1);
for i=1:length(tv);
    t=t+dt;
    if ([x;dx;F]'*Pcyl*[x;dx;F]<BBtrigg&&x<emax)||(CylTrigg==1)%;1/2*[x;dx]'*P*[x;dx]<c0&&abs([-Km 1]*[x;dx;F])<vmax*dt%
        Fn=-[kp kd]*[x;dx];
        CylTrigg=1;contF=1;
    end
    if ([x;dx;F]'*Pcyl*[x;dx;F]>BBtrigg2||x>emax)
        CylTrigg=0;
        contF=0;
        if floor((t-t0)/T)+1>M%(t0+timQP+0.3)<t
        t
        vpoints=z(3*N+4:end);t0p=t0;
        tic
        MPCPart1;
        c=[[x;dx;F];zeros(3*N,1)];
        [z,fval,exitflag] = quadprog(H,f,A,b,B,c,[],[],[],options);
        fval
        vpoints=z(3*N+4:end);xpoints=reshape(z(1:3*N+3),[3 N+1]);
        t0=t;timtQP=toc;t0p=t0;
        end
    end
    if contF==0
        n=floor((t-t0p)/T)+1;
        v=vpoints(n);
        v=max(min(v,vmax),-vmax);
        F=F+v*dt;
    else
        v=(Fn-F)/dt;
        v=max(min(v,vmax),-vmax);
        F=F+v*dt;
    end

    ddx=F-bsys*dx-ksys*x;
    dx=dx+ddx*dt;
    x=x+dx*dt;
    Xv(:,i)=[x;dx;F];
    ContFv(:,i)=contF;
    Nv(:,i)=n;
    dFv(:,i)=v;
end

%%
figure(1)
subplot(211)
plot(tv,Xv(3,:)')
legend('F')
subplot(212)
plot(tv,dFv')
legend('v')
figure(2)
subplot(211)
plot(tv,Xv'*[-Km 1]')
subplot(212)
plot(tv,Xv(1:2,:)')
legend('x','dx')
figure(3)
plot(tv,Nv)
legend('N')
figure(3)
subplot(211)
plot(tv,ContFv)
legend('Contv')
subplot(212)
plot(tv,sum((Xv'*Pcyl)'.*Xv))
legend('PROA')