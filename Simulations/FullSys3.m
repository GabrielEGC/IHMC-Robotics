clc, clear all, close all
N=100;
M=10;%update
ksys=39.478;
bsys=6.283;
kp=300;kd=20;
vmax=200;
dt=10^-5;
P=[   0.004485757051635   0.000294369065356
   0.000294369065356   0.000019317722835];

Km=-[kp kd];
Aclsys=[0 1;-(ksys+ksys*kp) -(bsys+ksys*kd)];
nbar=(Km*Aclsys)';
contF=0;
c0=vmax^2/2*(nbar'*P^-1*nbar)^-1;
timQP=-Inf;
t=0;t0=0;
x=-1;dx=0;F=0;
    MPCPart1;
	c=[[x;dx;F];zeros(3*N,1)];
    [z,fval,exitflag] = quadprog(H,f,A,b,B,c,[],[],[],options);
    fval;t0p=t0;vpoints=z(3*N+4:end);t0p=t0;xpoints=reshape(z(1:3*N+3),[3 N+1]);
%z=zeros(1,4*N+3);
tmax=1;
tv=dt*[1:tmax/dt];
for i=1:100000;
    t=t+dt;
    if 1/2*[x;dx]'*P*[x;dx]<c0
        Fn=-[kp kd]*[x;dx];
        contF=1;
    else
        if floor((t-t0)/T)+1>M%(t0+timQP+0.3)<t
        t
        vpoints=z(3*N+4:end);t0p=t0;
        tic
        MPCPart1;
        c=[[x;dx;F];zeros(3*N,1)];
        [z,fval,exitflag] = quadprog(H,f,A,b,B,c,[],[],[],options);
        fval
        vpoints=z(3*N+4:end);xpoints=reshape(z(1:3*N+3),[3 N+1]);
        t0=t;contF=0;timtQP=toc;
        e4n=xpoints(1,:)-[x]*ones(1,N+1);[~,n]=min((e4n.^2));t0p=t0+T*(n-1);n
        end
    end
    if ~contF
        n=floor((t-t0p)/T)+1;
        v=vpoints(min(n,N));
        v=max(min(v,vmax),-vmax);
        F=F+v*dt;
    else
        v=(Fn-F)/dt;
        v=max(min(v,vmax),-vmax);
        F=F+v*dt;
    end

    ddx=ksys*F-bsys*dx-ksys*x;
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
plot(tv,Xv')
legend('x','dx','F')
subplot(212)
plot(tv,dFv')
legend('v')
figure(2)
plot(tv,Xv'*[-Km 1]')
legend('F+Km*x')
figure(3)
plot(tv,Nv)
legend('N')
figure(4)
plot(tv,ContFv)
legend('Contv')