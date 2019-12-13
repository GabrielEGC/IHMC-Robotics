clc, clear all, close all
N=100;
M=50;%update
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
    fval;t0p=t0;vpoints=z(3*N+4:end);t0p=t0; nactflag=1;
%z=zeros(1,4*N+3);
tact=0;tactp=0;
tmax=1;
tv=dt*[1:tmax/dt];
for i=1:100000;
    t=t+dt;
    if 1/2*[x;dx]'*P*[x;dx]<c0
        Fn=-[kp kd]*[x;dx];
        contF=1;
    else
        contF=0;
        if (tact<t&&nactflag)
        vpoints=z(3*N+4:end);tactp=tact;
        nactflag=0;
        end

        if tact+M*T<t
        tic
        MPCPart1;
        c=[[x;dx;F];zeros(3*N,1)];
        [z,fval,exitflag] = quadprog(H,f,A,b,B,c,[],[],[],options);
        %vpoints=z(3*N+4:end);
        t0=t;timtQP=toc; tact=t+timtQP;nactflag=1;
        end
        
    end
    if ~contF
        n=floor((t-tactp)/T)+1;
        v=vpoints(n);
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
figure(3)
plot(tv,ContFv)
legend('Contv')