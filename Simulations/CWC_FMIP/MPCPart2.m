clc, clear all, close all
% PROBLEM SETUP, DO NOT CHANGE
A_sys = [1 0 .1 0;0 1 0 .1; 2 -2 1.2 0; 0 1.5 -1 .7];
B_sys = [0 0;0 0;1 -2;0 1];
N = 10;
x0 = [5;3;0;0];

% QP SETUP HERE
H = 2*[0.01*eye(4*N+4), zeros(4*N+4,2*N);zeros(4*N+4,2*N)', eye(2*N)];

f = zeros(6*N+4,1);


%%
fval_opt=Inf;
for m=1:N;
A1=kron(eye(m+1),[0,1,0,0]);
A2=[A1 zeros(m+1,6*N+4-4*(m+1))];
A3=kron(eye(N-m+1),[1,0,0,0]);
A4=[zeros(N-m+1,4*m) A3 zeros(N-m+1,2*N)];
Ap=[A2;A4];

b1=4*ones(m+1,1);
b2=2*ones(N-m+1,1);
bp=[b1;b2];
b3=-(2)*ones(m+1,1);
b4=-(-2)*ones(N-m+1,1);
bp2=[b3;b4];


A = [Ap;-Ap];
b = [bp;bp2];
%

B1=[kron(eye(N),A_sys) zeros(4*N,4) kron(eye(N),B_sys)];
BMPC1=[eye(4*N+4) zeros(4*N+4,2*N)]-[zeros(4,6*N+4);B1];
BMPC2=[zeros(4,4*N) eye(4) zeros(4,2*N)];
c1 = [x0;zeros(4*N,1)];
c2 = zeros(4,1);

%
B = [BMPC1;BMPC2];
c = [c1;c2];

% SOLVE QP
options = optimset('Display','Off');
      [z,fval,exitflag] = quadprog(H,f,A,b,B,c,[],[],[],options);
x=z(1:(4*N+4));
x=reshape(x,[4 N+1]);
u=z(4*N+5:end);
u=reshape(u,[2 N]);
exitflag
if (exitflag ==1)&&(fval<fval_opt)
   x_opt=x; 
   u_opt=u;
   m_opt=m;
   fval_opt=fval;
end

end
%%
% x(:,1)=x0;
% for i=1:N
%     x(:,i+1)=A_sys*x(:,i)+B_sys*u(:,i)
% end
x=x_opt;
plot(x')
%%
plot(x(1,:),x(2,:),'x'), hold on;
plot([-2 -2],[-1 5],'b')
plot([2 2],[-1 5],'b')
plot([-5.5 5.5],[2 2],'r')
plot([-5.5 5.5],[4 4],'r')