clc, clear all, close all
% PROBLEM SETUP, DO NOT CHANGE
A_sys = [1 0 .1 0;0 1 0 .1; 2 -2 1.2 0; 0 1.5 -1 .7];
B_sys = [0 0;0 0;1 -2;0 1];
N = 10;
x0 = 10*randn(4,1);
% QP SETUP HERE
H = 2*[0.01*eye(4*N+4), zeros(4*N+4,2*N);zeros(4*N+4,2*N)', eye(2*N)];

f = zeros(6*N+4,1);
A = zeros(0,6*N+4);
b = zeros(0);

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
u=z(4*N+5:end);
u=reshape(u,[2 N]);
%%
x(:,1)=x0;
for i=1:N
    x(:,i+1)=A_sys*x(:,i)+B_sys*u(:,i);
end
plot(x')