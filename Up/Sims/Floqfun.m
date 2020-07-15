function [maxFm,A]=Floqfun(x0,xParam,delx)
nsta=9;
M=zeros(nsta);
xm=eye(nsta);
Ninv=xm/delx;
x0p = Stride(x0, NaN,xParam);
if norm(x0p-x0)>0.001
    norm(x0p-x0)
    warning('No limit cycle');
end
for it=1:nsta
    x01 = x0+delx*xm(:,it);
    xf = Stride(x01, NaN,xParam);
    if isnan(xf)
        error('WARNING isnanxf');
    end
    M(:,it)=xf-x0p;
end
A=M*Ninv;
%%Part IV
efA=abs(eig(A));
maxFm=max(efA);
%[maxFm,im]=max(efA);
%tol1=1e-2;
% if maxFm>(1+tol1)
%     return
% else
%     efA(im)=[];
%     [~,im]=max(efA);efA(im)=[];
%     [~,im]=max(efA);efA(im)=[];
%     maxFm=max(efA);
% end
end