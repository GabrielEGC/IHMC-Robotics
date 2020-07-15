function [maxFm,A]=Floqfunopt(x0,xParam,delx,x0p,alr)
nsta=9;
M=zeros(nsta);
xm=eye(nsta);
Ninv=xm/delx;
% if norm(x0p-x0)>0.001
%     norm(x0p-x0)
%     warning('No limit cycle');
% end
for it=1:nsta
    x01 = x0+delx*xm(:,it);
    xf = Stride(x01, NaN,xParam,alr);
    if isnan(xf)
        error('WARNING isnanxf');
    end
    M(:,it)=xf-x0p;
end
A=M*Ninv;
%assignin('base','Ab',A);
%%Part IV
maxFm=abs(eig(A));
% efA=sort(abs(eig(A)))
% maxFm=efA(nsta);
%maxFm=sum(efA(end-1:end))
%+10*efA(nsta)-efA(1);%efA(nsta-1)+efA(nsta)
%maxFm=max(abs(eig(A)));
% if maxFm>(1+1e-3)
%     return
% else
%end

% [maxFm,im]=max(efA);
% tol1=1e-2;
% if maxFm>(1+tol1)
%     return
% else
%     efA(im)=[];
%     [~,im]=max(efA);efA(im)=[];
%     [~,im]=max(efA);efA(im)=[];
%     maxFm=max(efA);
% end
end