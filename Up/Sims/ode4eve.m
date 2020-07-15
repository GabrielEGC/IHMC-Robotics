function [tspan,Y] = ode4eve(odefun,tspan,y0,optode4eve,varargin)
%[tspan,Y,Ydot]
h = diff(tspan);

try
  f0 = feval(odefun,tspan(1),y0,varargin{:});
catch
  msg = ['Unable to evaluate the ODEFUN at t0,y0. ',lasterr];
  error(msg);  
end  

y0 = y0(:);   % Make a column vector.

neq = length(y0);
N = length(tspan);
Y = zeros(neq,N);
%Ydot = zeros(N,neq);
F = zeros(neq,4);

Y(:,1) = y0;

evefun=optode4eve.Events;
[posim1,~,dirim1]=feval(evefun,tspan(1),Y(:,1));%%CSTES
evet=false;

for i = 2:N
  ti = tspan(i-1);
  hi = h(i-1);
  yi = Y(:,i-1);
  F(:,1) = feval(odefun,ti,yi,varargin{:});  %Ydot(i,:) = F(:,1).';
  F(:,2) = feval(odefun,ti+0.5*hi,yi+0.5*hi*F(:,1),varargin{:});
  F(:,3) = feval(odefun,ti+0.5*hi,yi+0.5*hi*F(:,2),varargin{:});  
  F(:,4) = feval(odefun,tspan(i),yi+hi*F(:,3),varargin{:});
  Y(:,i) = yi + (hi/6)*(F(:,1) + 2*F(:,2) + 2*F(:,3) + F(:,4));
  [posi,~,~]=feval(evefun,tspan(i),Y(:,i));
  if any(posi.*posim1<=0 & (dirim1.*(posi-posim1))>0)
      evet=true;break;
  end
  posim1=posi;
end
if evet
    ein=find(posi.*posim1<=0 & (dirim1.*(posi-posim1))>0);
    delti_f=posim1(ein)/(posim1(ein)-posi(ein));
    ti_f=ti+delti_f*hi;
    Y(:,i)=Y(:,i-1)+delti_f*(Y(:,i)-Y(:,i-1));%1st order
    tspan=[tspan(1:i-1),ti_f];
    Y = Y(:,1:i).';
else
    Y = Y.';
end
%Ydot=Ydot(1:i,:);
end