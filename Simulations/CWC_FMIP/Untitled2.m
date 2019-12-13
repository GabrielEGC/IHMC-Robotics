clear all
load('rVtb1')
eFlag1=eFlag;
rVt1=rVt;load('rVtb2')

%
figure(3);
rVt1ef1=find(eFlag1>-3);
%rVt1ef2=find(eFlag>0);
%rVt1ef1=find(eFlag.*eFlag1==1)
%
plot3([rVt1(1,rVt1ef1);rVt(1,rVt1ef1)],[rVt1(2,rVt1ef1);rVt(2,rVt1ef1)],[rVt1(3,rVt1ef1);rVt(3,rVt1ef1)]), hold on, grid on
axis([-1 1 -1 1 0 1.6])
%%

figure(3);
rVt1ef=find(eFlag>0);
plot3(rVt(1,rVt1ef),rVt(2,rVt1ef),rVt(3,rVt1ef),'ob'), hold on, grid on