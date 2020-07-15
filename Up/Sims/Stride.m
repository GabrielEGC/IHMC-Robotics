function z=Stride(xIn,wrV,x,alr)
    InitSetUpStride;
%    xSt0=[xSt0(1:3);xIn];
%    [Fz0,~,~] = LiftOffEventsFcn(0,xSt0,@Fextfun,Lr,Aznom,kpLe,kpAz,tau0,l0,I_bod,kdtau);
%     if (norm(xSt0)>1e+12||Fz0(1)<1e-12)
%         %Fz0(1)
%         z=wrV;disp('ret. wrV1');%return;
%     end
    [tStF,xStF] = StancePhaseOut(xSt0,I_bod,Lrv,Aznom,kpLe,kpAz,tau0,l0,gv,m,tpHM0,kdtau,kdLe,kdAz);
    phasphiStfin=castangphi(xStF,Lrv);
   [Fz0,~,~]=LiftOffEventsFcn(tStF,xStF,@Fextfun,Lrv,Aznom,kpLe,kpAz,tau0,l0,I_bod,kdtau,kdLe,kdAz);
    if ((tStF==tpHM0)||(tStF<0.001)||xStF(3)<0||abs(Fz0(1))>1e-2)
        xStF(3)
        abs(Fz0(1))
        tStF
        z=wrV;disp('ret. wrV2');return;
    end
    [tFlF,xFlF,phasphiTD] = FlightPhaseOutStride(xStF,Aznom,l0,phasphi,Lrv,gv,I_bod,tpHM0,alr,phasphiStfin);
    if ((tFlF==tpHM0))%||(abs(rFlFfoot(3))>1e-4)
        tFlF
        z=wrV;disp('ret. wrV4');return;
    end
    z=[phasphiTD;xFlF];
end