function [z,alr,tFlF]=Striderest(xSt0,x)
    run('InitSetParam.m')
    run('InitSetParamVar.m')
    [tStF,xStF] = StancePhaseOut(xSt0,I_bod,Lrv,Aznom,kpLe,kpAz,tau0,l0,gv,m,tpHM0,kdtau,kdLe,kdAz);
    phiStfin=castangphi(xStF,Lrv);
    [tFlF,xFlF] = FlightPhaseOutNoStride(xStF,Aznom,l0,phasphi,Lrv,gv,I_bod,tpHM0);
    alr=-phiStfin+phasphi(1)-phasphi(2)*tFlF;
    z=[phasphi(1);xFlF];%%%POSSIBLE REASON OF NO CONV: EVAL CHANGING TO CASTPHI
end