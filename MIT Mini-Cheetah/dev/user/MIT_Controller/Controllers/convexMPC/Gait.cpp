#include "Gait.h"

// Offset - Duration Gait
ParkourGait::ParkourGait(int jumpConf, int _Np1, int _Np2, int _Np3, const std::string &name) //:_nIterations(nSegment)
{
  ////printf("HIUGHHH17");
  _name = name;
  setvNHacH(jumpConf, _Np1, _Np2, _Np3);
  _mpc_table = new int[_nIterations * 4];// allocate memory for MPC gait table
  _mpc_table_base = new int[_nIterations * 4];// allocate memory for MPC gait table
  getMpcTablebase();
  _CurrHybridMode=0;
  jumpConfGait=jumpConf;
  ////printf("HIUGHHH18");
}

OffsetDurationGait::OffsetDurationGait(int nSegment, Vec4<int> offsets, Vec4<int> durations, const std::string &name) :
  _offsets(offsets.array()),
  _durations(durations.array()),
  _nIterations(nSegment){
  _name = name;
  // allocate memory for MPC gait table
  _mpc_table = new int[nSegment * 4];
  _offsetsFloat = offsets.cast<float>() / (float) nSegment;
  _durationsFloat = durations.cast<float>() / (float) nSegment;
  _stance = durations[0];
  _swing = nSegment - durations[0];
  jumpConfGait=0;
}

ParkourGait::~ParkourGait() {
  delete[] _mpc_table;
  delete[] _mpc_table_base;
}

OffsetDurationGait::~OffsetDurationGait() {
  delete[] _mpc_table;
}

Vec4<float> ParkourGait::getis_slidcntct() {//
if(jumpConfGait==(int)3){
    if(_CurrHybridMode==0 || _CurrHybridMode==4){
      return Vec4<float>(0,0,0,0);
    }else if(_CurrHybridMode==1){
      return Vec4<float>(1,0,0,1);//02
    }else if(_CurrHybridMode==2){
      return Vec4<float>(0,0,0,0);
    }else if(_CurrHybridMode==3){
      return Vec4<float>(1,0,0,1);//02
    }else if(_CurrHybridMode==5){
      return Vec4<float>(0,1,1,0);//13
    }else if(_CurrHybridMode==6){
      return Vec4<float>(0,0,0,0);
    }else {
      return Vec4<float>(0,1,1,0);//13
    }
  }else{
    return Vec4<float>(0,0,0,0);//13
  }
}

Vec4<float> OffsetDurationGait::getis_slidcntct() {//
    return Vec4<float>(0,0,0,0);
}


void ParkourGait::recomputeVals(int jumpConf, int _Np1, int _Np2, int _Np3){
  ////printf("HIUGHHH17a");
  setvNHacH(jumpConf, _Np1, _Np2, _Np3);
  _mpc_table = new int[_nIterations * 4];// allocate memory for MPC gait table
  _mpc_table_base = new int[_nIterations * 4];// allocate memory for MPC gait table
  getMpcTablebase();//_mpc_table here?
  ////printf("HIUGHHH17b");
}

void OffsetDurationGait::recomputeVals(int jumpConf, int _Np1, int _Np2, int _Np3){
  (void) jumpConf;
  (void) _Np1;
  (void) _Np2;
  (void) _Np3;//
  _mpc_table = new int[_nIterations * 4];
  _offsetsFloat = _offsets.cast<float>() / (float) _nIterations;
  _durationsFloat = _durations.cast<float>() / (float) _nIterations;
  _stance = _durations[0];
  _swing = _nIterations - _durations[0];
}


void ParkourGait::setvNHacH(int jumpConf, int _Np1, int _Np2, int _Np3){
  MatDx1<int> NullMat;
  NullMat.conservativeResize(0,Eigen::NoChange);//MatDx1<int> FliMacNH;
  Vec4<int> HM1C(0,1,2,3);
  if(jumpConf==1){
    Eigen::Matrix<int,2,1> vNHdef;//NModes
    vNHdef<<_Np3-1,_Np1-1;
    NModes=vNHdef.rows();
    _nIterations=vNHdef(NModes-1)+1;
    vNH.resize(NModes, Eigen::NoChange);
    vNH=vNHdef;
    acH.conservativeResize(NModes, Eigen::NoChange);
    acH(0).resize(HM1C.rows(),Eigen::NoChange);acH(0)=HM1C;
    acH(1).resize(NullMat.rows(),Eigen::NoChange);acH(1)=NullMat;
    FliMacNH.resize(2, Eigen::NoChange);
    FliMacNH<<_Np3,_Np1;//+1
  }else if(jumpConf==2){//complete veryfiiiiiii
    Eigen::Matrix<int,9,1> vNHdef;//NModes
    NModes=vNHdef.rows();
    vNHdef<<1,3,5,7,9,11,13,15,17;
    _nIterations=vNHdef(NModes-1)+1;
    Vec2<int> HM3C(0,1);
    Vec2<int> HM2C(2,3);
    vNH.resize(NModes, Eigen::NoChange);
    vNH=vNHdef;
    acH.conservativeResize(NModes, Eigen::NoChange);
    acH(0).resize(HM1C.rows(),Eigen::NoChange);acH(0)=HM1C;
    acH(1).resize(HM2C.rows(),Eigen::NoChange);acH(1)=HM2C;
    acH(2).resize(NullMat.rows(),Eigen::NoChange);acH(2)=NullMat;
    acH(3).resize(HM3C.rows(),Eigen::NoChange);acH(3)=HM3C;
    acH(4).resize(HM1C.rows(),Eigen::NoChange);acH(4)=HM1C;
    acH(5).resize(HM2C.rows(),Eigen::NoChange);acH(5)=HM2C;
    acH(6).resize(NullMat.rows(),Eigen::NoChange);acH(6)=NullMat;
    acH(7).resize(HM3C.rows(),Eigen::NoChange);acH(7)=HM3C;
    acH(8).resize(HM1C.rows(),Eigen::NoChange);acH(8)=HM1C;
    FliMacNH.resize(4, Eigen::NoChange);
    FliMacNH<<4,6,12,14;
    ischained=true;//complete veryfiiiiiii
  }else if(jumpConf==3){
    Eigen::Matrix<int,8,1> vNHdef;//NModes
    NModes=vNHdef.rows();
    vNHdef<<5,8,9,12,18,21,22,25;//4,6,8,10,15,17,19,21;//vNHdef<<1,3,5,7,9,11,13,15;//7,10,12,15,23,26,28,31;//7,10,11,14,22,25,26,29;
    _nIterations=vNHdef(NModes-1)+1;
    Vec2<int> HM2C(1,2);//13
    Vec2<int> HM3C(0,3);//02
    vNH.resize(NModes, Eigen::NoChange);
    vNH=vNHdef;
    acH.conservativeResize(NModes, Eigen::NoChange);
    acH(0).resize(HM1C.rows(),Eigen::NoChange);acH(0)=HM1C;
    acH(1).resize(NullMat.rows(),Eigen::NoChange);acH(1)=NullMat;
    acH(2).resize(HM3C.rows(),Eigen::NoChange);acH(2)=HM3C;
    acH(3).resize(NullMat.rows(),Eigen::NoChange);acH(3)=NullMat;
    acH(4).resize(HM1C.rows(),Eigen::NoChange);acH(4)=HM1C;
    acH(5).resize(NullMat.rows(),Eigen::NoChange);acH(5)=NullMat;
    acH(6).resize(HM2C.rows(),Eigen::NoChange);acH(6)=HM2C;
    acH(7).resize(NullMat.rows(),Eigen::NoChange);acH(7)=NullMat;
    //FliMacNH.resize(4, Eigen::NoChange);
    //FliMacNH<<4,6,12,14;
  }else if(jumpConf==4){//Wrking double jump
    Eigen::Matrix<int,4,1> vNHdef;//NModes
    NModes=vNHdef.rows();
    vNHdef<<_Np3-1,_Np1-1,_Np1-1+_Np2+_Np3,_Np1*2-1+_Np2;//3,15,21,33
    _nIterations=vNHdef(NModes-1)+1;
    vNH.resize(NModes, Eigen::NoChange);
    vNH=vNHdef;
    acH.conservativeResize(NModes, Eigen::NoChange);
    acH(0).resize(HM1C.rows(),Eigen::NoChange);acH(0)=HM1C;
    acH(1).resize(NullMat.rows(),Eigen::NoChange);acH(1)=NullMat;
    acH(2).resize(HM1C.rows(),Eigen::NoChange);acH(2)=HM1C;
    acH(3).resize(NullMat.rows(),Eigen::NoChange);acH(3)=NullMat;
    FliMacNH.resize(4, Eigen::NoChange);
    FliMacNH<<_Np3,_Np1,_Np1+_Np2+_Np3,_Np1*2+_Np2;//(4,16,22,34)
  }else if(jumpConf==5){//BackFlipStab
    Eigen::Matrix<int,4,1> vNHdef;//NModes
    NModes=vNHdef.rows();
    vNHdef<<_Np3-1,_Np3-1+_Np2,_Np1-1,_Np1-1+4;//3,15,21,33
    _nIterations=vNHdef(NModes-1)+1;
    vNH.resize(NModes, Eigen::NoChange);
    vNH=vNHdef;
    Vec2<int> HM2C(2,3);//2,3
    acH.conservativeResize(NModes, Eigen::NoChange);
    acH(0).resize(HM1C.rows(),Eigen::NoChange);acH(0)=HM1C;
    acH(1).resize(HM2C.rows(),Eigen::NoChange);acH(1)=HM2C;
    acH(2).resize(NullMat.rows(),Eigen::NoChange);acH(2)=NullMat;
    acH(3).resize(HM1C.rows(),Eigen::NoChange);acH(3)=HM1C;
    FliMacNH.resize(3, Eigen::NoChange);
    FliMacNH<<_Np2+_Np3,_Np1,_Np1+4;//(4,16,22,34)
  }else if(jumpConf==6){//SideFlipStab
    Eigen::Matrix<int,4,1> vNHdef;//NModes
    NModes=vNHdef.rows();
    vNHdef<<_Np3-1,_Np3-1+_Np2,_Np1-1,_Np1-1+4;//3,15,21,33
    _nIterations=vNHdef(NModes-1)+1;
    vNH.resize(NModes, Eigen::NoChange);
    vNH=vNHdef;
    Vec2<int> HM2C(1,3);//2,3
    acH.conservativeResize(NModes, Eigen::NoChange);
    acH(0).resize(HM1C.rows(),Eigen::NoChange);acH(0)=HM1C;
    acH(1).resize(HM2C.rows(),Eigen::NoChange);acH(1)=HM2C;
    acH(2).resize(NullMat.rows(),Eigen::NoChange);acH(2)=NullMat;
    acH(3).resize(HM1C.rows(),Eigen::NoChange);acH(3)=HM1C;
    FliMacNH.resize(3, Eigen::NoChange);
    FliMacNH<<_Np2+_Np3,_Np1,_Np1+4;//(4,16,22,34)
  }else{//Wrking double jump
    Eigen::Matrix<int,5,1> vNHdef;//NModes
    NModes=vNHdef.rows();
    vNHdef<<_Np3-1,_Np1-1,_Np1-1+_Np2+_Np3,_Np1*2-1+_Np2,_Np1*2-1+_Np2+2;//3,15,21,33
    _nIterations=vNHdef(NModes-1)+1;
    vNH.resize(NModes, Eigen::NoChange);
    vNH=vNHdef;
    acH.conservativeResize(NModes, Eigen::NoChange);
    acH(0).resize(HM1C.rows(),Eigen::NoChange);acH(0)=HM1C;
    acH(1).resize(NullMat.rows(),Eigen::NoChange);acH(1)=NullMat;
    acH(2).resize(HM1C.rows(),Eigen::NoChange);acH(2)=HM1C;
    acH(3).resize(NullMat.rows(),Eigen::NoChange);acH(3)=NullMat;
    acH(4).resize(HM1C.rows(),Eigen::NoChange);acH(4)=HM1C;
    FliMacNH.resize(4, Eigen::NoChange);
    FliMacNH<<_Np3,_Np1,_Np1+_Np2+_Np3,_Np1*2+_Np2;//(4,16,22,34)
    ischained=true;
  }
  ////printf("HIUGHHH2");
  /*Array<MatDx1<int>,9,1> acHdef;
  acHdef<< Vec4<int>(0,1,2,3),
  Vec2<int>(2,3),
  NullMat,
  Vec2<int>(0,1),
  Vec4<int>(0,1,2,3),
  Vec2<int>(2,3),
  NullMat,
  Vec2<int>(0,1),
  Vec4<int>(0,1,2,3);
  vNH.resize(NModes, Eigen::NoChange);
  vNH=vNHdef;
  acH.conservativeResize(NModes, Eigen::NoChange);
  for(int HCi = 0; HCi<NModes; HCi++){
        acH(HCi).resize(acHdef(HCi).rows(),Eigen::NoChange);//conservativeResize
        acH(HCi)=acHdef(HCi);
  }*/
}

/*void ParkourGait::getvNHacH(){  
}*/

Vec4<float> ParkourGait::getContactState() {//
////printf("HIUGHHH3");
  Array4f progress(0,0,0,0);
  int NacH=acH(_CurrHybridMode).rows();//carefull in last dt before liftoff
  for(u16 legact=0;legact<NacH;legact++){
    u16 leg=acH(_CurrHybridMode)(legact);
    u16 lhm=1;u16 uhm=1;
    int _contdur;
    for(u16 lhmi=1;lhmi<(_CurrHybridMode+1);lhmi++){
      if(_mpc_table_base[4*vNH(_CurrHybridMode-lhm)+leg]<0.5){ break;} lhm++;}
    for(u16 uhmi=1;uhmi<(NModes-_CurrHybridMode);uhmi++){
      if(_mpc_table_base[4*vNH(_CurrHybridMode+uhm)+leg]<0.5){ break;} uhm++;}
    if(lhm!=(_CurrHybridMode+1)){
      _contdur=vNH(_CurrHybridMode+uhm-1)-vNH(_CurrHybridMode-lhm);
    }else{
      _contdur=vNH(_CurrHybridMode+uhm-1)+1;
    }
    ////printf("\n leg %i",leg);
    ////printf("uhm %i",uhm);
    ////printf("vNH %i",vNH(_CurrHybridMode+uhm-1));
    ////printf("CH %i",_CurrHybridMode);
    ////printf("_phase %f",_phase);
    ////printf("_contdur %i",_contdur);
    ////printf("_nIterations %i \n",_nIterations);
    progress[leg]=1-(vNH(_CurrHybridMode+uhm-1)+1-_phase*_nIterations)/_contdur;//8vNH(NModes-1)verifyyyyy
  }
  ////printf("HIUGHHH4");
  return progress.matrix();
}

Vec4<float> ParkourGait::getSwingState() {//assume acH is ordered
////printf("HIUGHHH5");
  Array4f progress(0,0,0,0);
  u16 legact=0;
  for(u16 leg=0;leg<4;leg++){
    if(legact<acH(_CurrHybridMode).rows()){
      if(leg==acH(_CurrHybridMode)(legact)){legact++; continue;}
    }

    u16 lhm=1;u16 uhm=1;
    float _contdur;
    for(u16 lhmi=1;lhmi<(_CurrHybridMode+1);lhmi++){
      if(_mpc_table_base[4*vNH(_CurrHybridMode-lhm)+leg]>0.5){ break;} lhm++;}
    for(u16 uhmi=1;uhmi<(NModes-_CurrHybridMode);uhmi++){
      if(_mpc_table_base[4*vNH(_CurrHybridMode+uhm)+leg]>0.5){ break;} uhm++;}
    if(lhm!=(_CurrHybridMode+1)){
      _contdur=vNH(_CurrHybridMode+uhm-1)-vNH(_CurrHybridMode-lhm);
    }else{
      _contdur=vNH(_CurrHybridMode+uhm-1)+1;
    }
    progress[leg]=1-(vNH(_CurrHybridMode+uhm-1)+1-_phase*_nIterations)/_contdur;//8vNH(NModes-1)
  }
  ////printf("HIUGHHH6");
  return progress.matrix();
}

Vec4<float> OffsetDurationGait::getContactState() {
  Array4f progress = _phase - _offsetsFloat;
  for(int i = 0; i < 4; i++){
    if(progress[i] < 0) progress[i] += 1.;
    if(progress[i] > _durationsFloat[i]){
      progress[i] = 0.;}
    else{progress[i] = progress[i] / _durationsFloat[i];}
  }
  //////printf("contact state: %.3f %.3f %.3f %.3f\n", progress[0], progress[1], progress[2], progress[3]);
  return progress.matrix();
}

Vec4<float> OffsetDurationGait::getSwingState(){
  Array4f swing_offset = _offsetsFloat + _durationsFloat;
  for(int i = 0; i < 4; i++)
    if(swing_offset[i] > 1) swing_offset[i] -= 1.;
  Array4f swing_duration = 1. - _durationsFloat;
  Array4f progress = _phase - swing_offset;
  for(int i = 0; i < 4; i++){
    if(progress[i] < 0) progress[i] += 1.f;
    if(progress[i] > swing_duration[i]||swing_duration[i]==0){
      progress[i] = 0.;}
    else{
      progress[i] = progress[i] / swing_duration[i];}
  }//////printf("swing state: %.3f %.3f %.3f %.3f\n", progress[0], progress[1], progress[2], progress[3]);
  return progress.matrix();
}

void ParkourGait::getMpcTablebase(){//////printf("MPC table:\n");
////printf("HIUGHHH7");
for (u16 N3 = 0; N3 < (4*_nIterations);N3++)
  _mpc_table_base[N3] = 0;
  
for (u16 iN = 0; iN < vNH(0)+1;iN++){  
  u16 nacHN=acH(0).rows();
  for (u16 N3 = 0; N3 < nacHN;N3++)
    _mpc_table_base[iN*4 + acH(0)(N3)] = 1;
}
for (u16 N = 1; N < NModes;N++){
  for (u16 iN = vNH(N-1)+1; iN < vNH(N)+1;iN++){  
    u16 nacHN=acH(N).rows();
    for (u16 N3 = 0; N3 < nacHN;N3++)
      _mpc_table_base[iN*4 + acH(N)(N3)] = 1; ////printf("HIUGHHH8");
}}}

int* ParkourGait::getMpcTable(){//////printf("MPC table:\n");
////printf("HIUGHHH9");
  for(int i = 0; i < _nIterations; i++){
    int iter = (i + _iteration) % _nIterations;//+1
    for(int j = 0; j < 4; j++){
      _mpc_table[i*4 + j] = _mpc_table_base[4*iter + j];//////printf("%d ", _mpc_table[i*4 + j]);
    }//////printf("\n");
  }  ////printf("HIUGHHH10");
  return _mpc_table;
}

int* OffsetDurationGait::getMpcTable(){
  //////printf("MPC table:\n");
  for(int i = 0; i < _nIterations; i++){
    int iter = (i + _iteration) % _nIterations;//+1
    Array4i progress = iter - _offsets;
    for(int j = 0; j < 4; j++){
      if(progress[j] < 0) progress[j] += _nIterations;
      if(progress[j] < _durations[j])
        _mpc_table[i*4 + j] = 1;
      else
        _mpc_table[i*4 + j] = 0;//////printf("%d ", _mpc_table[i*4 + j]);
    }//////printf("\n");
  }
  return _mpc_table;
}

void ParkourGait::setIterations(int iterationsPerMPC, int currentIteration){
  ////printf("HIUGHHH11");
  _phase = (float)((currentIteration+gaitphascorr)% (iterationsPerMPC*_nIterations))/(float)(iterationsPerMPC*_nIterations);
  _iteration = (int) (_phase*_nIterations);  //_iteration = ((currentIteration+gaitphascorr) / iterationsPerMPC) % _nIterations;
  ////printf("_iterationpast: %i \n",_iteration);  //_iteration = (int) (_phase*_nIterations);  ////printf("_iterationnew: %i \n",_iteration);  ////printf("_phase: %f \n",_phase);
  _CurrHybridMode=-1000;
  if (vNH(0)>(_iteration-1)){
    _CurrHybridMode=0;
  }else{
    for(u16 HM=1;HM<NModes;HM++){
      if ((vNH(HM)>(_iteration-1))&&(_iteration>vNH(HM-1))){
        _CurrHybridMode=HM;
          return;
}}}}

void ParkourGait::setItCorrected(int iterationsPerMPC, int currentIteration){
  ////printf("HIUGHHH11");
  gaitphascorr = (vNH(_CurrHybridMode)+1)*iterationsPerMPC-(currentIteration% (iterationsPerMPC*_nIterations));
 //// ///printf("gaitphascorr: %i \n",gaitphascorr);
 //// ///printf("_CurrHybridMode: %i \n",_CurrHybridMode);
  _phase = (float)((currentIteration+gaitphascorr)% (iterationsPerMPC*_nIterations))/(float)(iterationsPerMPC*_nIterations);
  _iteration = (int) (_phase*_nIterations);  //_iteration = ((currentIteration+gaitphascorr) / iterationsPerMPC) % _nIterations;
 //// ///printf("_iteration: %i \n",_iteration);  
printf("_iterationnewitcorr: %i \n",_iteration);  
printf("_phasenewitcorr: %f \n",_phase);
  _CurrHybridMode=-1000;
  if (vNH(0)>(_iteration-1)){
    _CurrHybridMode=0;
  }else{
    for(u16 HM=1;HM<NModes;HM++){
      if ((vNH(HM)>(_iteration-1))&&(_iteration>vNH(HM-1))){
        _CurrHybridMode=HM;
       //// ///printf("CurrentHM: %i",_CurrHybridMode);
          return;
}}}}

void OffsetDurationGait::setItCorrected(int iterationsPerMPC, int currentIteration){
  gaitphascorr = -(currentIteration% (iterationsPerMPC*_nIterations));
  _phase = (float)0;
  _iteration = (int) (0);  //_iteration = ((currentIteration+gaitphascorr) / iterationsPerMPC) % _nIterations;
  printf("gaitphascorr: %i \n",gaitphascorr);
}

void ParkourGait::setItCorrectedStart(int iterationsPerMPC, int currentIteration){
  ////printf("HIUGHHH11");
  gaitphascorr = -(currentIteration% (iterationsPerMPC*_nIterations));
  _phase = (float)0;
  _iteration = (int) (0);  //_iteration = ((currentIteration+gaitphascorr) / iterationsPerMPC) % _nIterations;
  _CurrHybridMode=0;
  printf("gaitphascorr: %i \n",gaitphascorr);
}

void OffsetDurationGait::setItCorrectedStart(int iterationsPerMPC, int currentIteration){
  gaitphascorr = -(currentIteration% (iterationsPerMPC*_nIterations));
  _phase = (float)0;
  _iteration = (int) (0);  //_iteration = ((currentIteration+gaitphascorr) / iterationsPerMPC) % _nIterations;
  printf("gaitphascorr: %i \n",gaitphascorr);
}

//vNHdef<<1,3,5,7,9,11,13,15,17;
void ParkourGait::setIterationssb(int _iteration1, float _phase1){
  ////printf("HIUGHHH21");
  _iteration = _iteration1;
  _phase = _phase1;
 //// ///printf("\n_iteration: %i",_iteration);  
  ////printf("_iterationnew: %i \n",_iteration);  
 //// ///printf("_phase: %f \n",_phase);
  _CurrHybridMode=-1000;
  if (vNH(0)>(_iteration-1)){
    _CurrHybridMode=0;
  }else{
    for(u16 HM=1;HM<NModes;HM++){
      if ((vNH(HM)>(_iteration-1))&&(_iteration>vNH(HM-1))){
        _CurrHybridMode=HM;
       //// ///printf("CurrentHM: %i \n",_CurrHybridMode);
        return;
}}}}////printf("HIUGHHH22");

void OffsetDurationGait::setIterations(int iterationsPerMPC, int currentIteration){
  //_iteration = ((currentIteration+gaitphascorr)/iterationsPerMPC) % _nIterations;
  _phase = (float)((currentIteration+gaitphascorr)% (iterationsPerMPC*_nIterations))/(float)(iterationsPerMPC*_nIterations);
  _iteration = (int)(_phase*_nIterations);
}

void OffsetDurationGait::setIterationssb(int _iteration1, float _phase1){
  _iteration = _iteration1;
  _phase = _phase1;
}

void ParkourGait::set_nIterations(int _nIt){
  (void)_nIt;//_nIterations=_nIt;
}

void OffsetDurationGait::set_nIterations(int _nIt){
  _nIterations=_nIt;
}

void ParkourGait::setdurations(int dur){
  //Array4i((int)data.userParameters->Npar3,(int)data.userParameters->Npar3,(int)data.userParameters->Npar3,(int)data.userParameters->Npar3);  
  //_durations=Array4i(dur,dur,dur,dur);
  (void) dur;
}

void OffsetDurationGait::setdurations(int dur){
  //Array4i((int)data.userParameters->Npar3,(int)data.userParameters->Npar3,(int)data.userParameters->Npar3,(int)data.userParameters->Npar3);  
  _durations=Array4i(dur,dur,dur,dur);
}

int ParkourGait::getCurrentGaitPhase() {
  return _iteration;
}

int OffsetDurationGait::getCurrentGaitPhase() {
  return _iteration;
}

int ParkourGait::getCurrentHybridMode() {
  return _CurrHybridMode;
}

int OffsetDurationGait::getCurrentHybridMode() {
  return 1;//VERIFYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYY
}

int ParkourGait::getItNxtHybMode() {
  return (vNH(_CurrHybridMode)+1);//it
}

int OffsetDurationGait::getItNxtHybMode() {
  return 0;//VERIFYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYY
}

float OffsetDurationGait::getCurrentSwingTime(float dtMPC, int leg) {//GET TOTAL SWING TIME
  (void)leg;
  return dtMPC * _swing;
}

float ParkourGait::getCurrentSwingTime(float dtMPC, int leg) {//ifsstance, next swing time, if swing, this swing time
////printf("HIUGHHH13");
  float _contdur=0;
  u16 lhm=1;u16 uhm=1;
  if(_mpc_table_base[4*vNH(_CurrHybridMode)+leg]<0.5){ 
    for(u16 lhmi=1;lhmi<(_CurrHybridMode+1);lhmi++){
      if(_mpc_table_base[4*vNH(_CurrHybridMode-lhm)+leg]>0.5){ break;} lhm++;}
    for(u16 uhmi=1;uhmi<(NModes-_CurrHybridMode);uhmi++){
      if(_mpc_table_base[4*vNH(_CurrHybridMode+uhm)+leg]>0.5){ break;} uhm++;}
    if(lhm!=(_CurrHybridMode+1)){
    _contdur=vNH(_CurrHybridMode+uhm-1)-vNH(_CurrHybridMode-lhm);
    }else{
      _contdur=vNH(_CurrHybridMode+uhm-1)+1;
    }
  }else{
    for(u16 lhmi=1;lhmi<(NModes+1);lhmi++){
      if(_mpc_table_base[4*vNH((_CurrHybridMode+lhm)%NModes)+leg]<0.5){ break;} lhm++;}//note change of sign
    uhm=lhm;
    for(u16 uhmi=lhm;uhmi<(NModes+1+lhm);uhmi++){
      if(_mpc_table_base[4*vNH((_CurrHybridMode+uhmi)%NModes)+leg]>0.5){ break;} uhm++;}//VERIFY LATERRRRRR

    if(((_CurrHybridMode+lhm)%NModes)!=0){
      _contdur=vNH((_CurrHybridMode+uhm-1)%NModes)-vNH((_CurrHybridMode+lhm)%NModes-1);
      if (_contdur<0)
        _contdur=vNH((_CurrHybridMode+uhm-1)%NModes)-vNH((_CurrHybridMode+lhm)%NModes-1)+_nIterations;
    }else{
      _contdur=vNH((_CurrHybridMode+uhm-1)%NModes)+1;
  }}
  ////printf("HIUGHHH14");
  return (dtMPC*_contdur);
}

float ParkourGait::getCurrentStanceTime(float dtMPC, int leg) {//ifswing, next stance time, if stance, this stance time
////printf("HIUGHHH15");
  float _contdur=0;
  u16 lhm=1;u16 uhm=1;
  if(_mpc_table_base[4*vNH(_CurrHybridMode)+leg]>0.5){ 
    for(u16 lhmi=1;lhmi<(_CurrHybridMode+1);lhmi++){
      if(_mpc_table_base[4*vNH(_CurrHybridMode-lhm)+leg]<0.5){ break;} lhm++;}
    for(u16 uhmi=1;uhmi<(NModes-_CurrHybridMode);uhmi++){
      if(_mpc_table_base[4*vNH(_CurrHybridMode+uhm)+leg]<0.5){ break;} uhm++;}
    if(lhm!=(_CurrHybridMode+1)){
      _contdur=vNH(_CurrHybridMode+uhm-1)-vNH(_CurrHybridMode-lhm);
    }else{
      _contdur=vNH(_CurrHybridMode+uhm-1)+1;
    }
  }else{
    for(u16 lhmi=1;lhmi<(NModes+1);lhmi++){
      if(_mpc_table_base[4*vNH((_CurrHybridMode+lhm)%NModes)+leg]>0.5){ break;} lhm++;}//note change of sign
    uhm=lhm;
    for(u16 uhmi=lhm;uhmi<(NModes+1+lhm);uhmi++){
      if(_mpc_table_base[4*vNH((_CurrHybridMode+uhmi)%NModes)+leg]<0.5){ break;} uhm++;}//VERIFY LATERRRRRR
    
    if(((_CurrHybridMode+lhm)%NModes)!=0){
      _contdur=vNH((_CurrHybridMode+uhm-1)%NModes)-vNH(((_CurrHybridMode+lhm)%NModes)-1);
      if (_contdur<0)
        _contdur=vNH((_CurrHybridMode+uhm-1)%NModes)-vNH(((_CurrHybridMode+lhm)%NModes)-1)+_nIterations;
    }else{
      _contdur=vNH((_CurrHybridMode+uhm-1)%NModes)+1;
  }}
  ////printf("HIUGHHH16");
  return (dtMPC*_contdur);
}

float OffsetDurationGait::getCurrentStanceTime(float dtMPC, int leg) {
  (void) leg;//warning error
  return dtMPC * _stance;
}

int ParkourGait::getHoriz(){
  ////printf("HIUGHHH18a, %i",_nIterations);
  ////printf("HIUGHHH18b");
  return _nIterations;
}
int OffsetDurationGait::getHoriz(){return _nIterations;}


bool ParkourGait::getFlightState(){
  int iter = _iteration % _nIterations;//+1
  for(int j = 0; j < 4; j++){
    if(_mpc_table_base[4*iter + j]>0) return false;
  }
  return true;
}
bool OffsetDurationGait::getFlightState(){
  //
  int iter = (_iteration) % _nIterations;
  Array4i progress = iter - _offsets;
  for(int j = 0; j < 4; j++){
  if(progress[j] < 0) progress[j] += _nIterations;
  if(progress[j] < _durations[j]) return false;
  }
  return true;
}
//vNHdef<<3,15,21,33;FliMacNH<<_Np3,_Np1,_Np1+_Np2+_Np3,_Np1*2+_Np2;//(4,16,22,34)
float ParkourGait::get_phaseLO(){
  if(countflystage==1){
    return FliMacNH(0)/(float)_nIterations;
  }else if(countflystage==2){
    return FliMacNH(2)/(float)_nIterations;
  }else{
    return 0;
  }
}

float ParkourGait::get_phasela(){
  if(countflystage==1){
    return FliMacNH(1)/(float)_nIterations;
  }else if(countflystage==2){
    return FliMacNH(3)/(float)_nIterations;//prev "1"
  }else{
    return 0;
  }
}

float OffsetDurationGait::get_phaseLO(){
  if(countflystage==1){
    return _durationsFloat[0];
  }else{
    return 0;
  }
}

float OffsetDurationGait::get_phasela(){
  if(countflystage==1){
    return 1;
  }else{
    return 0;
  }
}

void ParkourGait::debugPrint() {}
void OffsetDurationGait::debugPrint() {}

int ParkourGait::getFliMacNH(int NumFlightM) {
  if(ischained){
    if (NumFlightM==1){
      return FliMacNH(1);
    }else{
      return FliMacNH(3);
    }
  }else{
    return FliMacNH(1);
  }
}
int OffsetDurationGait::getFliMacNH(int NumFlightM) {
  (void) NumFlightM;
  return _nIterations;
}