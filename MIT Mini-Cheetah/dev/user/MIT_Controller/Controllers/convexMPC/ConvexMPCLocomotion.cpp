#include <iostream>
#include <Utilities/Timer.h>
#include <Utilities/Utilities_print.h>
#include <Math/orientation_tools.h>
#include "ConvexMPCLocomotion.h"
#include "convexMPC_interface.h"
#include "../../../../common/FootstepPlanner/GraphSearch.h"
#include "common_types.h"
#include "../Printl/printheader.h"
#include <Utilities/pseudoInverse.h>
#define DRAW_DEBUG_SWINGS
#define DRAW_DEBUG_PATH
ConvexMPCLocomotion::ConvexMPCLocomotion(float _dt, int _iterations_between_mpc, MIT_UserParameters* parameters) :
  iterationsBetweenMPC(_iterations_between_mpc),
  horizonLength(10),
  dt(_dt),
  jumpingpg(0, (int)parameters->Npar1, (int)parameters->Npar2, (int)parameters->Npar3,"JumpParkour"),
  walkingpg(3, (int)parameters->Npar1, (int)parameters->Npar2, (int)parameters->Npar3,"WalkP"),
  jumpingbfspg(5, (int)parameters->Npar1, (int)parameters->Npar2, (int)parameters->Npar3,"JumpBFStabParkour"),
  jumpingsfspg(6, (int)parameters->Npar1, (int)parameters->Npar2, (int)parameters->Npar3,"JumpSideStabParkour"),
  trotting(horizonLength, Vec4<int>(0,5,5,0), Vec4<int>(5,5,5,5),"Trotting"),
  bounding(horizonLength, Vec4<int>(0,0,4,4),Vec4<int>(4,4,4,4),"Bounding"),//bounding(horizonLength, Vec4<int>(5,5,0,0),Vec4<int>(3,3,3,3),"Bounding"),
  pronking(horizonLength, Vec4<int>(0,0,0,0),Vec4<int>(4,4,4,4),"Pronking"),
  jumpingodg((int)parameters->Npar1, Vec4<int>(0,0,0,0), Vec4<int>((int)parameters->Npar3,(int)parameters->Npar3,(int)parameters->Npar3,(int)parameters->Npar3), "Jumping"),
  galloping(horizonLength, Vec4<int>(0,2,7,9),Vec4<int>(4,4,4,4),"Galloping"),  //galloping(horizonLength, Vec4<int>(0,2,7,9),Vec4<int>(6,6,6,6),"Galloping"),//galloping(horizonLength, Vec4<int>(0,2,7,9),Vec4<int>(3,3,3,3),"Galloping"),
  standing(horizonLength, Vec4<int>(0,0,0,0),Vec4<int>(10,10,10,10),"Standing"),
  trotRunning(horizonLength, Vec4<int>(0,5,5,0),Vec4<int>(4,4,4,4),"Trot Running"),//trotRunning(horizonLength, Vec4<int>(0,5,5,0),Vec4<int>(3,3,3,3),"Trot Running"),
  walking(horizonLength, Vec4<int>(0,3,5,8), Vec4<int>(5,5,5,5), "Walking"),
  walking2(horizonLength, Vec4<int>(0,5,5,0), Vec4<int>(7,7,7,7), "Walking2"),
  pacing(horizonLength, Vec4<int>(5,0,5,0),Vec4<int>(5,5,5,5),"Pacing")
{
  //I_bod.setZero();
  //I_bod.diagonal() = Vec3<float>(.07f, 0.26f, 0.242f);
  _parameters = parameters;
  dtMPC = dt * iterationsBetweenMPC;
  default_iterations_between_mpc = iterationsBetweenMPC;
  ////printf("[Convex MPC] dt: %.3f iterations: %d, dtMPC: %.3f\n", dt, iterationsBetweenMPC, dtMPC);
  Rpla[0].setIdentity();
  Rpla[1].setIdentity();
  Rpla[2].setIdentity();
  Rpla[3].setIdentity();
  setup_problem(dtMPC, horizonLength, muv, 120, -1, Rpla);
  for(int i = 0; i < 4; i++)
    firstSwing[i] = true;
  pBody_des.setZero();
  vBody_des.setZero();
  aBody_des.setZero();
  Wpla.resize(4,3);
  Wpla.setZero();
  Wpla(0,0)=1;Wpla(1,0)=1;Wpla(2,0)=1;Wpla(3,0)=1;
  zfeet.setZero();
  apla.setZero();
  aplatarget=&apla;
  Rplatarget[0]=&Rpla[0];
  Rplatarget[1]=&Rpla[1];
  Rplatarget[2]=&Rpla[2];
  Rplatarget[3]=&Rpla[3];
  //apla2<<0,0,-0.269353;
  apla2b<<-0.0764,0.349066,-0.00427825;
  apla2<<-0.0764,0.349066,-0.00427825;
  apla2.setZero();//Uncomment for horizontal planes when noncopjum is 1 DEFINE
  apla2b.setZero();//Uncomment for horizontal planes when noncopjum is 1 DEFINE
  /*apla2b<<0.1016,0,0;
  apla2<<0.1016,0,0;*/
  rpydesWBC.setZero();
  drpydesWBC.setZero();
  Rpla2[0]=aplatoRpla(apla2);Rpla2[1]=Rpla2[0];Rpla2[2]=Rpla2[0];Rpla2[3]=Rpla2[0];
  v_des_robot.setZero();
  v_des_worldWBC.setZero();
  Np1=parameters->Npar1;
  Np2=parameters->Npar2;
  Np3=parameters->Npar3;
}

void ConvexMPCLocomotion::initialize(){
  for(int i = 0; i < 4; i++) firstSwing[i] = true;
  firstRun = true;
}

void ConvexMPCLocomotion::recompute_timing(int iterations_per_mpc) {
  iterationsBetweenMPC = iterations_per_mpc;
  dtMPC = dt * iterations_per_mpc;
}

void ConvexMPCLocomotion::_SetupCommand(ControlFSMData<float> & data){
  if(data._quadruped->_robotType == RobotType::MINI_CHEETAH){
    _body_height = heightoffCoM + (float)data.userParameters->bo_height; // JUMPINGUNCOM
  }else if(data._quadruped->_robotType == RobotType::CHEETAH_3){
    _body_height = 0.45;
  }else{
    assert(false);
  }
  float x_vel_cmd, y_vel_cmd;
  float filter(0.9);
  if(data.controlParameters->use_rc){
    const rc_control_settings* rc_cmd = data._desiredStateCommand->rcCommand;
    data.userParameters->cmpc_gait = rc_cmd->variable[0];
    _yaw_turn_rate = -rc_cmd->omega_des[2];
    x_vel_cmd = rc_cmd->v_des[0];
    y_vel_cmd = rc_cmd->v_des[1] * 0.5;
    _body_height += rc_cmd->height_variation * 0.08;
  }else{
    _yaw_turn_rate = data._desiredStateCommand->rightAnalogStick[0];
    _roll_turn_rate = data._desiredStateCommand->rightAnalogStick[1];
    x_vel_cmd = data._desiredStateCommand->leftAnalogStick[1];
    y_vel_cmd = data._desiredStateCommand->leftAnalogStick[0];
    x_vel_cmd = x_vel_cmd-fminf(fmaxf(x_vel_cmd,-0.1f),0.1f);
    y_vel_cmd = y_vel_cmd-fminf(fmaxf(y_vel_cmd,-0.1f),0.1f);
    _yaw_turn_rate = _yaw_turn_rate - fminf(fmaxf(_yaw_turn_rate,-0.1f),0.1f);
    _yaw_turn_rate = (float)(data.userParameters->yawmaxfs)*_yaw_turn_rate;//
    _roll_turn_rate = _roll_turn_rate - fminf(fmaxf(_roll_turn_rate,-0.1f),0.1f);
    _roll_turn_rate = (float)(data.userParameters->rollmaxfs)*_roll_turn_rate;//
  }
  _x_vel_des = _x_vel_des*(1-filter) + (float)(data.userParameters->xmaxfs)*x_vel_cmd*filter;//
  _y_vel_des = _y_vel_des*(1-filter) + (float)(data.userParameters->ymaxfs)*y_vel_cmd*filter;//
  //_yaw_des = data._stateEstimator->getResult().rpy[2] + dt * iterationsBetweenMPC * _yaw_turn_rate;
  //_roll_des = data._stateEstimator->getResult().rpy[0] + dt * iterationsBetweenMPC * _roll_turn_rate;
  //_yaw_des = fmaxf(fminf(_yaw_des + dt * _yaw_turn_rate,0.74f),-0.74f);
  //_roll_des = fmaxf(fminf(_roll_des + dt * _roll_turn_rate,0.74f),-0.74f);
  _yaw_des = _yaw_des + dt * _yaw_turn_rate;
  _roll_des = _roll_des + dt * _roll_turn_rate;
  //cout << _roll_des << endl;//cout << _yaw_des << endl;//_roll_des = 0.;
  _pitch_des = (float)(data.userParameters->pitchdesin);//
  yawrotjump = (float)data.userParameters->ywrtjmp;
  //_roll_turn_rate = 0.;
}

template<>
void ConvexMPCLocomotion::run(ControlFSMData<float>& data) {
  bool omniMode = false;
  _SetupCommand(data);
  gaitNumber = data.userParameters->cmpc_gait;
  if(gaitNumber >= 10) {
    gaitNumber -= 10;
    omniMode = true;
  }
  auto& seResult = data._stateEstimator->getResult();
  // Check if transition to standing//UGHHHH
  if(((gaitNumber == 4) && current_gait != 4) || firstRun){//consider recast//CAREFUL WHEN CHAINIG MOTIONS
    //_yaw_des
    world_position_desiredWBC[0] = seResult.position[0];
    world_position_desiredWBC[1] = seResult.position[1];
    world_position_desiredWBC[2] = _body_height;
    //_roll_des = seResult.rpy[0];//RECONSIDER THIS IF GENERALIZATION OF BALANCE IS IMPLEMENTED (invert wrt Rpla)
    //_pitch_des = seResult.rpy[1];
    _yaw_des = seResult.rpy[2];//TEST?//CORRECTTTTTT UNDER FLIGHT//new code verify
    for(int i = 0; i < 4; i++){
      pFoot_des[i] = pFoot[i];//WELL YES BUT ACTUALLY NO, carefull with heights
      vFoot_des[i].Zero();
      aFoot_des[i].Zero();
    }
    /*for(int i = 0; i < 4; i++){//TEST LATER
      footSwingTrajectories[i].computeSwingTrajectoryBezier(1,1);//second argument?
      pFoot_des[i] = footSwingTrajectories[i].getPosition();//smoothize it??? consider Pf may be??
      vFoot_des[i].Zero();
      aFoot_des[i].Zero();
    } */
  }

  if((float)data.userParameters->noncopjump!=prevnoncopjump){
    prevnoncopjump=(float)data.userParameters->noncopjump;
    if(prevnoncopjump>0.5){
      aplatarget=&apla2;
      Rplatarget[0]=&Rpla2[0];
      Rplatarget[1]=&Rpla2[1];
      Rplatarget[2]=&Rpla2[2];
      Rplatarget[3]=&Rpla2[3];}
    else{
      aplatarget=&apla;
      Rplatarget[0]=&Rpla[0];
      Rplatarget[1]=&Rpla[1];
      Rplatarget[2]=&Rpla[2];
      Rplatarget[3]=&Rpla[3];}
  }
  if((float)data.userParameters->ifleap!=isjleap){//VERIFY
    isjleap=(float)data.userParameters->ifleap;
  }
  if((float)data.userParameters->noncopjump!=isnoncopj){
    isnoncopj=(float)data.userParameters->noncopjump;
  }
  /*if((float)data.userParameters->bckflstab!=isbckflstab){//VERIFY PRINT
    isbckflstab=(float)data.userParameters->bckflstab;
  }*/

  previous_gait = current_gait;
  Gait* gait = &trotting;//9
  if((int)data.userParameters->jumppg!=isparkjump || (float)data.userParameters->bckflstab!=isbckflstab || (float)data.userParameters->sideflstab!=issideflstab){
    isparkjump=(int)data.userParameters->jumppg;
    isbckflstab=(float)data.userParameters->bckflstab;
    issideflstab=(float)data.userParameters->sideflstab;
    if(isparkjump>0.5){jumping=&jumpingpg;}
    else if(isbckflstab<0.5f){jumping=&jumpingodg;}
    else if(issideflstab<0.5f){jumping=&jumpingbfspg;}
    else{jumping=&jumpingsfspg;}//
    jumping->recomputeVals(jumping->jumpConfGait, Np1, Np2, Np3);
  }
  ////cout<<"HI1";
  Gait* nxtgait = &trotting;//9
  if(gaitNumber == 1)       gait = &bounding;
  else if(gaitNumber == 2)  gait = &pronking;
  else if(gaitNumber == 3)  gait = &walking;
  else if(gaitNumber == 4)  gait = &standing;
  else if(gaitNumber == 5)  gait = &trotRunning;
  else if(gaitNumber == 6)  gaitNumber=9;//gait = jumping;//almost never using this
  else if(gaitNumber == 7)  gait = &galloping;
  else if(gaitNumber == 8)  gait = &pacing;
  else if(gaitNumber == 0)  gait = &walkingpg;
  nxtgait = gait;
  current_gait = gaitNumber;
  
  jumping->setIterations(27/2, iterationCounter);
  if(Np1!=(float)data.userParameters->Npar1){
    Np1=(float)data.userParameters->Npar1;
    jumping->set_nIterations(Np1);
    jumping->recomputeVals(jumping->jumpConfGait, Np1, Np2, Np3);
  }
  if(Np2!=(float)data.userParameters->Npar2){
    Np2=(float)data.userParameters->Npar2;
    jumping->recomputeVals(jumping->jumpConfGait, Np1, Np2, Np3);
  }
  if(Np3!=(float)data.userParameters->Npar3){
    Np3=(float)data.userParameters->Npar3;
    jumping->setdurations(Np3);//modify
    jumping->recomputeVals(jumping->jumpConfGait, Np1, Np2, Np3);
  }
  jump_state.trigger_pressed(data._desiredStateCommand->trigger_pressed);//modif7

  gait->setIterations(iterationsBetweenMPC, iterationCounter);
  nxtgait->setIterations(iterationsBetweenMPC, iterationCounter);
  bool earlyjumpstop=false;
  bool PhasebasedTD=false;
  bool m4lgcont=false;
  ////cout<<"HI3";
  if(currently_jumping){
    if(flyingjump){
      float heightoffCoM_fest=(*aplatarget)(0)+(*aplatarget)(1)*seResult.position[0]+(*aplatarget)(2)*seResult.position[1];
      if (data.userParameters->erlylnd>0.5){
        earlyjumpstop=(seResult.position[2]-heightoffCoM_fest<(float)data.userParameters->bo_height) && (seResult.vWorld[2]<0);//CAREFUL WITH DOWNSLOPS
      }
      if (data.userParameters->ltelnd>0.5){
        float dheightoffCoM_fest=(*aplatarget)(1)*seResult.vWorld[0]+(*aplatarget)(2)*seResult.vWorld[1];
        int segcj=jumping->getCurrentGaitPhase();
        deltestlndpre = deltestlnd;
        float deltzovg=(seResult.position[2]-heightoffCoM_fest-(float)data.userParameters->bo_height)/9.81f;
        float ddeltzovg=(seResult.vWorld[2]-dheightoffCoM_fest)/9.81f;
        deltestlnd=ddeltzovg+sqrt((ddeltzovg*ddeltzovg)+2*deltzovg);//deltestlnd/dtMPC//OBTAIN FROM MPC LOL // CAREFUL WHEN DIFF INITIAL FINAL HEIGHT
        float _phaseLO=jumping->get_phaseLO();//define
        float _phasela=jumping->get_phasela();//define 1 FOR LIM _phase->1-
        tFlcount+=1;
        if (deltestlnd0==-10){deltestlnd0=deltestlnd;tFlcount0=tFlcount;cout<<"deltestlnd0:"<<deltestlnd0<<endl;}//APEX?????
        else if((deltestlnd0+tFlcount0*dt)<(deltestlnd+tFlcount*dt)){deltestlnd0=deltestlnd;tFlcount0=tFlcount;}
        //if (deltestlnd0==-10){deltestlnd0=deltestlnd;cout<<"deltestlnd0:"<<deltestlnd0<<endl;}//WRITE PROBABALISTIC FUNCTION HERE, ENHANCE BLEDT'S PAPER WITH CENTROIDAL DYN AND 0 CONTACT.
        int jNp1=jumping->getHoriz();//SCALED
        float _phase1=0;
        int _iteration1=0;
        float _relphase1 = (tFlcount*dt)/(deltestlnd0+tFlcount0*dt);
        _relphase1=fmin(fmax(_relphase1,0),1);
        if(deltestlndpre>deltestlnd){
          counttrust++;
          _relphase1=fminf(_relphase1,0.9999f);//correct//_phase1=fminf(_phase1,0.9f);//CORRECT EQUIVALENCE
        }else if(counttrust>40){//0.08msafter counter????????
          _relphase1=1.0001;//_phase1=1.001f;
        }
        if (_relphase1<1){
          _phase1=fmaxf(_phaseLO+(_phasela-_phaseLO)*_relphase1,_phaseLO);//SCALED
          _iteration1=(max((int)(_phase1*jNp1),(int)round(_phaseLO*jNp1)))%jNp1;
          _phase1=_phase1-(int)_phase1;
          jumping->setIterationssb(_iteration1,_phase1);
        }
        else{
          _relphase1=fminf(_relphase1,0.9999f);
          _phase1=fmaxf(_phaseLO+(_phasela-_phaseLO)*_relphase1,_phaseLO);//SCALED
          _iteration1=(max((int)(_phase1*jNp1),(int)round(_phaseLO*jNp1)))%jNp1;
          _phase1=_phase1-(int)_phase1;
          jumping->setIterationssb(_iteration1,_phase1);
          _relphase1=1.0001;

          _phase1=_phasela;
          _iteration1=((int)(_phasela*jNp1))%jNp1;
          //jumping->setIterationssb(_iteration1,_phase1);
          jumping->setItCorrected(27/2, iterationCounter);//// ///cout<<"Hijumpearlystop"<<endl;
          flyingjump=false;deltestlnd0=-10;tFlcount=0;tFlcount0=0;deltestlndpre=0;deltestlnd=0;counttrust=0;
        }//CONSIDER RECAST AS HYDRIMDE VNACH
       printf("counttrust: %i",counttrust);cout<<"tlndclosed: "<<deltestlnd<<" tlndphasebased: "<<(Np1-segcj)*dtMPC<<" segcj: "<<segcj;//<<(seResult.position[2]-heightoffCoM_fest-(float)data.userParameters->bo_height);
       cout<<" _iteration1: "<<_iteration1<<" _phase1: "<<_phase1<<" _relphase1: "<<_relphase1<<endl;
        //float _phase1=fmaxf(1-deltestlnd/(Np1*dtMPC),_phaseLO);
        //float _phase1=fmaxf(1-deltestlnd*(1-_phaseLO)/(deltestlnd0),_phaseLO);//SCALED
      }////cout<<"HI5";
      PhasebasedTD=jumping->getFlightState();
    }
    if(earlyjumpstop&&PhasebasedTD){
      jumping->setItCorrected(27/2, iterationCounter);
     //// ///cout<<"Hijumpearlystop"<<endl;
      flyingjump=false;deltestlnd0=-10;tFlcount=0;tFlcount0=0;deltestlndpre=0;deltestlnd=0;counttrust=0;
    }
  }else{
    cout<<"gait->getCurrentGaitPhase()"<<gait->getCurrentGaitPhase();
    if(jump_state.jump_pending){
      if(!pre0init){
        if(gait->getCurrentGaitPhase()!=0){
          pre0init=true;
        }
      }else if(gait->getCurrentGaitPhase()==0){
        cout<<"HI6a"<<endl;
        m4lgcont=true;
        jumping->setItCorrectedStart(iterationsBetweenMPC,iterationCounter);
      }
    }
  }
  if(jump_state.should_jump(jumping->getCurrentGaitPhase(),m4lgcont)){
    gait = jumping;
    current_gait = 6;
    recompute_timing(27/2);
    _body_height = heightoffCoM+(float)data.userParameters->bo_height;//bo_hei_jum
    currently_jumping = true;
    if(isbckflstab>.5){
      chSfOri = false;
    }else{
      chSfOri = true;
    }
  }else{
    chSfOri = true;
    recompute_timing(default_iterations_between_mpc);//switch beetwen iBMPC and 27/2 for jumping...
    currently_jumping = false;flyingjump=false;deltestlnd0=-10;tFlcount=0;tFlcount0=0;deltestlndpre=0;deltestlnd=0;counttrust=0;
    //jumping->gaitphascorr=0;
    jumping->countflystage=0;
   //// ///cout<<"HIRIP";
  }
  if(_body_height < 0.02) _body_height = 0.29;
  if (current_gait!=previous_gait){
    bswitchgait = true;
    if (previous_gait==6){//DEBUG/CONSIDER USE A LITTLE OF FORCE THERE
    pre0init = false;
    world_position_desiredWBC[0] = seResult.position[0];
    world_position_desiredWBC[1] = seResult.position[1];
    world_position_desiredWBC[2] = _body_height;
    _yaw_des = seResult.rpy[2];
      for(u16 foot=0;foot<4;foot++){
        firstSwing[foot]=true;
  }}}
  ////cout<<"HI7";
  int* mpcTable = gait->getMpcTable();
  //print_named_arraytr("MPCTablepre",mpcTable,4,jumping->getFliMacNH(2));
 //// ///printf("jumping->getCurrentGaitPhase() %i",jumping->getCurrentGaitPhase());
  ////cout<<"HI8";
  int horlen;
  if (gait == jumping){
    int gaitHormcurp = gait->getHoriz()-gait->getCurrentGaitPhase();
    int nxtgaithor = nxtgait->getHoriz();
    horlen=fmin(fmax(gaitHormcurp,nxtgaithor),11.0f);
    if(horlen-gaitHormcurp>0){//for(int i = 0; i < gaitcurp; i++)for(int j = 0; j < 4; j++)mpcTable[4*(gaitHor-gaitcurp+i)+j]=nxtmpcTable[4*((nxtintpick+i)% nxtgaithor)+j];
      int* nxtmpcTable = nxtgait->getMpcTable();
      int nxtintpick = (-(nxtgait->getCurrentGaitPhase())) % nxtgaithor;// //++ +1?
      nxtintpick = (nxtintpick+nxtgaithor);
      for(int i = 0; i < horlen-gaitHormcurp; i++){//min(gaitcurp,gaitcurp+horlen-gaitHor)
        for(int j = 0; j < 4; j++)
          mpcTable[4*(gaitHormcurp+i)+j]=nxtmpcTable[4*((nxtintpick+i)% nxtgaithor)+j];
      }
    }
  }
  else{
    horlen=fmin(gait->getHoriz(),11.0f);
    if(previous_gait == 6){
    //gait->gaitphascorr = (gait->gaitphascorr+(nxtgait->getHoriz()-(nxtgait->getCurrentGaitPhase())+1)*iterationsBetweenMPC) % (nxtgait->getHoriz()*iterationsBetweenMPC);// //
    //gait->setIterations(iterationsBetweenMPC, iterationCounter);
    gait->setItCorrectedStart(iterationsBetweenMPC, iterationCounter);//Start?
    mpcTable = gait->getMpcTable();
  }}
  //printf("gaitgaitphascorr: %i \n",gait->gaitphascorr);
  //printf("jumpinggaitphascorr: %i \n",jumping->gaitphascorr);
  if(gait==jumping && jumping->getFliMacNH(2)==jumping->getCurrentGaitPhase()){//CAREFUL WHEN PRKR JUMP IS ONLY 1 JUMP
    world_position_desiredWBC[0] = seResult.position[0];
    world_position_desiredWBC[1] = seResult.position[1];
    world_position_desiredWBC[2] = _body_height;
    _yaw_des = seResult.rpy[2];
  }
  //print_named_arraytr("MPCTable",mpcTable,4,jumping->getFliMacNH(1)); 
  ////cout<<"HI9";
  v_des_robot << _x_vel_des, _y_vel_des, 0;//v_des_worldWBC = omniMode ? v_des_robot : seResult.rBody.transpose() * v_des_robot;
  if (gait!=jumping){
    v_des_worldWBC = omniMode ? v_des_robot : coordinateRotation(CoordinateAxis::Z, -seResult.rpy(2)) * v_des_robot;//ASSUMING YAW SOLVES THIS
  }/*else{
    v_des_robot = coordinateRotation(CoordinateAxis::Z, seResult.rpy(2)) * v_des_worldWBC;
  }*/ //yes for chained

  for(int i = 0; i < 4; i++) {
    pFoot[i] = seResult.position + seResult.rBody.transpose() * (data._quadruped->getHipLocation(i) + data._legController->datas[i].p);
  }
  if(gait != &standing) {
    world_position_desiredWBC += dt * Vec3<float>(v_des_worldWBC[0], v_des_worldWBC[1],0);//MODIFY-CHECK RUNTIME//CHANGE FOR SLOPES//ADD BOUNDS FOR WBC
    const float max_pos_error = .08;
    float xStart = world_position_desiredWBC[0];
    float yStart = world_position_desiredWBC[1];
    if(xStart - seResult.position[0] > max_pos_error){
      xStart = seResult.position[0] + max_pos_error;//v_des_world[0]=max_pos_error/dt;
    }
    if(seResult.position[0] - xStart > max_pos_error){
      xStart = seResult.position[0] - max_pos_error;
    }
    if(yStart - seResult.position[1] > max_pos_error){
      yStart = seResult.position[1] + max_pos_error;
    }
    if(seResult.position[1] - yStart > max_pos_error){
      yStart = seResult.position[1] - max_pos_error;
    }
    world_position_desiredWBC[0] = xStart;
    world_position_desiredWBC[1] = yStart;
  }else{
    world_position_desiredWBC[0] = (pFoot_des[0](0)+pFoot_des[1](0)+pFoot_des[2](0)+pFoot_des[3](0))/4;
    world_position_desiredWBC[1] = (pFoot_des[0](1)+pFoot_des[1](1)+pFoot_des[2](1)+pFoot_des[3](1))/4;
    world_position_desiredWBC[2] = _body_height;
  }
  if(firstRun){
    world_position_desiredWBC[0] = seResult.position[0];
    world_position_desiredWBC[1] = seResult.position[1];
    world_position_desiredWBC[2] = _body_height;
    for(int i = 0; i < 4; i++){
      footSwingTrajectories[i].setHeight(0.05);
      footSwingTrajectories[i].setInitialPosition(pFoot[i]);footSwingTrajectoriesInit[i]=pFoot[i];
      footSwingTrajectories[i].setFinalPosition(pFoot[i]);
      pFoot_des[i] = pFoot[i];
      vFoot_des[i].Zero();
      aFoot_des[i].Zero();
    }
    firstRun = false;
  }////cout<<"HI10";
  if (gait==&walkingpg){
    is_slidcntct = walkingpg.getis_slidcntct();;// END of WBC Update
  }else{
    is_slidcntct = Vec4<float>(0,0,0,0);// END of WBC Update
  }
  // foot placement
  for(int l = 0; l < 4; l++)
    swingTimes[l] = gait->getCurrentSwingTime(dtMPC, l);
  //cout<<"swingTimes: \t"<<swingTimes[0]<<", "<<swingTimes[1]<<", "<<swingTimes[2]<<", "<<swingTimes[3]<<" \t";
  float side_sign[4] = {-1, 1, -1, 1};//.....................................?
  float interleave_y[4] = {-0.08, 0.08, 0.02, -0.02};
  float interleave_gain = -0.2;//float interleave_gain = -0.13;
  float v_abs = std::fabs(v_des_robot[0]);//float v_abs = std::fabs(seResult.vBody[0]);
  //cout<<"stance_time: ";
  if(gait==jumping && (float)data.userParameters->ywrtjmp>0.05){//VERIFY STUFF
    _yaw_turn_rate+=seResult.omegaBody(2);//drpydesWBC(2)
    /*if(gait==&jumpingpg && (jumping->getFliMacNH(1)+1)>jumping->getCurrentGaitPhase()){//CAREFUL WHEN PRKR JUMP IS ONLY 1 JUMP
      _yaw_turn_rate+=seResult.omegaBody(2)*0.5;//drpydesWBC(2)
    }*/
    /////cout<<"{_yaw_turn_rate: "<<_yaw_turn_rate<<"}";
  }
  for(int i = 0; i < 4; i++){
    if(firstSwing[i]) {//footSwingTrajectories[i].setHeight(.05);}
      swingTimeRemaining[i] = swingTimes[i];
    } else {
      swingTimeRemaining[i] -= dt;
    }
    Vec3<float> offset(0, side_sign[i] * (.065+data.userParameters->shiftleg), 0);//VERIFY
    if (gait==&walkingpg){//&& is_slidcntct(i)>0.5 ADD
      float shlShf=data.userParameters->slidshift;
      if(walkingpg.getCurrentHybridMode()<4){
        if(i==0||i==3){//02
          offset(0)=-shlShf*((walkingpg.getCurrentHybridMode()+1)/2)+shlShf;
          /////cout<<"offset("<<i<<"): "<<offset(0)<<endl;
        }else{
          offset(0)=+shlShf;//-shlShf*(walkingpg.getCurrentHybridMode())
          /////cout<<"offset("<<i<<"): "<<offset(0)<<endl;
        }
      }else{
        if(i==1||i==2){//13
          offset(0)=-shlShf*((walkingpg.getCurrentHybridMode()-4+1)/2)+shlShf;
          /////cout<<"offset("<<i<<"): "<<offset(0)<<endl;
        }else{
          offset(0)=shlShf;//-shlShf*(walkingpg.getCurrentHybridMode()-4)+
          /////cout<<"offset("<<i<<"): "<<offset(0)<<endl;
        }
      }
    }
    Vec3<float> pRobotFrame = (data._quadruped->getHipLocation(i) + offset);
    pRobotFrame[1] += interleave_y[i] * v_abs * interleave_gain;
    float stance_time = gait->getCurrentStanceTime(dtMPC, i); //for heuristics //cout<<stance_time<<", ";
    /////cout<<"stance_time["<<i<<"]"<<stance_time;
    Vec3<float> pYawCorrected = coordinateRotation(CoordinateAxis::Z, -_yaw_turn_rate* stance_time/ 2) * pRobotFrame;
    Vec3<float> Pf = seResult.position + seResult.rBody.transpose() * (pYawCorrected + v_des_robot * swingTimeRemaining[i]);//this produces variations//CONSIDER RECAST OF v_des_robot IN MPC//+ seResult.vWorld * swingTimeRemaining[i];//float p_rel_max = 0.35f;
    //vtrajCorr[i] = seResult.omegaWorld.cross(Pf - seResult.position)*0;
    float p_rel_max = 0.3f*1.4142f;//0.3f;
    float pfx_rel = seResult.vWorld[0] * (.5 + _parameters->cmpc_bonus_swing) * stance_time + .03f*(seResult.vWorld[0]-v_des_worldWBC[0]) +
      (0.5f*seResult.position[2]/9.81f) * (seResult.vWorld[1]*_yaw_turn_rate);//.1f
    float pfy_rel = seResult.vWorld[1] * .5 * stance_time + .03f*(seResult.vWorld[1]-v_des_worldWBC[1]) +
      (0.5f*seResult.position[2]/9.81f) * (-seResult.vWorld[0]*_yaw_turn_rate);//.1f
    if(gait==jumping && isjleap>0.5){ //(float)data.userParameters->ifleap
      if(isnoncopj){//add lat cond
        pfx_rel+=.2f*(seResult.vWorld[0]);
        pfy_rel+=.2f*(seResult.vWorld[1]);
      }else if(gait==&jumpingpg && (jumping->getFliMacNH(1)+1)>jumping->getCurrentGaitPhase()){
        pfx_rel+=.0f*(seResult.vWorld[0]);//not helping in double jumps
        pfy_rel+=.0f*(seResult.vWorld[1]);
      }else {
        pfx_rel+=.03f*(seResult.vWorld[0]);//not helping in double jumps
        pfy_rel+=.03f*(seResult.vWorld[1]);
      }
    }
    float pf_norm=sqrt(pfx_rel*pfx_rel+pfy_rel*pfy_rel);
    if (pf_norm>p_rel_max){pfx_rel=pfx_rel*p_rel_max/pf_norm;pfy_rel=pfy_rel*p_rel_max/pf_norm;}
    Pf[0] +=  pfx_rel;
    Pf[1] +=  pfy_rel;
    ////cout<<"HI11";
    //Vec3<float> plegfut= seResult.rBody * Vec3<float>(Pf[0]-seResult.position[0],Pf[1]-seResult.position[1],0)-data._quadruped->getHipLocation(i);//coordinateRotation(CoordinateAxis::Z, -seResult.rpy(2)) * Vec3<float>(Pf[0]-seResult.position[0],Pf[1]-seResult.position[1],0);
    //pFoot[i] = seResult.position + seResult.rBody.transpose() * (data._quadruped->getHipLocation(i) + data._legController->datas[i].p);
    if(gait==jumping){Pf[2] = (*aplatarget)(0)+(*aplatarget)(1)*Pf[0]+(*aplatarget)(2)*Pf[1]-0.29+seResult.position[2]-heightoffCoM_prev;//cout<<"Pf2"<<Pf[2]<<endl;//_body_height
    }else{Pf[2] = apla(0)+apla(1)*Pf[0]+apla(2)*Pf[1]-0.003;}//Pf[2] = 0.0;
    footSwingTrajectories[i].setFinalPosition(Pf);//cout << "pfx_rel: "<<pfx_rel<< endl;
    if(gait != jumping){
      heightoffCoM_prev=heightoffCoM;
      float dvzf = apla(1)*(Pf[0]-footSwingTrajectoriesInit[i](0))+apla(2)*(Pf[1]-footSwingTrajectoriesInit[i](1));
      float shfi = fminf(fmaxf(0.06f+fmaxf(fminf(dvzf,0.0f),-0.06f),dvzf),0.12f);
      if (gait==&walkingpg){shfi=(shfi+0.06)*(1-is_slidcntct(i))+Pf[2]*is_slidcntct(i);}
      footSwingTrajectories[i].setHeight(shfi);//cout<<","<<shfi;
    }else{
      footSwingTrajectories[i].setHeight(seResult.position[2]-heightoffCoM_prev-0.18);//footSwingTrajectories[i].setHeight(seResult.position[2]-heightoffCoM-0.2);//BAD WAY: BEETER TRACKING//
    }
  }
  //calc gait//load LCM leg swing gains
  /*Kp << 700, 0, 0,
     0, 700, 0,
     0, 0, 150; */
  Kp << 2000, 0, 0,
     0, 2000, 0,
     0, 0, 400;
  Kp_stance = 0*Kp;
  Kd << 7, 0, 0,
     0, 7, 0,
     0, 0, 7;
  Kd_stance = Kd;// gait
  contact_state = gait->getContactState();
  Vec4<float> swingStates = gait->getSwingState();
  //cout<<"swingStates: \t"<<swingStates.transpose()<<"\t contact_state: \t"<<contact_state.transpose()<<endl;
  //////printf("jumping->get_current_gait_phase(): %i\n", jumping->getCurrentGaitPhase());
  /*for(int i=0;i<4;i++){
    //pFoot[i] = seResult.position + seResult.rBody.transpose() * (data._quadruped->getHipLocation(i) + data._legController->datas[i].p);
    Vec3<float> vFoot = seResult.vWorld + seResult.rBody.transpose() * (data._legController->datas[i].v)+seResult.omegaWorld.cross(seResult.rBody.transpose() * (data._quadruped->getHipLocation(i) + data._legController->datas[i].p));
    //cout<<"dataslegcontrollerv"<<i<<data._legController->datas[i].v<<endl;
    cout<<"dataslegcontrollerv"<<i<<vFoot<<endl;
  }*/
  Timer updMPCTimer;
  updateMPCIfNeeded(gait, mpcTable, data, horlen);//15.0f lat frontal leap//horlen=fmin(fmax(gait->getHoriz()-gaitcurp,nxtgait->getHoriz()),11.0f);
  //updateMPCIfNeeded(mpcTable, data, omniMode,gait->getHoriz()-fmin(gaitcurp,gait->getHoriz()-nxtgait->getHoriz()));
  if (updMPCTimer.getMs()>0.02){//cout << "\nRpla: " << Rpla[0]<<endl;
   printf("TIME update MPC If Needed: %.3f\n", updMPCTimer.getMs());
  }
  //  StateEstimator* se = hw_i->state_estimator;
  pestCoMx=0;pestCoMy=0;
  for(s16 foot = 0; foot < 4; foot++){
      float trst = 1;
      float phse = fmin(contact_state(foot), float(1));
      float trust_wdw = 0.2;
      if (phse < trust_wdw) {trst = phse / trust_wdw;}
      else if (phse > (float(1) - trust_wdw)) {trst = (float(1) - phse) / trust_wdw;}
      if(trst==float(1) && gait != jumping){//SWING??
          Wpla(foot,1)=pFoot[foot](0);
          Wpla(foot,2)=pFoot[foot](1);
          zfeet(foot)=pFoot[foot](2);
          pseudoInverse(Wpla,0.001,Wplainv);
          float coef_filsat=fminf(fmaxf((float)data.userParameters->coef_fil,0.01),0.99);
          apla=(1-coef_filsat)*apla+coef_filsat*Wplainv*zfeet;//cout<<"apla: "<<apla<<endl;//cout << "zfeet: \n" << zfeet;//cout << "\nheighoffCoM: " << heightoffCoM;
      }
    pestCoMx+=(pFoot[foot](0))/4;
    pestCoMy+=(pFoot[foot](1))/4;
  }
  ////cout<<"HI13";
  if(data.userParameters->apla2align>0.5){
    float s = std::sin(seResult.rpy(2));
    float c = std::cos(seResult.rpy(2));
    Matrix<float,2,2> Rmy;Rmy<<c, s,-s, c;
    Matrix<float,1,2>apla2con=Vec2<float>(apla2b(1),apla2b(2)).transpose()*Rmy;
    apla2(0)=apla2b(0)-apla2con*Vec2<float>(seResult.position[0],seResult.position[1]);
    apla2(1)=apla2con(0);
    apla2(2)=apla2con(1);
    Rpla2[0]=aplatoRpla(apla2);Rpla2[1]=Rpla2[0];Rpla2[2]=Rpla2[0];Rpla2[3]=Rpla2[0];
  }
  if(data.userParameters->det_terrn>0.5){
    Rpla[0]=aplatoRpla(apla);Rpla[1]=Rpla[0];Rpla[2]=Rpla[0];Rpla[3]=Rpla[0];
    heightoffCoM=apla(0)+apla(1)*pestCoMx+apla(2)*pestCoMy;
  }
  else{/*Vec3<float> RPYPlane1(0.8726,0,0);Vec3<float> RPYPlane2(-0.8726,0,0);Rpla[0]=rpyToRotMat(RPYPlane1);Rpla[1]=rpyToRotMat(RPYPlane2);Rpla[2]=rpyToRotMat(RPYPlane1);Rpla[3]=rpyToRotMat(RPYPlane2);*/ //Vshaped
    Rpla[0].Identity();Rpla[1].Identity();Rpla[2].Identity();Rpla[3].Identity();
    heightoffCoM=0;
  }
  ////cout<<"HI14";
#ifdef DRAW_DEBUG_PATH
  auto* trajectoryDebug = data.visualizationData->addPath();
    trajectoryDebug->num_points = 10;
    trajectoryDebug->color = {0.2, 0.2, 0.7, 0.5};
    for(int i = 0; i < gait->getHoriz(); i++) {
      trajectoryDebug->position[i][0] = trajAll[13*i + 4];
      trajectoryDebug->position[i][1] = trajAll[13*i + 5];
      trajectoryDebug->position[i][2] = trajAll[13*i + 6];
      auto* ball = data.visualizationData->addSphere();
      ball->radius = 0.01;
      ball->position = trajectoryDebug->position[i];
      ball->color = {0.2, 0.2, 0.7, 0.5};
    }
    update_trajh(trajPrev,10);//10?gait->getHoriz()
    auto* trajectoryPrevDebug = data.visualizationData->addPath();
    trajectoryPrevDebug->num_points = 10;//10?gait->getHoriz()
    trajectoryPrevDebug->color = {0.9, 0.2, 0.1, 0.5};
    for(int i = 0; i < 10; i++) {//gait->getHoriz()
      trajectoryPrevDebug->position[i][0] = trajPrev[13*i + 4];
      trajectoryPrevDebug->position[i][1] = trajPrev[13*i + 5];
      trajectoryPrevDebug->position[i][2] = trajPrev[13*i + 6];
      //cout<<"trajectoryPrevDebug->position{"<<i<<"}"<<trajectoryPrevDebug->position[i]<<endl;
      auto* ball = data.visualizationData->addSphere();
      ball->radius = 0.01;
      ball->position = trajectoryPrevDebug->position[i];
      ball->color = {0.9, 0.2, 0.1, 0.5};
    }
auto* trajheurDebug = data.visualizationData->addPath();
    trajheurDebug->num_points = jumping->getHoriz();
    trajheurDebug->color = {0.9, 0.8, 0.1, 0.5};
    for(int i = 0; i < jumping->getHoriz(); i++) {
      trajheurDebug->position[i][0] = trajHeurJump[13*i + 4];
      trajheurDebug->position[i][1] = trajHeurJump[13*i + 5];
      trajheurDebug->position[i][2] = trajHeurJump[13*i + 6];
      auto* ball = data.visualizationData->addSphere();
      ball->radius = 0.01;
      ball->position = trajheurDebug->position[i];
      ball->color = {0.1, 0.8, 0.1, 0.5};//cout<<"debugtraj"<<i;
    }
#endif
  for(int foot = 0; foot < 4; foot++){
    float swingState = swingStates[foot];//
    if(swingState > 0){// foot is in swing
      if(firstSwing[foot]){
        firstSwing[foot] = false;
        footSwingTrajectories[foot].setInitialPosition(pFoot[foot]);//TRIGGER THIS AT THE END OF JUMPING/BEG OF NEW BEH
        footSwingTrajectoriesInit[foot]=pFoot[foot];
      }
#ifdef DRAW_DEBUG_SWINGS
      auto* debugPath = data.visualizationData->addPath();
      if(debugPath) {
        debugPath->num_points = 100;
        debugPath->color = {0.2,1,0.2,0.5};
        float step = (1.f - swingState) / 100.f;
        for(int i = 0; i < 100; i++) {
          footSwingTrajectories[foot].computeSwingTrajectoryBezier(swingState + i * step, swingTimes[foot]);
          debugPath->position[i] = footSwingTrajectories[foot].getPosition();
      }}
      auto* finalSphere = data.visualizationData->addSphere();
      if(finalSphere) {
        finalSphere->position = footSwingTrajectories[foot].getPosition();
        finalSphere->radius = 0.02;finalSphere->color = {0.6, 0.6, 0.2, 0.7};
      }
      footSwingTrajectories[foot].computeSwingTrajectoryBezier(swingState, swingTimes[foot]);
      auto* actualSphere = data.visualizationData->addSphere();
      auto* goalSphere = data.visualizationData->addSphere();
      goalSphere->position = footSwingTrajectories[foot].getPosition();
      actualSphere->position = pFoot[foot];
      goalSphere->radius = 0.02;actualSphere->radius = 0.02;
      goalSphere->color = {0.2, 1, 0.2, 0.7};actualSphere->color = {0.8, 0.2, 0.2, 0.7};
#endif
      footSwingTrajectories[foot].computeSwingTrajectoryBezier(swingState, swingTimes[foot]);
      //footSwingTrajectories[foot]->updateFF(hw_i->leg_controller->leg_datas[foot].q,
      //hw_i->leg_controller->leg_datas[foot].qd, 0); // velocity dependent friction compensation todo removed
      //hw_i->leg_controller->leg_datas[foot].qd, fsm->main_control_settings.variable[2]);
      Vec3<float> pDesLeg;
      Vec3<float> vDesLeg;
      if(isbckflstab>0.5f && gait==jumping){ //From permaflight 2 leg balance, write local WBC task
        //pDesLeg<<-0.1,side_sign[foot] * (.065+0.025),-0.25;//Vec3<float> pFlow(0,side_sign[foot] * .065,-0.2);//-0.16
        float phj4l=jumping->getCurrentGaitPhase();
        if(issideflstab>0.5f){
          pDesLeg<<0.0,side_sign[foot] * (.065+0.01),-0.2f - std::abs(2*(phj4l-(Np2+Np3))/(Np1-(Np2+Np3))-1)*0.05;//Vec3<float> pFlow(0,side_sign[foot] * .065,-0.2);//-0.16
          if((foot==1||foot==3)&&(phj4l<(Np2+Np3+3))){
            pDesLeg<<0.0,side_sign[foot] * (.065+0.155),-0.25f;//Vec3<float> pFlow(0,side_sign[foot] * .065,-0.2);//-0.16
          }
          vDesLeg.setZero();//seResult.rBody * (vDesFootWorld - seResult.vWorld)+seResult.omegaWorld.cross(pDesLeg);//   Derivatives?omegaWorld//verify
          pFoot_des[foot] = seResult.position + seResult.rBody.transpose() * (data._quadruped->getHipLocation(foot) + pDesLeg);
          vFoot_des[foot] = seResult.vWorld + 1.9*seResult.omegaWorld.cross(pFoot_des[foot] - seResult.position);//Vec3<float> vDesFootWorld;vDesFootWorld.setZero();//footSwingTrajectories[foot].getVelocity()
          aFoot_des[foot] = 1.9*seResult.omegaWorld.cross(vFoot_des[foot] - seResult.vWorld);
          aFoot_des[foot](2)-=9.81f;//cout<<"pFoot_des"<<pFoot_des[foot];
        }else{
          pDesLeg<<-0.1,side_sign[foot] * (.065+0.025),-0.25f;//Vec3<float> pFlow(0,side_sign[foot] * .065,-0.2);//-0.16
          if((foot==2||foot==3)&&(phj4l<(Np2+Np3+3))){
            pDesLeg<<-0.3,side_sign[foot] * (.065+.045),-0.0f;//Vec3<float> pFlow(0,side_sign[foot] * .065,-0.2);//-0.16
          }
          vDesLeg.setZero();//seResult.rBody * (vDesFootWorld - seResult.vWorld)+seResult.omegaWorld.cross(pDesLeg);//   Derivatives?omegaWorld//verify
          pFoot_des[foot] = seResult.position + seResult.rBody.transpose() * (data._quadruped->getHipLocation(foot) + pDesLeg);
          vFoot_des[foot] = seResult.vWorld + seResult.omegaWorld.cross(pFoot_des[foot] - seResult.position);//Vec3<float> vDesFootWorld;vDesFootWorld.setZero();//footSwingTrajectories[foot].getVelocity()
          aFoot_des[foot] = seResult.omegaWorld.cross(vFoot_des[foot] - seResult.vWorld);
          aFoot_des[foot](2)-=9.81f;//cout<<"pFoot_des"<<pFoot_des[foot];
        }
#ifdef DRAW_DEBUG_SWINGS
  auto* pfdesSphere = data.visualizationData->addSphere();
  pfdesSphere->position = pFoot_des[foot];
  pfdesSphere->radius = 0.02;pfdesSphere->color = {1.0, 0.0, 1.0, 0.7};
#endif
      }else{
        Vec3<float> pDesFootWorld = footSwingTrajectories[foot].getPosition();
        Vec3<float> vDesFootWorld = footSwingTrajectories[foot].getVelocity();//+vtrajCorr[foot]
        if (gait==jumping){
          vDesFootWorld(2)+=seResult.vWorld[2];
        }
        pDesLeg = seResult.rBody * (pDesFootWorld - seResult.position)- data._quadruped->getHipLocation(foot);//Rtranspose is seResult.rBody (casts from world to body)
        vDesLeg = seResult.rBody * (vDesFootWorld - seResult.vWorld)-seResult.omegaBody.cross(pDesLeg+data._quadruped->getHipLocation(foot));//   Derivatives?omegaWorld...-omegabody?
        // Update for WBC//cout <<"pDesFootWorld: "<<pDesFootWorld<<endl;
        pFoot_des[foot] = pDesFootWorld;
        vFoot_des[foot] = vDesFootWorld;
        aFoot_des[foot] = footSwingTrajectories[foot].getAcceleration();
        if (gait==jumping){
          aFoot_des[foot](2)-=9.81f;
        }
      }


      if(!data.userParameters->use_wbc){// Update leg control command regardless of the usage of WBIC
        data._legController->commands[foot].pDes = pDesLeg;
        data._legController->commands[foot].vDes = vDesLeg;
        data._legController->commands[foot].kpCartesian = Kp;
        data._legController->commands[foot].kdCartesian = Kd;
      }}else{ //foot is in stance
      firstSwing[foot] = true;
#ifdef DRAW_DEBUG_SWINGS
      auto* actualSphere = data.visualizationData->addSphere();
      actualSphere->position = pFoot[foot];
      actualSphere->radius = 0.02;actualSphere->color = {0.2, 0.2, 0.8, 0.7};
#endif
      Vec3<float> pDesLeg, vDesLeg;
      if(gait!=&standing){
        Vec3<float> pDesFootWorld = footSwingTrajectories[foot].getPosition();
        Vec3<float> vDesFootWorld = footSwingTrajectories[foot].getVelocity();//+vtrajCorr[foot]
        pDesLeg = seResult.rBody * (pDesFootWorld - seResult.position) - data._quadruped->getHipLocation(foot);
        vDesLeg = seResult.rBody * (vDesFootWorld - seResult.vWorld)-seResult.omegaBody.cross(pDesLeg+data._quadruped->getHipLocation(foot));;//?
      //cout << "Foot " << foot << " relative velocity desired: " << vDesLeg.transpose() << "\n";
      }else{
        Vec3<float> pDesFootWorld(pFoot_des[foot]);        //cout<<"pDesFootWorld"<<foot<<" : " <<pFoot_des[foot].transpose()<<"\n";
        Vec3<float> vDesFootWorld(0,0,0);
        pDesLeg = seResult.rBody * (pDesFootWorld - seResult.position) - data._quadruped->getHipLocation(foot);//MAYBE????
        vDesLeg = seResult.rBody * (vDesFootWorld - seResult.vWorld)-seResult.omegaBody.cross(pDesLeg+data._quadruped->getHipLocation(foot));;//?
      }
      if(!data.userParameters->use_wbc){
        data._legController->commands[foot].pDes = pDesLeg;
        data._legController->commands[foot].vDes = vDesLeg;
        data._legController->commands[foot].kpCartesian = Kp_stance;
        data._legController->commands[foot].kdCartesian = Kd_stance;
        data._legController->commands[foot].forceFeedForward = f_ff[foot];
        data._legController->commands[foot].kdJoint = Mat3<float>::Identity() * 0.2;
        //footSwingTrajectories[foot]->updateFF(hw_i->leg_controller->leg_datas[foot].q,
        //hw_i->leg_controller->leg_datas[foot].qd, 0); todo removed
        //hw_i->leg_controller->leg_commands[foot].tau_ff += 0*footSwingController[foot]->getTauFF();
      }else{ // Stance foot damping
        data._legController->commands[foot].pDes = pDesLeg;
        data._legController->commands[foot].vDes = vDesLeg;
        data._legController->commands[foot].kpCartesian = 0.*Kp_stance;
        data._legController->commands[foot].kdCartesian = Kd_stance;
      }
      //cout << "Foot " << foot << " force: " << f_ff[foot].transpose() << "\n";//Fr_des[foot] = -f_ff[foot];// Update for WBC
}}////cout<<"HI15";
  data._stateEstimator->setContactPhase(contact_state);// TODO SET ONLY ONCE MAYBE?

  update_trajh(trajPrev,1);// Update For WBC //TRUST ON MPC
  if(gait==jumping){
    v_des_worldWBC<<trajPrev[10],trajPrev[11],trajPrev[12];//CHANGE TO HEUR?//OVERTRUSTING HEURISTIC?//WRONG TRAJALL IS SHIFTED 
  }
  //cout<<endl;cout<<"v_des_worldWBC: "<<v_des_worldWBC<<endl;

  pBody_des[0] = world_position_desiredWBC[0];
  pBody_des[1] = world_position_desiredWBC[1];
  //pBody_des[2] = _body_height;//world_position_desiredWBC[2];//reconsider for task0123 456
  pBody_des[2] = trajPrev[6];
  if (gait==&walkingpg){
    pBody_des[2]=_body_height;
  }
  vBody_des[0] = v_des_worldWBC[0];
  vBody_des[1] = v_des_worldWBC[1];
  vBody_des[2] = v_des_worldWBC[2];
  //cout<<"v_des_worldWBCUP: "<<v_des_worldWBC.transpose();
  aBody_des.setZero();
  //if(gait==jumping){cout<<"rpydesWBC:"<<rpydesWBC.transpose()<<"\t \t drpydesWBC:"<<drpydesWBC.transpose()<<"\t \t omW:"<<seResult.omegaWorld.transpose()<<endl;}
  pBody_RPY_des[0] = rpydesWBC(0);//_roll_des;
  pBody_RPY_des[1] = rpydesWBC(1); 
  pBody_RPY_des[2] = rpydesWBC(2);//_yaw_des;//_yaw_des;
  /*vBody_Ori_des[0] = _roll_turn_rate;//;
  vBody_Ori_des[1] = 0.;
  vBody_Ori_des[2] = _yaw_turn_rate;*/
  vBody_Ori_des[0] = drpydesWBC(0);//_roll_turn_rate;//;
  vBody_Ori_des[1] = drpydesWBC(1);//0.;
  vBody_Ori_des[2] = drpydesWBC(2);//_yaw_turn_rate;*/
  // END of WBC Update
  iterationCounter++;

#ifdef DRAW_DEBUG_PATH
  auto* wpdWBCSphere = data.visualizationData->addSphere();
  if(wpdWBCSphere) {
    wpdWBCSphere->position << pBody_des;
    wpdWBCSphere->radius = 0.02;wpdWBCSphere->color = {0., 1.0, 0., 0.7};
  }
#endif

}
template<>
void ConvexMPCLocomotion::run(ControlFSMData<double>& data) {
  (void)data;
  ////printf("call to old CMPC with double!\n");
}
void ConvexMPCLocomotion::updateMPCIfNeeded(Gait* gait, int *mpcTable, ControlFSMData<float> &data, int horLength) {
  ////cout<<"HI16";
  //iterationsBetweenMPC = 30;
  if((iterationCounter % iterationsBetweenMPC) == 0){
    Timer MPCPrepTimer;
    if(gait==jumping){
      int summpcT=0;
      for (int indmp=0; indmp < 4; indmp++) summpcT+=mpcTable[indmp];//4*horLength
      if(summpcT!=0){
        flyingjump=false;
      }else if(!flyingjump){
        flyingjump=true; //CAREFUL WHEN REDIFINING JUMPS/LANDS
        jumping->countflystage+=1;cout<<"flyingjump true"<<endl;
        if(data.userParameters->noncopjump>0.5)apla=apla2;
    }}
    ////cout<<"HI18";
    drpydesWBC.setZero();
    auto seResult = data._stateEstimator->getResult();
    Vec3<float> orienvecdes(_roll_des,_pitch_des,_yaw_des);// + dtMPC * horLength * _yaw_turn_rate;/*if (currently_jumping)_pitch_des = 0;//.1415*2/5;else_pitch_des=0;*/
    Vec4<float> quatdes,quatdesWBC;
    if((current_gait == 4)||firstRun){//CONSIDER: && bswitchgait==false//// ///cout<<"StandMPC";
      if(data.userParameters->ori_pla_body<0.5){quatdes=rpyToQuat(orienvecdes);}
      else{quatdes=rotationMatrixToQuaternion(rpyToRotMat(orienvecdes)*Rpla[0].transpose());}
      if (0>seResult.orientation.transpose()*quatdes){quatdes=-quatdes;}
      float trajInitial[13] = {quatdes(0), quatdes(1), quatdes(2), quatdes(3),
        (float)world_position_desiredWBC[0],(float)world_position_desiredWBC[1],world_position_desiredWBC[2],//(float)_body_height,//world_position_desiredWBC[2];
        0,0,0,0,0,0};
      for(int i = 0; i < horLength; i++)
        for(int j = 0; j < 13; j++)
          trajAll[13*i+j] = trajInitial[j];
      quatdesWBC=quatdes;
    }else if(bswitchgait || gait!=jumping){//&& gait==&jumping1//TAKE OUT AND CONSIDER && gait==jumping (looks smoother to track ei previous traj for nonleap)//condition may glitch
    //////printf("Inside first switch\n");
      if(gait==jumping){
        cout<<"JumpStart";
        /////Timer T1Trajgen;
        float xFinal = 0;
        float yFinal = 0;
        float zFinal = 0;
        if(isbckflstab<0.5f){
          RecompHeurFl(data,quatdes,horLength,(int)0,(int)0,jumping->getFliMacNH(1),xFinal,yFinal,zFinal);//getFliMacNH//horLength is NOT, CORRECT(NP3-NP2)?
        }else{// IMPLEMENT MORPHED OFFLINE TRAJ HERE
          Vec4<float> quatprev(seResult.orientation);
          Eigen::Quaternionf quatprevq;quatprevq.w()=quatprev(0);quatprevq.vec()=quatprev.block(1,0,3,1);
          //if (0>quatprev(0)){quatprev=-quatprev;}//cout<<"quatdesmid"<<quatdesmid<<endl;
          float xStart = world_position_desiredWBC[0];
          float yStart = world_position_desiredWBC[1];
          float zStart = seResult.position[2];
          int NP1=0;
          int NP2=Np3+(Np2+1)/2;
          float tLO=NP2*dtMPC;//WELL YES BUT ACTUALLY NO
          

          if((float)data.userParameters->ifleap>0.5){//testlater
            v_des_worldWBC = coordinateRotation(CoordinateAxis::Z, -seResult.rpy(2)) * Vec3<float>((float)data.userParameters->dxjump, (float)data.userParameters->dyjump, 0);
          }
          if((float)data.userParameters->ifleap>0.5){//testlater
            v_des_worldWBC = coordinateRotation(CoordinateAxis::Z, -seResult.rpy(2)) * Vec3<float>((float)data.userParameters->dxjump, (float)data.userParameters->dyjump, 0);
          }
          //change lol
          float xMid = xStart;
          float yMid = yStart;
          float zMid = zStart;
          if(issideflstab>0.5f){
            xMid += 0.14 * v_des_worldWBC[0]/sqrt(v_des_worldWBC[0]*v_des_worldWBC[0]+v_des_worldWBC[1]*v_des_worldWBC[1]);//0.32rfeet-something
            yMid += 0.14 * v_des_worldWBC[1]/sqrt(v_des_worldWBC[0]*v_des_worldWBC[0]+v_des_worldWBC[1]*v_des_worldWBC[1]);//0.32
            zMid += 0.16;//0.1+0.32
          }else{
            xMid += 0.36 * v_des_worldWBC[0]/sqrt(v_des_worldWBC[0]*v_des_worldWBC[0]+v_des_worldWBC[1]*v_des_worldWBC[1]);//0.32rfeet-something
            yMid += 0.36 * v_des_worldWBC[1]/sqrt(v_des_worldWBC[0]*v_des_worldWBC[0]+v_des_worldWBC[1]*v_des_worldWBC[1]);//0.32
            zMid += 0.32;//0.1+0.32
          }
          
          Vec3<float> vH((xMid-xStart)/(tLO-NP1*dtMPC),(yMid-yStart)/(tLO-NP1*dtMPC),(zMid-zStart)/(tLO-NP1*dtMPC));//No gravity

          Vec4<float> quatdesmid;
          Vec3<float> orienvecdesmid(_roll_des,_pitch_des,_yaw_des);//7*M_PI/16.f
          if(issideflstab>0.5f){orienvecdesmid(0)-=M_PI_2;}
          else{orienvecdesmid(1)-=M_PI_2;}
          

          /*Vec4<float> qt4heuomid;
          qt4heuomid<<cos((-M_PI_2)/2),Vec3<float>(1,0,0)*sin((-M_PI_2)/2);
          Eigen::Quaternionf qt4heuoqmid;qt4heuoqmid.w()=qt4heuomid(0);qt4heuoqmid.vec()=qt4heuomid.block(1,0,3,1);
          qt4heuoqmid=quatprevq*qt4heuoqmid;     
          quatdesmid<<qt4heuoqmid.w(),qt4heuoqmid.vec();*/
          float shifanrot=0;
          if(data.userParameters->ori_pla_body<0.5){quatdesmid=rpyToQuat(orienvecdesmid);}
          else{quatdesmid=rotationMatrixToQuaternion(rpyToRotMat(orienvecdesmid)*(*Rplatarget)[0].transpose());}
          if (0>quatprev.transpose()*quatdesmid ){quatdesmid=-quatdesmid;shifanrot=0;}//2*M_PI


          /////printf("    NP1: %i, NP2: %i, NP3: %i\n", NP1, NP2, NP3);
          Vec4<float> quati;quati<< quatdesmid+quatprev;quati.normalize();//=;//cout<<"_pitch_des: "<<_pitch_des<<endl;//CAREFULLLL
          Eigen::Quaternionf quti;quti.w()=quati(0);quti.vec()=quati.block(1,0,3,1);
          Eigen::Matrix<float,3,3> Rqu = quti.toRotationMatrix();//VERIFY THIS CRAP
          Vec3<float> LHeur=4*Rqu*I_bod*Rqu.transpose()*QMATquatv(quati).transpose()*(quatdesmid-quatprev)/(tLO-2*dtMPC);//CORRECTTTTTTT NEGATIVE??? VERIFY UNDER YAW ROTS
          /*if(current_gait==6 && isjleap>0.5 && isnoncopj>0.5){LHeur+=Vec3<float>(0.18f,0,0);}//CORRECTTTTTTT FOR YAW ROTS nd frontal jumps*/
          if(issideflstab>0.5f){
            LHeur+=coordinateRotation(CoordinateAxis::Z, -seResult.rpy(2))*Vec3<float>(-0.5f,0,0);//NOTE SIGN
          }
          for(int i = NP1; i < NP2; i++){//(int)data.userParameters->Npar3
            quati<<quatdesmid*(i-NP1)*dtMPC+quatprev*(tLO-i*dtMPC);//cos(-M_PI*(i*dtMPC-tLO)/(tla-tLO)),Vec3<float>(1,0,0)*sin(-M_PI*(i*dtMPC-tLO)/(tla-tLO));
            quati<<(quatdesmid-quatprev)*(i-2)*dtMPC*(i>1)+quatprev*tLO-NP1*dtMPC*quatdesmid;//cos(-M_PI*(i*dtMPC-tLO)/(tla-tLO)),Vec3<float>(1,0,0)*sin(-M_PI*(i*dtMPC-tLO)/(tla-tLO));
            quati<<((quatdesmid-quatprev)*i*dtMPC-2*dtMPC*quatdesmid)*(i>1)+quatprev*tLO;//cos(-M_PI*(i*dtMPC-tLO)/(tla-tLO)),Vec3<float>(1,0,0)*sin(-M_PI*(i*dtMPC-tLO)/(tla-tLO));
            quati.normalize();//cout<<i<<", quati: "<<quati.transpose()<<endl;
            trajHeurJump[13*i + 0] = quati(0);
            trajHeurJump[13*i + 1] = quati(1);
            trajHeurJump[13*i + 2] = quati(2);
            trajHeurJump[13*i + 3] = quati(3);
            trajHeurJump[13*i + 4] = xStart + vH(0)*(i-NP1)*dtMPC;
            trajHeurJump[13*i + 5] = yStart + vH(1)*(i-NP1)*dtMPC;
            trajHeurJump[13*i + 6] = zStart + vH(2)*(i-NP1)*dtMPC;// - (9.81f/2)*(i*dtMPC-tLO)*(i*dtMPC-tLO)//dtMPC * (apla2(1)*v_des_world[0]+apla2(2)*v_des_world[1])
            for(int j = 7; j < 10; j++){
              trajHeurJump[13*i + j] = (i>1)*LHeur(j-7);}
            trajHeurJump[13*i + 10] = vH(0);
            trajHeurJump[13*i + 11] = vH(1);
            trajHeurJump[13*i + 12] = vH(2);//dtMPC * (apla2(1)*v_des_world[0]+apla2(2)*v_des_world[1])
          }
          int NP3=jumping->getFliMacNH(1);
          float tla=NP3*dtMPC;
          xFinal = xMid+(NP3-NP2)*dtMPC * v_des_worldWBC[0];
          yFinal = yMid+(NP3-NP2)*dtMPC * v_des_worldWBC[1];
          zFinal=(*aplatarget)(0)+(*aplatarget)(1)*xFinal+(*aplatarget)(2)*yFinal+(float)data.userParameters->bo_height+0.18;// CHANGE
          
          /*_yaw_des=seResult.rpy(2)+yawrotjump;//yawrotjump deactivated for now modify
          Vec3<float> orienvecdes1(_roll_des,_pitch_des,_yaw_des);
          if(data.userParameters->ori_pla_body<0.5){quatdes=rpyToQuat(orienvecdes1);}
          else{quatdes=rotationMatrixToQuaternion(rpyToRotMat(orienvecdes1)*(*Rplatarget)[0].transpose());}
          if (0>quatprev.transpose()*quatdes ){quatdes=-quatdes;}//cout<<"quatdes"<<quatdes<<endl;*/
          Vec4<float> qt4heuo;
          Vec3<float>axrot(0,0,0);
          if(issideflstab>0.5f){axrot(0)=1;}
          else{axrot(1)=1;}
          qt4heuo<<cos((shifanrot-M_PI_2-(2*M_PI-M_PI_2)/2)/2),axrot*sin((shifanrot-M_PI_2-(2*M_PI-M_PI_2)/2)/2);
          //if(quatprev(0)<0){qt4heuo=-qt4heuo;}
          //Eigen::Quaternionf quatprevq(quatprev);
          Eigen::Quaternionf qt4heuoq;qt4heuoq.w()=qt4heuo(0);qt4heuoq.vec()=qt4heuo.block(1,0,3,1);
          qt4heuoq=quatprevq*qt4heuoq;
          qt4heuo<<qt4heuoq.w(),qt4heuoq.vec();
          quati<< quatdesmid+qt4heuo;quati.normalize();//=quatdes+quatprev;//cout<<"_pitch_des: "<<_pitch_des<<endl;//CAREFULLLL
          quti.w()=quati(0);quti.vec()=quati.block(1,0,3,1);//quti=Eigen::Quaternionf(quati)
          Rqu = quti.toRotationMatrix();
          LHeur=2*4*Rqu*I_bod*Rqu.transpose()*QMATquatv(quati).transpose()*(qt4heuo-quatdesmid)/(tla-tLO);//4*I_wor*QMATquatv(quati).transpose()*(quatdes-quatprev)/(tla-tLO);//CORRECTTTTTTT NEGATIVE??? VERIFY UNDER YAW ROTS
          //verify if x2 is ok or we neeed x4
          /*if(current_gait==6 && isjleap>0.5 && isnoncopj>0.5){LHeur+=Vec3<float>(0.18f,0,0);}//CORRECTTTTTTT FOR YAW ROTS nd frontal jumps*/
          vH<<v_des_worldWBC[0],v_des_worldWBC[1],(zFinal-zMid)/(tla-tLO)+9.81f*(tla-tLO)/2;
          for(int i = NP2; i < NP3; i++){//(int)data.userParameters->Npar3
            quati<<cos((shifanrot-M_PI_2+0.001-(2*M_PI-M_PI_2)*(i*dtMPC-tLO)/(tla-tLO))/2),axrot*sin((shifanrot-M_PI_2+0.001-(2*M_PI-M_PI_2)*(i*dtMPC-tLO)/(tla-tLO))/2);//quatdes*(i*dtMPC-tLO)+quatprev*(tla-i*dtMPC);
            //if(quatprev(0)<0){quati=-quati;}
            Eigen::Quaternionf quatiq;quatiq.w()=quati(0);quatiq.vec()=quati.block(1,0,3,1);
            quatiq=quatprevq*quatiq;
            quati<<quatiq.w(),quatiq.vec();
            quati.normalize();//cout<<i<<", quati: "<<quati.transpose()<<endl;
            trajHeurJump[13*i + 0] = quati(0);
            trajHeurJump[13*i + 1] = quati(1);
            trajHeurJump[13*i + 2] = quati(2);
            trajHeurJump[13*i + 3] = quati(3);
            trajHeurJump[13*i + 4] = xMid + vH(0)*(i*dtMPC-tLO);
            trajHeurJump[13*i + 5] = yMid + vH(1)*(i*dtMPC-tLO);
            trajHeurJump[13*i + 6] = zMid + vH(2)*(i*dtMPC-tLO) - (9.81f/2)*(i*dtMPC-tLO)*(i*dtMPC-tLO);//dtMPC * (apla2(1)*v_des_world[0]+apla2(2)*v_des_world[1])
            for(int j = 7; j < 10; j++){
              trajHeurJump[13*i + j] = LHeur(j-7);}
            trajHeurJump[13*i + 10] = vH(0);
            trajHeurJump[13*i + 11] = vH(1);
            trajHeurJump[13*i + 12] = vH(2) - (i*dtMPC-tLO) * 9.81f;//dtMPC * (apla2(1)*v_des_world[0]+apla2(2)*v_des_world[1])
          }
        }
        /////printf("   T1Trajgen: %.3f\n", T1Trajgen.getMs());
        /////Timer T2Trajgen;
        if(jumping->ischained){
          Vec4<float> quatdes2;
          Vec3<float> orienvecdes2(_roll_des,_pitch_des,_yaw_des+yawrotjump);// + dtMPC * horLength * _yaw_turn_rate;
          if(data.userParameters->ori_pla_body<0.5){quatdes2=rpyToQuat(orienvecdes2);}
          else{quatdes2=rotationMatrixToQuaternion(rpyToRotMat(orienvecdes2)*(*Rplatarget)[0].transpose());}
          if (0>quatdes.transpose()*quatdes2){quatdes2=-quatdes2;}
          Vec3<float> v_des_worldHM2(coordinateRotation(CoordinateAxis::Z, -(float)yawrotjump)*v_des_worldWBC);
          float xFinal2 = xFinal+horLength*dtMPC*v_des_worldHM2[0];
          float yFinal2 = yFinal+horLength*dtMPC*v_des_worldHM2[1];
          float heightoffCoM_f2=(*aplatarget)(0)+(*aplatarget)(1)*xFinal2+(*aplatarget)(2)*yFinal2;//_body_height = heightoffCoM + (float)data.userParameters->bo_hei_jum;//JUMPINGUNCOM
          float zFinal2=heightoffCoM_f2+(float)data.userParameters->bo_height;
          float tLO2=(jumping->getFliMacNH(1)+Np2)*dtMPC;//Np2???? CORRECT !!!!!!!!!!!!!!!!!!
          float tla2=jumping->getFliMacNH(2)*dtMPC;
          /////printf("   T2Trajgen: %.3f\n", T2Trajgen.getMs());
          /////Timer T3Trajgen;
          HeurFlighComp(xFinal, yFinal, zFinal, xFinal2, yFinal2, zFinal2, quatdes, quatdes2, tLO2, tla2, jumping->getFliMacNH(1), jumping->getFliMacNH(1)+Np2, jumping->getFliMacNH(2)); // Np2 CORRECT 
        }
        /////printf("   T3Trajgen: %.3f\n", T3Trajgen.getMs());
        /////Timer T4Trajgen;
        for (int i = 0; i < (13*jumping->getFliMacNH(2)); i++){trajAll[i]=trajHeurJump[i];}//horLength only maybe?
        /////printf("   T4Trajgen: %.3f\n", T4Trajgen.getMs());
        /////Timer T5Trajgen;
        print_named_array("Heuristic trajectory: \n",trajHeurJump,jumping->getFliMacNH(2),13);
        /////printf("   T5Trajgen: %.3f\n", T5Trajgen.getMs());
      }else{
        if((float)data.userParameters->ifleap>0.5 && bswitchgait){//corrected VERIFYYYYYYY
          world_position_desiredWBC[0] = (pFoot_des[0](0)+pFoot_des[1](0)+pFoot_des[2](0)+pFoot_des[3](0))/4;
          world_position_desiredWBC[1] = (pFoot_des[0](1)+pFoot_des[1](1)+pFoot_des[2](1)+pFoot_des[3](1))/4;
        }
        ///////cout<<"JumpEnd"<<endl;
        if(data.userParameters->ori_pla_body<0.5){quatdes=rpyToQuat(orienvecdes);}
        else{quatdes=rotationMatrixToQuaternion(rpyToRotMat(orienvecdes)*Rpla[0].transpose());}
        if (0>seResult.orientation.transpose()*quatdes){quatdes=-quatdes;}
        //cout << "v_des_world: "<<v_des_world<< endl;
        float xStart = world_position_desiredWBC[0];
        float yStart = world_position_desiredWBC[1];
        float zStart = seResult.position[2];
        float xFinal = xStart+horLength*dtMPC * v_des_worldWBC[0];
        float yFinal = yStart+horLength*dtMPC * v_des_worldWBC[1];
        float zFinal = apla(0)+apla(1)*xFinal+apla(2)*yFinal+(float)data.userParameters->bo_height; // bo1->bo2->jump
        float _body_height_fdot=(zFinal-zStart)/horLength;
        v_des_worldWBC[2]=_body_height_fdot/dtMPC;
        float trajInitial[13] = {quatdes(0),quatdes(1),quatdes(2),quatdes(3),
        xStart,yStart,zStart,0,0,0,v_des_worldWBC[0],v_des_worldWBC[1],v_des_worldWBC[2]};
        for(int i = 0; i < horLength; i++){//cout<<i<<",";
          for(int j = 0; j < 13; j++)
            trajAll[13*i+j] = trajInitial[j];
          if(i != 0){
            trajAll[13*i + 4] = trajAll[13 * (i - 1) + 4] + dtMPC * v_des_worldWBC[0];
            trajAll[13*i + 5] = trajAll[13 * (i - 1) + 5] + dtMPC * v_des_worldWBC[1];
            trajAll[13*i + 6] = trajAll[13 * (i - 1) + 6] + _body_height_fdot;//dtMPC * (apla2(1)*v_des_world[0]+apla2(2)*v_des_world[1])
          }
        }
      }
      quatdesWBC<<trajAll[13],trajAll[14],trajAll[15],trajAll[16];//VERIFY WHAT'S FIRST COMP IN MPC
      bswitchgait = false;//may glitch with condition FOR NOT FORGETTING PREVIOUS TRAJ.,.I.E. CHAINING TRAJS.
    }else{
     //// ///cout<<"TrackPrevTrajMPC";//cout<<"TrackPrevTrajMPC,swfromjump:"<<swfromjump;
      //update_traj(trajAll);//cout<<"trajAll: "<< *trajAll<<endl;//////printf("Tracking previous MPC traj\n");
      if(gait==jumping){
        //cout<< "horLength"<<horLength<<endl;
        /*if (horLength>(-1+(int)data.userParameters->Npar1-(int)data.userParameters->Npar3)){
          for (int i = 13*0; i < 13*horLength; i++){//13*4
            trajAll[i]=(0.0*trajAll[i]+1.0*trajHeurJump[13*((int)data.userParameters->Npar1-horLength)+i]);
        }}*/
        if(jumping->ischained){
          if(jumping->getCurrentGaitPhase()==jumping->getFliMacNH(1)){//Careful with state-based phase
            RecompHeurFl(data,quatdes,horLength,jumping->getFliMacNH(1),jumping->getFliMacNH(1)+Np2,jumping->getFliMacNH(2));// Np2 CORRECT !!!
            orienvecdes<<_roll_des,_pitch_des,_yaw_des;
            //print_named_array("Heuristic trajectory: \n",trajHeurJump,jumping->getFliMacNH(2),13);
          }
        }
        
        for (int i = 13*0; i < 13*fmin(horLength,jumping->getFliMacNH(2)-jumping->getCurrentGaitPhase()); i++){//13*4 13*(int)34  horLength<(int)34-jumping->getCurrentGaitPhase()
          trajAll[i]=(0.0*trajAll[i]+1.0*trajHeurJump[13*(jumping->getCurrentGaitPhase())+i]);//(int)data.userParameters->Npar1-horLength
        }

        //cout<<endl;cout<<jumping->getCurrentGaitPhase()<<endl;
        if((float)data.userParameters->ifleap>0.5){//WELL YES BUT ACTUALLY NO
          v_des_worldWBC = coordinateRotation(CoordinateAxis::Z, -seResult.rpy(2)) * Vec3<float>((float)data.userParameters->dxjump, (float)data.userParameters->dyjump, 0);
        }
        v_des_worldWBC<<trajAll[10],trajAll[11],trajAll[12];//CHANGE TO HEUR?//OVERTRUSTING HEURISTIC?//WRONG TRAJALL IS SHIFTED 
        /////cout<<endl;cout<<"v_des_worldWBC: "<<v_des_worldWBC<<endl;
        quatdesWBC<<trajAll[13],trajAll[14],trajAll[15],trajAll[16];//CHANGE TO HEUR?//OVERTRUSTING HEURISTIC?
        //for(int ic=0; ic<17;ic++) cout<<"trajAll: "<< trajAll[ic]<<endl;//////printf("Tracking previous MPC traj\n");
        //quatdesWBC<<trajHeurJump[13],trajHeurJump[14],trajHeurJump[15],trajHeurJump[16];//CHANGE TO HEUR?//OVERTRUSTING HEURISTIC?
        Eigen::Quaternionf quti;quti.w()=seResult.orientation(0);quti.vec()=seResult.orientation.block(1,0,3,1);//CAREFULLLL
        Eigen::Matrix<float,3,3> Rqu = quti.toRotationMatrix();
        Eigen::Matrix<float,3,3> I_worinv,EomtodrpyMM;
        I_worinv = Rqu*(I_bod.inverse())*Rqu.transpose();
        //Vec3<float> omWBCHeur=I_worinv*Vec3<float>(trajHeurJump[7],trajHeurJump[8],trajHeurJump[9]);//WELL YES BUT ACTUALLY NO //add CMM
        Vec3<float> omWBCHeur=I_worinv*Vec3<float>(trajAll[7],trajAll[8],trajAll[9])*0.3;//WELL YES BUT ACTUALLY NO //add CMM
        if(jumping->getCurrentGaitPhase()>Np3 && issideflstab>0.5f && isbckflstab>0.5f){
          omWBCHeur=omWBCHeur*2;
        }
        EomtodrpyMM=Eomtodrpy(seResult.rpy(1), seResult.rpy(2));
        drpydesWBC=EomtodrpyMM*omWBCHeur;

        heightoffCoM=(*aplatarget)(0)+(*aplatarget)(1)*(world_position_desiredWBC[0] + dtMPC * horLength * v_des_worldWBC[0])+(*aplatarget)(2)*(world_position_desiredWBC[1] + dtMPC * horLength * v_des_worldWBC[1]);//VERIFY
        _body_height = heightoffCoM + (float)data.userParameters->bo_height; // JUMPINGUNCOM
        
        if(data.userParameters->ori_pla_body<0.5){quatdes=rpyToQuat(orienvecdes);}
        else{quatdes=rotationMatrixToQuaternion(rpyToRotMat(orienvecdes)*(*Rplatarget)[0].transpose());}
      }else{if(data.userParameters->ori_pla_body<0.5){quatdes=rpyToQuat(orienvecdes);}
        else{quatdes=rotationMatrixToQuaternion(rpyToRotMat(orienvecdes)*Rpla[0].transpose());}}

      if(gait!=jumping || (horLength>(jumping->getFliMacNH(2)-jumping->getCurrentGaitPhase()))){  //  todo  : change this
        Vec4<float> quatprev(trajAll[13*(horLength-1)+0],trajAll[13*(horLength-1)+1],trajAll[13*(horLength-1)+2],trajAll[13*(horLength-1)+3]);
        if (0>quatprev.transpose()*quatdes){quatdes=-quatdes;}
        float trajEnd[13] = {quatdes(0),quatdes(1),quatdes(2),quatdes(3),
          world_position_desiredWBC[0] + dtMPC * horLength * v_des_worldWBC[0],
          world_position_desiredWBC[1] + dtMPC * horLength * v_des_worldWBC[1],
          (float)_body_height + dtMPC * horLength * ((*aplatarget)(1)*v_des_worldWBC[0]+(*aplatarget)(2)*v_des_worldWBC[1])+0.18f,0,0,0,
          v_des_worldWBC[0],v_des_worldWBC[1],(*aplatarget)(1)*v_des_worldWBC[0]+(*aplatarget)(2)*v_des_worldWBC[1]};
        /*for(s16 itra = 0; itra < 13; itra++)
          trajAll[13*(horLength-1)+itra] = trajEnd[itra];//becuase kjoyst gives a roughr*/

        for(int i = horLength-1; i > max(jumping->getFliMacNH(2)-jumping->getCurrentGaitPhase()-1,-1); i--){//cout<<i<<",";
          for(s16 j = 0; j < 13; j++)
            trajAll[13*i+j] = trajEnd[j];
          if(i != horLength-1){
            trajAll[13*i + 4] = trajAll[13 * (i + 1) + 4] - dtMPC * v_des_worldWBC[0];
            trajAll[13*i + 5] = trajAll[13 * (i + 1) + 5] - dtMPC * v_des_worldWBC[1];
            trajAll[13*i + 6] = trajAll[13 * (i + 1) + 6] - dtMPC * ((*aplatarget)(1)*v_des_worldWBC[0]+(*aplatarget)(2)*v_des_worldWBC[1]);//dtMPC * (apla2(1)*v_des_world[0]+apla2(2)*v_des_world[1])
          }
        }
      }
    }
    if(gait!=jumping){
      quatdesWBC=quatdes;
    }
    rpydesWBC=quatToRPY(quatdesWBC);//PREVIOUSLY:quatdesWBC=quatdes////MAY NOT WORK IN NONCOPLANAR//CHECCKKKKKKKKKKK//SOLVE, THIS IS THE ERROR
    printf("  Prep MPC Time: %.3f\n", MPCPrepTimer.getMs());
    Timer solveTimerSDMPC;
    solveDenseMPC(mpcTable, data, horLength);
    printf("  SOLVE TIME DenseMPC: %.3f\n", solveTimerSDMPC.getMs());
    /*cout<<"Fr_des1: \n"<<Fr_des[0]<<endl;
   //// ///cout<<"Fr_des2: \n"<<Fr_des[1]<<endl;
   //// ///cout<<"Fr_des3: \n"<<Fr_des[2]<<endl;
   //// ///cout<<"Fr_des4: \n"<<Fr_des[3]<<endl;*/
  }
}
void ConvexMPCLocomotion::solveDenseMPC(int *mpcTable, ControlFSMData<float> &data, int horLength) {
  auto seResult = data._stateEstimator->getResult();
  Matrix<double,12,1> QMatr;
  QMatr.block(0,0,3,1)=data.userParameters->Qwei2;
  QMatr.block(3,0,3,1)=data.userParameters->Qwei3;
  QMatr.block(6,0,3,1)=data.userParameters->Qwei4;
  QMatr.block(9,0,3,1)=data.userParameters->Qwei5;
  double* weights = QMatr.data();
  double* uweights = data.userParameters->Qweiu.data();
  float* p = seResult.position.data();
  float* v = seResult.vWorld.data();
  float* w = seResult.omegaWorld.data();
  /*Vec4<float> qpreor(seResult.orientation);
  if (0>seResult.orientation(0)){qpreor=-qpreor;}//cout<<"quatdesmid"<<quatdesmid<<endl;
  float* q = qpreor.data();*/

  float* q = seResult.orientation.data();
  float seLcomp[3]={0.,0.,0.};
  
  if((current_gait==6 && isjleap>0.5 && isnoncopj>0.5)||(isbckflstab>0.5f && issideflstab>0.5f && current_gait==6)){
    Vec3<float> seLcompv = coordinateRotation(CoordinateAxis::Z, -seResult.rpy(2))*Vec3<float>(-0.5f,0.f,0.f);
    seLcomp[0]=seLcompv(0);seLcomp[1]=seLcompv(1);seLcomp[2]=0.;//sign changed
  }else{
  }
  float r[12];
  for(int i = 0; i < 12; i++) r[i] = pFoot[i%4][i/4]-seResult.position[i/4];
  int summpcT=0;
  for (int indmp=0; indmp < 4*horLength; indmp++) summpcT+=mpcTable[indmp];
  dtMPC = dt * iterationsBetweenMPC;
  muv[0]=(float)data.userParameters->muest;muv[1]=muv[0];muv[2]=muv[0];muv[3]=muv[0];
  ////printf("presetup_problem");is_slidcntct
  Timer t1;
  setup_problem(dtMPC,horLength,muv,120,summpcT, Rpla, data.userParameters->savdat);
  ////printf("postsetup_problem");
  printf("    setup_problem: %.3f\n", t1.getMs());
  Timer t2;
  update_problem_data_floats(p,v,q,w,r,weights,trajAll,uweights,mpcTable,seLcomp);
  printf("    update_problem_data_floats: %.3f\n", t2.getMs());
  for(int leg = 0; leg < 4; leg++){
    Vec3<float> f;
    for(int axis = 0; axis < 3; axis++)
      f[axis] = get_solution(leg*3 + axis);
    printf("[%d] %7.3f %7.3f %7.3f\n", leg, f[0], f[1], f[2]);
    f_ff[leg] = -seResult.rBody * f;
    Fr_des[leg] = f;// Update for WBC
}}

void ConvexMPCLocomotion::HeurFlighComp(float xStart,float yStart,float zStart, float xFinal, float yFinal, float zFinal,const Vec4<float>& quatprev,const Vec4<float>& quatdes, float tLO, float tla, int hhL1, int hhL2, int hhL3){//implement offline-obtained trajs here
  /////printf("    hhL1: %i, hhL2: %i, hhL3: %i\n", hhL1, hhL2, hhL3);
  Vec3<float> vH((xFinal-xStart)/(tla-tLO),(yFinal-yStart)/(tla-tLO),(zFinal-zStart)/(tla-tLO)+9.81f*(tla-tLO)/2);
  float traj0tLO[13] = {quatprev(0),quatprev(1),quatprev(2),quatprev(3),
  xStart,yStart,zStart,0,0,0,0,0,0};//consider initial velocity MODIFYYYYYY
  for(int i = hhL1; i < hhL2; i++){    //cout<<i<<",";
    for(int j = 0; j < 13; j++)
      trajHeurJump[13*i+j] = traj0tLO[j];//CAREFUL WITH TVLQR
  }
  Vec4<float> quati=quatdes+quatprev;
  if(isbckflstab>0.5f){
    quati<< 1, 0, -1, 0;//cout<<"_pitch_des: "<<_pitch_des<<endl;
  }
  quati.normalize(); //CAREFULLLL
  Eigen::Quaternionf quti;quti.w()=quati(0);quti.vec()=quati.block(1,0,3,1);//q.w() = q_[0];//this->q.x() = q_[1];//this->q.y() = q_[2];//this->q.z() = q_[3];
  Eigen::Matrix<float,3,3> Rqu = quti.toRotationMatrix();
  Eigen::Matrix<float,3,3> I_wor;
  I_wor = Rqu*I_bod*Rqu.transpose();
  /*if(isbckflstab>0.5f){
    _pitch_des=_pitch_des-7*M_PI/8;
    //cout<<"_pitch_des: "<<_pitch_des<<endl;
  }*/
  Vec3<float> LHeur=4*I_wor*QMATquatv(quati).transpose()*(quatdes-quatprev)/(tla-tLO);//CORRECTTTTTTT NEGATIVE??? VERIFY UNDER YAW ROTS
  if(isbckflstab>0.5f){
    LHeur=2*4*I_wor*QMATquatv(quati).transpose()*(Vec4<float>(0,0,-1,0)-quatprev)/(tla-tLO);//CORRECTTTTTTT FOR YAW ROTS
  }
  if(current_gait==6 && isjleap>0.5 && isnoncopj>0.5){//CORRECTTTTTTT FOR YAW ROTS nd frontal jumps
    LHeur+=Vec3<float>(0.18f,0,0);
  }
  quati.setZero();
  for(int i = hhL2; i < hhL3; i++){//(int)data.userParameters->Npar3
    if(isbckflstab<0.5f){
      quati=quatdes*(i*dtMPC-tLO)+quatprev*(tla-i*dtMPC);
      quati.normalize();
    }else{
      quati<<cos(-M_PI*(i*dtMPC-tLO)/(tla-tLO)),Vec3<float>(1,0,0)*sin(-M_PI*(i*dtMPC-tLO)/(tla-tLO));//quatdes*(i*dtMPC-tLO)+quatprev*(tla-i*dtMPC);
      quati.normalize();
    }
    //cout<<i<<", quati: "<<quati.transpose()<<endl;
    trajHeurJump[13*i + 0] = quati(0);
    trajHeurJump[13*i + 1] = quati(1);
    trajHeurJump[13*i + 2] = quati(2);
    trajHeurJump[13*i + 3] = quati(3);
    trajHeurJump[13*i + 4] = xStart + vH(0)*(i*dtMPC-tLO);
    trajHeurJump[13*i + 5] = yStart + vH(1)*(i*dtMPC-tLO);
    trajHeurJump[13*i + 6] = zStart + vH(2)*(i*dtMPC-tLO) - (9.81f/2)*(i*dtMPC-tLO)*(i*dtMPC-tLO);//dtMPC * (apla2(1)*v_des_world[0]+apla2(2)*v_des_world[1])
    for(int j = 7; j < 10; j++){
      trajHeurJump[13*i + j] = LHeur(j-7);}
    trajHeurJump[13*i + 10] = vH(0);
    trajHeurJump[13*i + 11] = vH(1);
    trajHeurJump[13*i + 12] = vH(2) - (i*dtMPC-tLO) * 9.81f;//dtMPC * (apla2(1)*v_des_world[0]+apla2(2)*v_des_world[1])
  }
}