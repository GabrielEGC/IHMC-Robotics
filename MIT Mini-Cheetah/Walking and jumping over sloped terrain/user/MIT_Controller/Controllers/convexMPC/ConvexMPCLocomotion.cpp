#include <iostream>
#include <Utilities/Timer.h>
#include <Utilities/Utilities_print.h>
#include <Math/orientation_tools.h>
#include <eigen3/Eigen/Dense>
#include "ConvexMPCLocomotion.h"
#include "convexMPC_interface.h"
#include "../../../../common/FootstepPlanner/GraphSearch.h"
#include "Gait.h"
#include <Utilities/pseudoInverse.h>
using Eigen::Matrix;
using Eigen::Array4i;
using Eigen::Dynamic;
using Eigen::Map;
using std::cout;
using std::endl;
//using Eigen::Quaternionf;
#define DRAW_DEBUG_SWINGS
#define DRAW_DEBUG_PATH
////////////////////
// Controller
////////////////////
ConvexMPCLocomotion::ConvexMPCLocomotion(float _dt, int _iterations_between_mpc, MIT_UserParameters* parameters) :
  iterationsBetweenMPC(_iterations_between_mpc),
  horizonLength(10),
  dt(_dt),
  trotting(horizonLength, Vec4<int>(0,5,5,0), Vec4<int>(5,5,5,5),"Trotting"),
  bounding(horizonLength, Vec4<int>(0,0,4,4),Vec4<int>(4,4,4,4),"Bounding"),
  //bounding(horizonLength, Vec4<int>(5,5,0,0),Vec4<int>(3,3,3,3),"Bounding"),
  pronking(horizonLength, Vec4<int>(0,0,0,0),Vec4<int>(4,4,4,4),"Pronking"),
  jumping((int)parameters->Npar1, Vec4<int>(0,0,0,0), Vec4<int>((int)parameters->Npar3,(int)parameters->Npar3,(int)parameters->Npar3,(int)parameters->Npar3), "Jumping"),
  //galloping(horizonLength, Vec4<int>(0,2,7,9),Vec4<int>(6,6,6,6),"Galloping"),
  //galloping(horizonLength, Vec4<int>(0,2,7,9),Vec4<int>(3,3,3,3),"Galloping"),
  galloping(horizonLength, Vec4<int>(0,2,7,9),Vec4<int>(4,4,4,4),"Galloping"),
  standing(horizonLength, Vec4<int>(0,0,0,0),Vec4<int>(10,10,10,10),"Standing"),
  //trotRunning(horizonLength, Vec4<int>(0,5,5,0),Vec4<int>(3,3,3,3),"Trot Running"),
  trotRunning(horizonLength, Vec4<int>(0,5,5,0),Vec4<int>(4,4,4,4),"Trot Running"),
  walking(horizonLength, Vec4<int>(0,3,5,8), Vec4<int>(5,5,5,5), "Walking"),
  walking2(horizonLength, Vec4<int>(0,5,5,0), Vec4<int>(7,7,7,7), "Walking2"),
  pacing(horizonLength, Vec4<int>(5,0,5,0),Vec4<int>(5,5,5,5),"Pacing"),
  oneaerial(12, Vec4<int>(0,3,6,9), Vec4<int>(9,9,9,9), "Oneaerial"),
  random(horizonLength, Vec4<int>(9,13,13,9), 0.4, "Flying nine thirteenths trot"),
  random2(horizonLength, Vec4<int>(8,16,16,8), 0.5, "Double Trot")
{
  _parameters = parameters;
  dtMPC = dt * iterationsBetweenMPC;
  default_iterations_between_mpc = iterationsBetweenMPC;
  printf("[Convex MPC] dt: %.3f iterations: %d, dtMPC: %.3f\n", dt, iterationsBetweenMPC, dtMPC);
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
   apla2(0)=-0.00357554;
   apla2(1)=0;
   apla2(2)=-0.269353;
   rpyplane.setZero();
   rpydesfin.setZero();

    Vec3<float> norpla(-apla(1),-apla(2),1);
    Vec3<float> nor2(apla(2),-apla(1),0);
    norpla.normalize();
    nor2.normalize();
    Vec3<float> nor3;
    nor3 = nor2.cross(norpla);
    Rpla2[0].block(0,0,3,1)=nor3;
    Rpla2[0].block(0,1,3,1)=nor2;
    Rpla2[0].block(0,2,3,1)=norpla;
    rpyplane=rotationMatrixToRPY(Rpla[0]);
    rpyplane(2)=0;
    Rpla2[0]=rpyToRotMat(rpyplane);
    Rpla2[1]=rpyToRotMat(rpyplane);
    Rpla2[2]=rpyToRotMat(rpyplane);
    Rpla2[3]=rpyToRotMat(rpyplane);

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
    y_vel_cmd = y_vel_cmd-fminf(fmaxf(y_vel_cmd,-0.3f),0.3f);
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
  //cout << _roll_des << endl;
  //cout << _yaw_des << endl;
  //_roll_des = 0.;
  _pitch_des = (float)(data.userParameters->pitchdesin);//*rpyplane(1);
  //_roll_turn_rate = 0.;
}

template<>
void ConvexMPCLocomotion::run(ControlFSMData<float>& data) {
  bool omniMode = false;
  // Command Setup
  _SetupCommand(data);
  gaitNumber = data.userParameters->cmpc_gait;
  if(gaitNumber >= 10) {
    gaitNumber -= 10;
    omniMode = true;
  }
  auto& seResult = data._stateEstimator->getResult();
  // Check if transition to standing
  if(((gaitNumber == 4) && current_gait != 4) || firstRun){
    stand_traj[0] = seResult.position[0];
    stand_traj[1] = seResult.position[1];
    stand_traj[2] = 0.21;
    stand_traj[3] = 0;
    stand_traj[4] = seResult.rpy[0];
    stand_traj[5] = seResult.rpy[2];
    world_position_desired[0] = stand_traj[0];
    world_position_desired[1] = stand_traj[1];
  }
  previous_gait = current_gait;
  Gait* gait = &trotting;//9
  Gait* nxtgait = &trotting;//9
  if(gaitNumber == 1)
    gait = &bounding;
  else if(gaitNumber == 2)
    gait = &pronking;
  else if(gaitNumber == 3)
    gait = &random;
  else if(gaitNumber == 4)
    gait = &standing;
  else if(gaitNumber == 5)
    gait = &oneaerial;//trotRunning;
  else if(gaitNumber == 6)
    gait = &jumping;
    //gait = &random2;
  else if(gaitNumber == 7)
    gait = &galloping;
  else if(gaitNumber == 8)
    gait = &pacing;

  nxtgait = gait;
  current_gait = gaitNumber;
  this->jumping._durations=Array4i((int)data.userParameters->Npar3,(int)data.userParameters->Npar3,(int)data.userParameters->Npar3,(int)data.userParameters->Npar3);
  this->jumping._nIterations=(int)data.userParameters->Npar1;
  gait->setIterations(iterationsBetweenMPC, iterationCounter);
  nxtgait->setIterations(iterationsBetweenMPC, iterationCounter);
  jumping.setIterations(iterationsBetweenMPC, iterationCounter);
  jumping.setIterations(27/2, iterationCounter);
  //printf("[%d] [%d]\n", jumping.get_current_gait_phase(), gait->get_current_gait_phase());
  jump_state.trigger_pressed(jump_state.should_jump(jumping.getCurrentGaitPhase(),(int)data.userParameters->Npar2),
      data._desiredStateCommand->trigger_pressed);//modif7
  // bool too_high = seResult.position[2] > 0.29;
  if(jump_state.should_jump(jumping.getCurrentGaitPhase(),(int)data.userParameters->Npar2)) {
    if(jumping.getCurrentGaitPhase() != (int)data.userParameters->Npar2 || jump_state.seen_end_count != (jump_state.END_COUNT-1)){
      gait = &jumping;
      current_gait = 6;
      recompute_timing(27/2);
      _body_height = heightoffCoM+(float)data.userParameters->bo_hei_jum;
      currently_jumping = true;
    }
  } else {
    recompute_timing(default_iterations_between_mpc);//switch beetwen iBMPC and 27/2 for jumping...
    currently_jumping = false;
  }
  if(_body_height < 0.02) {
    _body_height = 0.29;
  }
  if (current_gait!=previous_gait)//&& previous_gait!=6smooth transitions available
    bswitchgait = true;
  /*printf("Currently Jumping: %i\n", (int)currently_jumping);
  printf("Jump pending: %i\n", (int)jump_state.jump_pending);
  printf("Jump in progress: %i\n", (int)jump_state.jump_in_progress);
  printf("Current Gait: %i\n", (int)current_gait);*/
  int* mpcTable = gait->getMpcTable();
  jump_state.cshift=0;
  if (gait == &jumping){
    if (jump_state.seen_end_count==(jump_state.END_COUNT-1)){
      int gaitHor = gait->getHoriz();
      int gaitcurp = gait->getCurrentGaitPhase();
      int* nxtmpcTable = nxtgait->getMpcTable();
      int nxtgaithor = nxtgait->getHoriz();
      int nxtintpick = (-(nxtgait->getCurrentGaitPhase()) + 1) % nxtgaithor;// //++
      nxtintpick = (nxtintpick+nxtgaithor);
      for(int i = 0; i < gaitcurp+1; i++) {
        jump_state.cshift+=1;
        for(int j = 0; j < 4; j++) {
          mpcTable[4*(gaitHor-1-gaitcurp+i)+j]=nxtmpcTable[4*((nxtintpick+i)% nxtgaithor)+j];
  }}}}
  else if(previous_gait == 6){
      gait->gaitphascorr = (gait->gaitphascorr+nxtgait->getHoriz()-(nxtgait->getCurrentGaitPhase())+1) % nxtgait->getHoriz();// //++
      gait->setIterations(iterationsBetweenMPC, iterationCounter);
      mpcTable = gait->getMpcTable();
  }
  /*printf("Previous: %i\n", previous_gait);
  printf("Current:%i\n", gaitNumber);
  printf("Gait Phase Correction:%i\n", gait->gaitphascorr);*/
  // integrate position setpoint
  Vec3<float> v_des_robot(_x_vel_des, _y_vel_des, 0);
  Vec3<float> v_des_world = omniMode ? v_des_robot : seResult.rBody.transpose() * v_des_robot;
  //Vec3<float> v_robot = seResult.vWorld;
  //pretty_print(v_des_world, std::cout, "v des world");
  //pretty_print(v_robot, std::cout, "v_robot");

  for(int i = 0; i < 4; i++) {
    pFoot[i] = seResult.position + seResult.rBody.transpose() * (data._quadruped->getHipLocation(i) + data._legController->datas[i].p);
  }
  if(gait != &standing) {
    world_position_desired += dt * Vec3<float>(v_des_world[0], v_des_world[1], 0);//MODIFY-CHECK RUNTIME
    //world_position_desired =Vec3<float>(seResult.position[0], seResult.position[1], 0) + dt *iterationsBetweenMPC * Vec3<float>(v_des_world[0], v_des_world[1], 0);
    //if(gait != &jumping){
       //world_position_desired =Vec3<float>(seResult.position[0], seResult.position[1], 0) + dt *iterationsBetweenMPC * Vec3<float>(v_des_world[0], v_des_world[1], 0);
    //}
  }
  // some first time initialization
  if(firstRun){
    world_position_desired[0] = seResult.position[0];
    world_position_desired[1] = seResult.position[1];
    world_position_desired[2] = seResult.rpy[2];
    for(int i = 0; i < 4; i++){
      footSwingTrajectories[i].setHeight(0.05);
      footSwingTrajectories[i].setInitialPosition(pFoot[i]);
      footSwingTrajectories[i].setFinalPosition(pFoot[i]);
    }
    firstRun = false;
  }
  // foot placement
  for(int l = 0; l < 4; l++)
    swingTimes[l] = gait->getCurrentSwingTime(dtMPC, l);
  float side_sign[4] = {-1, 1, -1, 1};//.....................................?
  float interleave_y[4] = {-0.08, 0.08, 0.02, -0.02};
  //float interleave_gain = -0.13;
  float interleave_gain = -0.2;
  //float v_abs = std::fabs(seResult.vBody[0]);
  float v_abs = std::fabs(v_des_robot[0]);
  for(int i = 0; i < 4; i++){
    if(firstSwing[i]) {
      swingTimeRemaining[i] = swingTimes[i];
    } else {
      swingTimeRemaining[i] -= dt;
    }
    //if(firstSwing[i]) {//footSwingTrajectories[i].setHeight(.05);
    if(gait != &jumping){
      footSwingTrajectories[i].setHeight(0.06+0.06*sqrt(apla(0)*apla(0)+apla(1)*apla(1)));//JUMPINGUNCOM
      //footSwingTrajectories[i].setHeight(seResult.position[2]-0.23);
    }else{
      footSwingTrajectories[i].setHeight(seResult.position[2]-heightoffCoM-0.2);
    }
    Vec3<float> offset(0, side_sign[i] * (.065+data.userParameters->shiftleg), 0);
    Vec3<float> pRobotFrame = (data._quadruped->getHipLocation(i) + offset);
    pRobotFrame[1] += interleave_y[i] * v_abs * interleave_gain;
    float stance_time = gait->getCurrentStanceTime(dtMPC, i); //for heuristics
    Vec3<float> pYawCorrected = coordinateRotation(CoordinateAxis::Z, -_yaw_turn_rate* stance_time / 2) * pRobotFrame;
    Vec3<float> des_vel;
    des_vel[0] = _x_vel_des;
    des_vel[1] = _y_vel_des;
    des_vel[2] = 0.0;
    Vec3<float> Pf = seResult.position + seResult.rBody.transpose() * (pYawCorrected + des_vel * swingTimeRemaining[i]);//this produces variations
    //+ seResult.vWorld * swingTimeRemaining[i];//float p_rel_max = 0.35f;
    float p_rel_max = 0.3f*1.4142f;//0.3f;
    // Using the estimated velocity is correct//Vec3<float> des_vel_world = seResult.rBody.transpose() * des_vel;
    float pfx_rel = seResult.vWorld[0] * (.5 + _parameters->cmpc_bonus_swing) * stance_time +
      .03f*(seResult.vWorld[0]-v_des_world[0]) +
      (0.5f*seResult.position[2]/9.81f) * (seResult.vWorld[1]*_yaw_turn_rate);
    float pfy_rel = seResult.vWorld[1] * .5 * stance_time +
      .03f*(seResult.vWorld[1]-v_des_world[1]) +
      (0.5f*seResult.position[2]/9.81f) * (-seResult.vWorld[0]*_yaw_turn_rate);
    //pfx_rel = fminf(fmaxf(pfx_rel, -p_rel_max), p_rel_max);
    //pfy_rel = fminf(fmaxf(pfy_rel, -p_rel_max), p_rel_max);
    /*cout << "Leg: \n\n"<<i<< endl;
    cout << "pfx_relprev: "<<pfx_rel<< endl;
    cout << "pfy_relprev: "<<pfy_rel<< endl;*/

    float pf_norm=sqrt(pfx_rel*pfx_rel+pfy_rel*pfy_rel);
    if (pf_norm>p_rel_max){
      pfx_rel=pfx_rel*p_rel_max/pf_norm;
      pfy_rel=pfy_rel*p_rel_max/pf_norm;}
      /*cout << "pf_norm: "<<pf_norm<< endl;
      cout << "pfx_relpost: "<<pfx_rel<< endl;
    cout << "pfy_relpost: "<<pfy_rel<< endl;*/

    Pf[0] +=  pfx_rel;
    Pf[1] +=  pfy_rel;
    Pf[2] = apla(0)+apla(1)*Pf[0]+apla(2)*Pf[1]-0.003;//JUMPINGUNC
    //cout << "apla: "<<apla<< endl;
    //Pf[2] = 0.0;
    footSwingTrajectories[i].setFinalPosition(Pf);
    //cout << "pfx_rel: "<<pfx_rel<< endl;
  }
  // calc gait// load LCM leg swing gains
  Kp << 700, 0, 0,
     0, 700, 0,
     0, 0, 150;
  Kp_stance = 0*Kp;
  Kd << 7, 0, 0,
     0, 7, 0,
     0, 0, 7;
  Kd_stance = Kd;
  // gait
  Vec4<float> contactStates = gait->getContactState();
  Vec4<float> swingStates = gait->getSwingState();
  Timer updMPCTimer;
  //printf("jumping.get_current_gait_phase(): %i\n", jumping.getCurrentGaitPhase());

  updateMPCIfNeeded(mpcTable, data, omniMode, fmax(gait->getHoriz()-jump_state.cshift,nxtgait->getHoriz()));
  //updateMPCIfNeeded(mpcTable, data, omniMode,gait->getHoriz()-fmin(jump_state.cshift,gait->getHoriz()-nxtgait->getHoriz()));
  if (updMPCTimer.getMs()>0.01)
  printf("TIME update MPC If Needed: %.3f\n", updMPCTimer.getMs());
  //  StateEstimator* se = hw_i->state_estimator;
  Vec4<float> se_contactState(0,0,0,0);
  float pestCoMx=0;
  float pestCoMy=0;
  for(s16 foot = 0; foot < 4; foot++){
      if(swingStates[foot]==0 && gait != &jumping){//SWING??
          Wpla(foot,1)=pFoot[foot](0);
          Wpla(foot,2)=pFoot[foot](1);
          zfeet(foot)=pFoot[foot](2);
          pseudoInverse(Wpla,0.001,Wplainv);
          apla=(1-(float)data.userParameters->coef_fil)*apla+((float)data.userParameters->coef_fil)*Wplainv*zfeet;
      }
    pestCoMx+=(pFoot[foot](0))/4;
    pestCoMy+=(pFoot[foot](1))/4;
  }

  //cout << "zfeet: \n" << zfeet;
  //cout << "\nheighoffCoM: " << heightoffCoM;
  
  
  if(data.userParameters->det_terrn>0.5){
    Vec3<float> norpla(-apla(1),-apla(2),1);
    Vec3<float> nor2(apla(2),-apla(1),0);
    norpla.normalize();
    nor2.normalize();
    Vec3<float> nor3;
    nor3 = nor2.cross(norpla);
    Rpla[0].block(0,0,3,1)=nor3;
    Rpla[0].block(0,1,3,1)=nor2;
    Rpla[0].block(0,2,3,1)=norpla;
    /*cout << "\nnorpla: " << norpla;
    cout << "\nnor2: " << nor2;
    cout << "\nEye: " << Rpla*Rpla.transpose();
    cout << "\nEye2: " << Rpla.transpose()*Rpla;
    cout << "\nR: " << Rpla;*/
    rpyplane=rotationMatrixToRPY(Rpla[0]);
    rpyplane(2)=0;
    //rpyplane=rotationMatrixToRPY(Rpla);
    //float rollpla
    //float pitchpla
    //cout << "\nrollpla: " << rpyplane(0);
    //cout << "\npitchpla: " << rpyplane(1);
    Rpla[0]=rpyToRotMat(rpyplane);
    Rpla[1]=rpyToRotMat(rpyplane);
    Rpla[2]=rpyToRotMat(rpyplane);
    Rpla[3]=rpyToRotMat(rpyplane);
    heightoffCoM=apla(0)+apla(1)*pestCoMx+apla(2)*pestCoMy;
  }
  else{
    Matrix<float,3,1> RPYPlane1(0.8726,0,0);
    Matrix<float,3,1> RPYPlane2(-0.8726,0,0);
    Rpla[0]=rpyToRotMat(RPYPlane1);
    Rpla[1]=rpyToRotMat(RPYPlane2);
    Rpla[2]=rpyToRotMat(RPYPlane1);
    Rpla[3]=rpyToRotMat(RPYPlane2);
    /*Rpla[0].Identity();
    Rpla[1].Identity();
    Rpla[2].Identity();
    Rpla[3].Identity();*/
    heightoffCoM=0;
  }
  
  //Rpla[0].setIdentity();// JUMPINGUNCOM

  cout << "heightoffCoM: "<<heightoffCoM<< endl;
  cout << "pestCoMx: "<< pestCoMx<< endl;
  cout << "pestCoMy: "<< pestCoMy<< endl;
  cout << "seResult.position[2]: "<< seResult.position[2]<< endl;



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
      ball->color = {1.0, 0.2, 0.2, 0.5};
    }
#endif
  for(int foot = 0; foot < 4; foot++){
    float contactState = contactStates[foot];
    float swingState = swingStates[foot];
    if(swingState > 0){// foot is in swing
      if(firstSwing[foot]){
        firstSwing[foot] = false;
        footSwingTrajectories[foot].setInitialPosition(pFoot[foot]);
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
        }
      }
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
      Vec3<float> pDesFootWorld = footSwingTrajectories[foot].getPosition();
      Vec3<float> vDesFootWorld = footSwingTrajectories[foot].getVelocity();
      Vec3<float> pDesLeg = seResult.rBody * (pDesFootWorld - seResult.position) 
        - data._quadruped->getHipLocation(foot);
      Vec3<float> vDesLeg = seResult.rBody * (vDesFootWorld - seResult.vWorld);
      // Update for WBC
      pFoot_des[foot] = pDesFootWorld;
      vFoot_des[foot] = vDesFootWorld;
      aFoot_des[foot] = footSwingTrajectories[foot].getAcceleration();
      if(!data.userParameters->use_wbc){
        // Update leg control command regardless of the usage of WBIC
        data._legController->commands[foot].pDes = pDesLeg;
        data._legController->commands[foot].vDes = vDesLeg;
        data._legController->commands[foot].kpCartesian = Kp;
        data._legController->commands[foot].kdCartesian = Kd;
      }
    }
    else{ //foot is in stance
      firstSwing[foot] = true;
#ifdef DRAW_DEBUG_SWINGS
      auto* actualSphere = data.visualizationData->addSphere();
      actualSphere->position = pFoot[foot];
      actualSphere->radius = 0.02;actualSphere->color = {0.2, 0.2, 0.8, 0.7};
#endif
      Vec3<float> pDesFootWorld = footSwingTrajectories[foot].getPosition();
      Vec3<float> vDesFootWorld = footSwingTrajectories[foot].getVelocity();
      Vec3<float> pDesLeg = seResult.rBody * (pDesFootWorld - seResult.position) - data._quadruped->getHipLocation(foot);
      Vec3<float> vDesLeg = seResult.rBody * (vDesFootWorld - seResult.vWorld);
      //cout << "Foot " << foot << " relative velocity desired: " << vDesLeg.transpose() << "\n";
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
      //cout << "Foot " << foot << " force: " << f_ff[foot].transpose() << "\n";
      se_contactState[foot] = contactState;
      // Update for WBC
      //Fr_des[foot] = -f_ff[foot];
    }
  }
  // se->set_contact_state(se_contactState); todo removed
  data._stateEstimator->setContactPhase(se_contactState);
  // Update For WBC
  pBody_des[0] = world_position_desired[0];
  pBody_des[1] = world_position_desired[1];
  pBody_des[2] = _body_height;
  vBody_des[0] = v_des_world[0];
  vBody_des[1] = v_des_world[1];
  vBody_des[2] = 0.;
  aBody_des.setZero();
  pBody_RPY_des[0] = rpydesfin(0);//_roll_des;
  pBody_RPY_des[1] = rpydesfin(1); 
  pBody_RPY_des[2] = rpydesfin(2);//_yaw_des;//_yaw_des;
  vBody_Ori_des[0] = _roll_turn_rate;//;
  vBody_Ori_des[1] = 0.;
  vBody_Ori_des[2] = _yaw_turn_rate;//_yaw_turn_rate;
  contact_state = gait->getContactState();// END of WBC Update
  iterationCounter++;
}

template<>
void ConvexMPCLocomotion::run(ControlFSMData<double>& data) {
  (void)data;
  printf("call to old CMPC with double!\n");
}

void ConvexMPCLocomotion::updateMPCIfNeeded(int *mpcTable, ControlFSMData<float> &data, bool omniMode, int horLength) {
  //iterationsBetweenMPC = 30;
  if((iterationCounter % iterationsBetweenMPC) == 0){
    Timer MPCPrepTimer;
    auto seResult = data._stateEstimator->getResult();
    Vec3<float> v_des_robot(_x_vel_des, _y_vel_des,0);
    Vec3<float> v_des_world = omniMode ? v_des_robot : seResult.rBody.transpose() * v_des_robot;
    Matrix<float,3,1> orienvecdes;
    /*if (currently_jumping)_pitch_des = 0;//.1415*2/5;else_pitch_des=0;*/

    if(current_gait == 4 || firstRun){
      orienvecdes <<_roll_des,_pitch_des,_yaw_des;
      Matrix<float,4,1> quatdes;
      if(data.userParameters->ori_pla_body<0.5){quatdes=rpyToQuat(orienvecdes);}
      else{quatdes=rotationMatrixToQuaternion(rpyToRotMat(orienvecdes)*Rpla[0].transpose());}
      Matrix<float,4,1> quatprev(seResult.orientation.data());
      if ((quatprev-quatdes).norm()>(quatprev+quatdes).norm()){quatdes=-quatdes;}
      rpydesfin=quatToRPY(quatdes);
      float trajInitial[13] = {quatdes(0), quatdes(1), quatdes(2), quatdes(3),
        (float)stand_traj[0],(float)stand_traj[1],(float)_body_height,
        0,0,0,0,0,0};
      for(int i = 0; i < horLength; i++)
        for(int j = 0; j < 13; j++)
          trajAll[13*i+j] = trajInitial[j];
    }
    else if(bswitchgait){//1
      orienvecdes <<_roll_des,_pitch_des,_yaw_des;
      Matrix<float,4,1> quatdes;
      if(data.userParameters->ori_pla_body<0.5){quatdes=rpyToQuat(orienvecdes);}
      else{quatdes=rotationMatrixToQuaternion(rpyToRotMat(orienvecdes)*Rpla[0].transpose());}
      //printf("Inside first switch\n");
      float xStart = world_position_desired[0];
      float yStart = world_position_desired[1];
      Matrix<float,4,1> quatprev(seResult.orientation.data());
      if ((quatprev-quatdes).norm()>(quatprev+quatdes).norm()){quatdes=-quatdes;}
      rpydesfin=quatToRPY(quatdes);
      float trajInitial[13] = {quatdes(0),quatdes(1),quatdes(2),quatdes(3),
        xStart,yStart,(float)_body_height,
        0,0,0,
        v_des_world[0],v_des_world[1],0};
      for(int i = 0; i < horLength; i++){
        for(int j = 0; j < 13; j++)
          trajAll[13*i+j] = trajInitial[j];
        if(i != 0){
          trajAll[13*i + 4] = trajAll[13 * (i - 1) + 4] + dtMPC * v_des_world[0];
          trajAll[13*i + 5] = trajAll[13 * (i - 1) + 5] + dtMPC * v_des_world[1];
          //trajAll[13*i + 2] = trajAll[13 * (i - 1) + 2] + dtMPC * _yaw_turn_rate;
        }
      }
      bswitchgait = false;
    }
   else {
      //printf("Tracking previous MPC traj\n");
      update_traj(trajAll);
      //std::cout<<"trajAll: "<< *trajAll<<std::endl;
      orienvecdes <<_roll_des,_pitch_des,_yaw_des + dtMPC * horLength * _yaw_turn_rate;
      Matrix<float,4,1> quatdes;
      if(data.userParameters->ori_pla_body<0.5){quatdes=rpyToQuat(orienvecdes);}
      else{quatdes=rotationMatrixToQuaternion(rpyToRotMat(orienvecdes)*Rpla[0].transpose());}
      
      //quatdes=rotationMatrixToQuaternion(rpyToRotMat(orienvecdes)*Rpla[0].transpose());
      Matrix<float,4,1> quatprev(trajAll[13*(horLength-1)+0],trajAll[13*(horLength-1)+1],trajAll[13*(horLength-1)+2],trajAll[13*(horLength-1)+3]);
      if ((quatprev-quatdes).norm()>(quatprev+quatdes).norm()){quatdes=-quatdes;}
      rpydesfin=quatToRPY(quatdes);
      float trajEnd[13] = {quatdes(0),quatdes(1),quatdes(2),quatdes(3),
        world_position_desired[0] + dtMPC * horLength * v_des_world[0],
        world_position_desired[1] + dtMPC * horLength * v_des_world[1],
        (float)_body_height + dtMPC * horLength * (apla(1)*v_des_world[0]+apla(2)*v_des_world[1]),//(float)_body_height
        0,0,0,
        v_des_world[0],v_des_world[1],apla(1)*v_des_world[0]+apla(2)*v_des_world[1]};
      for(s16 itra = 0; itra < 13; itra++)
        trajAll[13*(horLength-1)+itra] = trajEnd[itra];//becuase kjoyst gives a roughr
      rpydesfin=quatToRPY(quatdes);
      /*std::cout<< "trajEnd: "<< *trajEnd<<std::endl;
      std::cout<< "quatprev: "<< quatprev <<std::endl;
      std::cout << "quatdes: " << quatdes << std::endl << "_yaw_des" << _yaw_des << std::endl;*/
      }

    //printf("Prep MPC Time: %.3f\n", MPCPrepTimer.getMs());
    Timer solveTimerSDMPC;
    solveDenseMPC(mpcTable, data, horLength);
    //printf("SOLVE TIME DenseMPC: %.3f\n", solveTimerSDMPC.getMs());
    /*cout<<"Fr_des1: \n"<<Fr_des[0]<<endl;
    cout<<"Fr_des2: \n"<<Fr_des[1]<<endl;
    cout<<"Fr_des3: \n"<<Fr_des[2]<<endl;
    cout<<"Fr_des4: \n"<<Fr_des[3]<<endl;*/
  }
}

void ConvexMPCLocomotion::solveDenseMPC(int *mpcTable, ControlFSMData<float> &data, int horLength) {
  auto seResult = data._stateEstimator->getResult();
  Matrix<double,13,1> QMatr;
  QMatr(0,0)=data.userParameters->Qwei1;
  QMatr.block(1,0,3,1)=data.userParameters->Qwei2;
  QMatr.block(4,0,3,1)=data.userParameters->Qwei3;
  QMatr.block(7,0,3,1)=data.userParameters->Qwei4;
  QMatr.block(10,0,3,1)=data.userParameters->Qwei5;
  double* Q = QMatr.data();
  double* weights = Q;
  float alpha = 4e-5; // make setting eventually
  float* p = seResult.position.data();
  float* v = seResult.vWorld.data();
  float* w = seResult.omegaWorld.data();
  float* q = seResult.orientation.data();
  float r[12];
  for(int i = 0; i < 12; i++)
    r[i] = pFoot[i%4][i/4]  - seResult.position[i/4];
  //printf("current posistion: %3.f %.3f %.3f\n", p[0], p[1], p[2]);
  if(alpha > 1e-4) {
    std::cout << "Alpha was set too high (" << alpha << ") adjust to 1e-5\n";
    alpha = 1e-5;
  }
  Vec3<float> pxy_act(p[0], p[1], 0);
  Vec3<float> pxy_des(world_position_desired[0], world_position_desired[1], 0);
  //Vec3<float> pxy_err = pxy_act - pxy_des;
  Vec3<float> vxy(seResult.vWorld[0], seResult.vWorld[1], 0);
  int summpcT=0;
  for (int indmp=0; indmp < 4*horLength; indmp++) summpcT+=mpcTable[indmp];
  Timer t1;
  dtMPC = dt * iterationsBetweenMPC;
  muv[0]=(float)data.userParameters->muest;
  muv[1]=(float)data.userParameters->muest;
  muv[2]=(float)data.userParameters->muest;
  muv[3]=(float)data.userParameters->muest;
  setup_problem(dtMPC,horLength,muv,120,summpcT, Rpla);
  update_solver_settings(_parameters->jcqp_max_iter, _parameters->jcqp_rho,
      _parameters->jcqp_sigma, _parameters->jcqp_alpha, _parameters->jcqp_terminate, _parameters->use_jcqp);
  //t1.stopPrint("Setup MPC");
  Timer t2;
  update_problem_data_floats(p,v,q,w,r,weights,trajAll,alpha,mpcTable);
  for(int leg = 0; leg < 4; leg++){
    Vec3<float> f;
    for(int axis = 0; axis < 3; axis++)
      f[axis] = get_solution(leg*3 + axis);
    //printf("[%d] %7.3f %7.3f %7.3f\n", leg, f[0], f[1], f[2]);
    f_ff[leg] = -seResult.rBody * f;
    Fr_des[leg] = f;// Update for WBC
  }
}