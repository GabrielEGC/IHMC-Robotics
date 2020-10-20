#ifndef CHEETAH_SOFTWARE_CONVEXMPCLOCOMOTION_H
#define CHEETAH_SOFTWARE_CONVEXMPCLOCOMOTION_H

#include <Controllers/FootSwingTrajectory.h>
#include <FSM_States/ControlFSMData.h>
#include "Gait.h"

#include <cstdio>

#define K_MAX_TRAJ_SEGMENTS 38//verify uinder 36

using std::cout;
using std::endl;

const Eigen::Matrix<float,3,3> I_bod((Eigen::Matrix<float,3,3>()<<.07f*1.4, 0, 0, 0, 0.26f, 0, 0, 0, 0.242f).finished());

template<typename T>
struct CMPC_Result {
  LegControllerCommand<T> commands[4];
  Vec4<T> contactPhase;
};

struct CMPC_Jump {
  int START_SEG = 0;
  int END_SEG = 0;
  bool jump_pending = false;
  bool jump_in_progress = false;
  bool pressed = false;
  int last_seg_seen = 0;
  int jump_wait_counter = 0;
  void debug(int seg) {
    (void)seg;
    //printf("[%d] pending %d running %d\n", seg, jump_pending, jump_in_progress);
  }

  void trigger_pressed(bool trigger) {
    if(!pressed && trigger) {
      if(!jump_pending && !jump_in_progress) {
        jump_pending = true;
        //printf("jump pending @ %d\n", seg);//potencial for give the robot time...
      }
    }
    pressed = trigger;
  }

  bool should_jump(int seg, bool m4lgcont) {
    debug(seg);
    START_SEG = 0;//Prev: Np2=0
    if(jump_pending && seg == START_SEG && m4lgcont) {//m4lgcont
      jump_pending = false;
      jump_in_progress = true;
      last_seg_seen = seg;
      return true;
    }
    if(jump_in_progress) {
      if(seg == END_SEG && seg != last_seg_seen) {
        jump_in_progress = false;
        last_seg_seen = seg;
        return false;
      }
      last_seg_seen = seg;
      return true;
    }
    last_seg_seen = seg;
    return false;
  }
};

class ConvexMPCLocomotion {
public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
  ConvexMPCLocomotion(float _dt, int _iterations_between_mpc, MIT_UserParameters* parameters);
  Eigen::Matrix<float,4,3> QMATquatv(Vec4<float> q){
    Eigen::Matrix<float,4,3> ouQ;
    ouQ << -q(1), -q(2), -q(3),
      q(0),q(3),-q(2),
      -q(3),q(0),q(1),
      q(2),-q(1),q(0);
    return (ouQ/2);};

  Eigen::Matrix<float,3,3> Eomtodrpy(float th, float psi){
    Eigen::Matrix<float,3,3> EomtodrpyM;
    EomtodrpyM << cos(psi)/cos(th), sin(psi)/cos(th), 0,
                  -sin(psi), cos(psi), 0,
                  cos(psi)*tan(th), sin(psi)*tan(th), 1;
    return (EomtodrpyM);};

  Eigen::Matrix<float,3,3> aplatoRpla(Vec3<float> apla2R){
    float s1,s2,c1,c2;
    Vec3<float> norpla(-apla2R(1),-apla2R(2),1);norpla.normalize();
    s2=-norpla(0);c2=sqrt(1-s2*s2);s1=norpla(1)/c2;c1=norpla(2)/c2;
    Eigen::Matrix<float,3,3> R;
    R << c2, 0, -s2,
        s2*s1, c1, c2*s1,
        s2*c1, -s1, c1*c2;
    return (R);};
  void RecompHeurFl(ControlFSMData<float> &data,Vec4<float>& quatdes,int horLength,int NP1, int NP2, int NP3, float& xFinal,float& yFinal,float& zFinal){
    auto seResult = data._stateEstimator->getResult();
    if((float)data.userParameters->ifleap>0.5){//testlater
      v_des_worldWBC = coordinateRotation(CoordinateAxis::Z, -seResult.rpy(2)) * Vec3<float>((float)data.userParameters->dxjump, (float)data.userParameters->dyjump, 0);
    }
    float _yaw_prev=seResult.rpy(2);
    Vec4<float> quatprev(seResult.orientation);
    _yaw_des=_yaw_prev+yawrotjump;//
    Vec3<float> orienvecdes1(_roll_des,_pitch_des,_yaw_des);
    if(data.userParameters->ori_pla_body<0.5){quatdes=rpyToQuat(orienvecdes1);}
    else{quatdes=rotationMatrixToQuaternion(rpyToRotMat(orienvecdes1)*(*Rplatarget)[0].transpose());}
    if (0>quatprev.transpose()*quatdes ){quatdes=-quatdes;}//cout<<"quatdes"<<quatdes<<endl;
    float xStart = world_position_desiredWBC[0];
    float yStart = world_position_desiredWBC[1];
    float zStart = seResult.position[2];
    xFinal = xStart+horLength*dtMPC * v_des_worldWBC[0];
    yFinal = yStart+horLength*dtMPC * v_des_worldWBC[1];
    zFinal=(*aplatarget)(0)+(*aplatarget)(1)*xFinal+(*aplatarget)(2)*yFinal+(float)data.userParameters->bo_height;
    float tLO=NP2*dtMPC;
    float tla=NP3*dtMPC;
    HeurFlighComp(xStart, yStart, zStart, xFinal, yFinal, zFinal, quatprev, quatdes, tLO, tla, NP1, NP2, NP3);
  }

  void RecompHeurFl(ControlFSMData<float> &data,Vec4<float>& quatdes,int horLength, int NP1, int NP2, int NP3){
    auto seResult = data._stateEstimator->getResult();
    if((float)data.userParameters->ifleap>0.5){//testlater
      v_des_worldWBC = coordinateRotation(CoordinateAxis::Z, -seResult.rpy(2)) * Vec3<float>((float)data.userParameters->dxjump, (float)data.userParameters->dyjump, 0);
    }
    float _yaw_prev=seResult.rpy(2);
    Vec4<float> quatprev(seResult.orientation);
    _yaw_des=_yaw_prev+yawrotjump;//
    Vec3<float> orienvecdes1(_roll_des,_pitch_des,_yaw_des);
    if(data.userParameters->ori_pla_body<0.5){quatdes=rpyToQuat(orienvecdes1);}
    else{quatdes=rotationMatrixToQuaternion(rpyToRotMat(orienvecdes1)*(*Rplatarget)[0].transpose());}
    if (0>quatprev.transpose()*quatdes){quatdes=-quatdes;}
    float xStart = (pFoot_des[0](0)+pFoot_des[1](0)+pFoot_des[2](0)+pFoot_des[3](0))/4;//world_position_desiredWBC[0]
    float yStart = (pFoot_des[0](1)+pFoot_des[1](1)+pFoot_des[2](1)+pFoot_des[3](1))/4;//world_position_desiredWBC[1]
    float zStart = (float)data.userParameters->bo_height+apla(0)+apla(1)*xStart+apla(2)*yStart;;//seResult.position[2];//I DON'T LIKE THIS.....
    float xFinal = xStart+horLength*dtMPC * v_des_worldWBC[0];
    float yFinal = yStart+horLength*dtMPC * v_des_worldWBC[1];
    float zFinal=(*aplatarget)(0)+(*aplatarget)(1)*xFinal+(*aplatarget)(2)*yFinal+(float)data.userParameters->bo_height;
    float tLO=NP2*dtMPC;
    float tla=NP3*dtMPC;
    HeurFlighComp(xStart, yStart, zStart, xFinal, yFinal, zFinal, quatprev, quatdes, tLO, tla, NP1, NP2, NP3);
    v_des_worldWBC<<(xFinal-xStart)/(tla-tLO),(yFinal-yStart)/(tla-tLO),(zFinal-zStart)/(tla-tLO)+9.81f*(tla-tLO)/2;//if00
  }

  void initialize();
  template<typename T>
  void run(ControlFSMData<T>& data);
  bool currently_jumping = false;
  Vec3<float> pBody_des;
  Vec3<float> vBody_des;
  Vec3<float> aBody_des;
  Vec3<float> pBody_RPY_des;
  Vec3<float> vBody_Ori_des;
  Vec3<float> pFoot_des[4];
  Vec3<float> vFoot_des[4];
  Vec3<float> aFoot_des[4];
  Vec3<float> Fr_des[4];
  Vec4<float> contact_state;
  Vec4<float> is_slidcntct;
  Vec3<float> v_des_worldWBC;
  Eigen::Matrix<float,Dynamic,Dynamic> Wpla;
  Eigen::Matrix<float,Dynamic,Dynamic> Wplainv;
  float heightoffCoM = 0;
  float heightoffCoM_prev = 0;
  Vec4<float> zfeet;
  Vec3<float> apla;
  Vec3<float>* aplatarget;
  Vec3<float> apla2;
  Vec3<float> apla2b;
  Eigen::Matrix<float,3,3> Rpla[4];
  Eigen::Matrix<float,3,3> Rpla2[4];
  Eigen::Matrix<float,3,3>* Rplatarget[4];
  float muv[4] = {0.4f,0.4f,0.4f,0.4f};
  Vec3<float> rpydesWBC;
  Vec3<float> drpydesWBC;
  float pestCoMx;
  float pestCoMy;
  Vec3<float> v_des_robot;
  float deltestlnd0=-10;
  int tFlcount0=0;
  int tFlcount=0;
  float deltestlndpre=0;
  float deltestlnd=0;
  int counttrust=0;
  bool chSfOri=true;

private:
  void _SetupCommand(ControlFSMData<float> & data);
  int Np1;
  int Np2;
  int Np3;
  int isparkjump=2;
  float isbckflstab=2.;
  float issideflstab=2.;
  float isjleap=0.;
  float isnoncopj=0.;
  int prevnoncopjump=2;
  float yawrotjump=0.;
  float _yaw_turn_rate;
  float _yaw_des = 0.; 
  float _roll_turn_rate;
  float _roll_des = 0.;
  float _pitch_des = 0.;
  float _x_vel_des = 0.;
  float _y_vel_des = 0.;
  float _body_height = 0.29;
  void HeurFlighComp(float xStart,float yStart,float zStart, float xFinal, float yFinal, float zFinal,const Vec4<float>& quatprev,const Vec4<float>& quatdes, float tLO, float tla, int hhL1, int hhL2, int hhL3);
  void recompute_timing(int iterations_per_mpc);
  void updateMPCIfNeeded(Gait *gait, int* mpcTable, ControlFSMData<float>& data, int horLength);
  void solveDenseMPC(int *mpcTable, ControlFSMData<float> &data, int horLength);
  int iterationsBetweenMPC;
  int horizonLength;
  int default_iterations_between_mpc;
  float dt;
  float dtMPC;
  int iterationCounter = 0;
  Vec3<float> f_ff[4];
  Vec4<float> swingTimes;
  FootSwingTrajectory<float> footSwingTrajectories[4];
  Vec3<float> footSwingTrajectoriesInit[4];
  ParkourGait jumpingpg, walkingpg, jumpingbfspg, jumpingsfspg;
  OffsetDurationGait trotting, bounding, pronking, jumpingodg, galloping, standing, trotRunning, walking, walking2, pacing;
  Gait* jumping;
  Mat3<float> Kp, Kd, Kp_stance, Kd_stance;
  bool firstRun = true;
  bool firstSwing[4];
  float swingTimeRemaining[4];
  int current_gait;
  int previous_gait;
  int gaitNumber;
  bool bswitchgait=false;
  bool flyingjump=false;
  Vec3<float> world_position_desiredWBC;
  Vec3<float> pFoot[4];
  //Vec3<float> vtrajCorr[4];
  CMPC_Result<float> result;
  float trajAll[13*K_MAX_TRAJ_SEGMENTS];
  float trajPrev[13*K_MAX_TRAJ_SEGMENTS];
  float trajHeur[13*K_MAX_TRAJ_SEGMENTS];
  float trajHeurJump[13*K_MAX_TRAJ_SEGMENTS];
  //u16 swfromjump=0;//CODE FOR NOT FORGETTING PREVIOUS TRAJ.,.I.E. CHAINING TRAJS.
  MIT_UserParameters* _parameters = nullptr;
  CMPC_Jump jump_state;
  bool pre0init = false;
  //vectorAligned<Vec12<double>> _sparseTrajectory;
};
#endif //CHEETAH_SOFTWARE_CONVEXMPCLOCOMOTION_H