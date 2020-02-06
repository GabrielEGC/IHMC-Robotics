#ifndef CHEETAH_SOFTWARE_CONVEXMPCLOCOMOTION_H
#define CHEETAH_SOFTWARE_CONVEXMPCLOCOMOTION_H

#include <Controllers/FootSwingTrajectory.h>
#include <FSM_States/ControlFSMData.h>
#include <SparseCMPC/SparseCMPC.h>
#include "cppTypes.h"
#include "Gait.h"

#include <cstdio>

using Eigen::Array4f;
using Eigen::Array4i;
using Eigen::Dynamic;

template<typename T>
struct CMPC_Result {
  LegControllerCommand<T> commands[4];
  Vec4<T> contactPhase;
};

struct CMPC_Jump {
  int START_SEG = 0;
  int END_SEG = 0;
  int END_COUNT = 2;
  bool jump_pending = false;
  bool jump_in_progress = false;
  bool pressed = false;
  int seen_end_count = 0;
  int last_seg_seen = 0;
  int jump_wait_counter = 0;
  int cshift = 0;
  void debug(int seg) {
    (void)seg;
    //printf("[%d] pending %d running %d\n", seg, jump_pending, jump_in_progress);
  }

  void trigger_pressed(int seg, bool trigger) {
    (void)seg;//interesting... run function?
    if(!pressed && trigger) {
      if(!jump_pending && !jump_in_progress) {
        jump_pending = true;
        //printf("jump pending @ %d\n", seg);//potencial for give the robot time...
      }
    }
    pressed = trigger;
  }

  bool should_jump(int seg, int Np2) {
    debug(seg);
    START_SEG = Np2;
    if(jump_pending && seg == START_SEG) {
      jump_pending = false;
      jump_in_progress = true;
      //printf("jump begin @ %d\n", seg);
      seen_end_count = 0;
      last_seg_seen = seg;
      return true;
    }
    if(jump_in_progress) {
      if(seg == END_SEG && seg != last_seg_seen) {
        seen_end_count++;
        if(seen_end_count == END_COUNT) {
          seen_end_count = 0;
          cshift = 0;
          jump_in_progress = false;
          //printf("jump end @ %d\n", seg);
          last_seg_seen = seg;
          return false;
        }
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
  Eigen::Matrix<float,Dynamic,Dynamic> Wpla;
  Eigen::Matrix<float,Dynamic,Dynamic> Wplainv;
  float heightoffCoM = 0;
  Eigen::Matrix<float,4,1> zfeet;
  Eigen::Matrix<float,3,1> apla;
  Eigen::Matrix<float,3,1> apla2;
;
  Eigen::Matrix<float,3,3> Rpla[4];
  Eigen::Matrix<float,3,3> Rpla2[4];
  float muv[4] = {0.4f,0.4f,0.4f,0.4f};
  Vec3<float> rpyplane;
  Vec3<float> rpydesfin;

private:
  void _SetupCommand(ControlFSMData<float> & data);
  float _yaw_turn_rate;
  float _yaw_des = 0.;
  float _roll_turn_rate;
  float _roll_des = 0.;
  float _pitch_des = 0.;
  float _x_vel_des = 0.;
  float _y_vel_des = 0.;
  float _body_height = 0.29;
  float _body_height_running = 0.35;
  float _body_height_jumping = 0.45;
  void recompute_timing(int iterations_per_mpc);
  void updateMPCIfNeeded(int* mpcTable, ControlFSMData<float>& data, bool omniMode, int horLength);
  void solveDenseMPC(int *mpcTable, ControlFSMData<float> &data, int horLength);
  //void solveSparseMPC(int *mpcTable, ControlFSMData<float> &data);
  //void initSparseMPC();
  int iterationsBetweenMPC;
  int horizonLength;
  int default_iterations_between_mpc;
  float dt;
  float dtMPC;
  int iterationCounter = 0;
  Vec3<float> f_ff[4];
  Vec4<float> swingTimes;
  FootSwingTrajectory<float> footSwingTrajectories[4];
  OffsetDurationGait trotting, bounding, pronking, jumping, galloping, standing, trotRunning, walking, walking2, pacing, oneaerial;
  MixedFrequncyGait random, random2;
  Mat3<float> Kp, Kd, Kp_stance, Kd_stance;
  bool firstRun = true;
  bool firstSwing[4];
  float swingTimeRemaining[4];
  float stand_traj[6];
  int current_gait;
  int previous_gait;
  int gaitNumber;
  bool bswitchgait=false;
  Vec3<float> world_position_desired;
  Vec3<float> pFoot[4];
  CMPC_Result<float> result;
  float trajAll[13*36];
  MIT_UserParameters* _parameters = nullptr;
  CMPC_Jump jump_state;
  vectorAligned<Vec12<double>> _sparseTrajectory;
  SparseCMPC _sparseCMPC;
};
#endif //CHEETAH_SOFTWARE_CONVEXMPCLOCOMOTION_H