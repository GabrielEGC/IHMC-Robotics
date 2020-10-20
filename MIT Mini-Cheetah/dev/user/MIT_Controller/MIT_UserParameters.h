#ifndef PROJECT_MITUSERPARAMETERS_H
#define PROJECT_MITUSERPARAMETERS_H

#include "ControlParameters/ControlParameters.h"

class MIT_UserParameters : public ControlParameters {
public:
  MIT_UserParameters()
      : ControlParameters("user-parameters"),
        INIT_PARAMETER(cmpc_gait),
        INIT_PARAMETER(use_wbc),
        INIT_PARAMETER(cmpc_bonus_swing),
        INIT_PARAMETER(Kp_body),
        INIT_PARAMETER(Kd_body),
        INIT_PARAMETER(Kp_ori),
        INIT_PARAMETER(Kd_ori),
        INIT_PARAMETER(Kp_foot),
        INIT_PARAMETER(Kd_foot),
        INIT_PARAMETER(Kp_joint),
        INIT_PARAMETER(Kd_joint),
        INIT_PARAMETER(stance_legs),
        INIT_PARAMETER(gait_type),
        INIT_PARAMETER(gait_period_time),
        INIT_PARAMETER(gait_switching_phase),
        INIT_PARAMETER(gait_override),
        INIT_PARAMETER(xmaxfs),
        INIT_PARAMETER(ymaxfs),
        INIT_PARAMETER(yawmaxfs),
        INIT_PARAMETER(quasig),
        INIT_PARAMETER(Qwei2),
        INIT_PARAMETER(Qwei3),
        INIT_PARAMETER(Qwei4),
        INIT_PARAMETER(Qwei5),
        INIT_PARAMETER(Qweiu),
        INIT_PARAMETER(Npar1),
        INIT_PARAMETER(Npar2),
        INIT_PARAMETER(Npar3),
        INIT_PARAMETER(bo_hei_jum),
        INIT_PARAMETER(bo_height),
        INIT_PARAMETER(rollmaxfs),
        INIT_PARAMETER(pitchdesin),
        INIT_PARAMETER(coef_fil),
        INIT_PARAMETER(muest),
        INIT_PARAMETER(shiftleg),
        INIT_PARAMETER(slidshift),
        INIT_PARAMETER(det_terrn),
        INIT_PARAMETER(ori_pla_body),
        INIT_PARAMETER(ifleap),
        INIT_PARAMETER(dxjump),
        INIT_PARAMETER(dyjump),
        INIT_PARAMETER(noncopjump),
        INIT_PARAMETER(erlylnd),
        INIT_PARAMETER(ltelnd),
        INIT_PARAMETER(apla2align),
        INIT_PARAMETER(jumppg),
        INIT_PARAMETER(ywrtjmp),
        INIT_PARAMETER(savdat),
        INIT_PARAMETER(bckflstab),
        INIT_PARAMETER(sideflstab)
  {}
  DECLARE_PARAMETER(double, cmpc_gait);
  DECLARE_PARAMETER(double, use_wbc);
  DECLARE_PARAMETER(double, cmpc_bonus_swing);
  DECLARE_PARAMETER(Vec3<double>, Kp_body);
  DECLARE_PARAMETER(Vec3<double>, Kd_body);
  DECLARE_PARAMETER(Vec3<double>, Kp_ori);
  DECLARE_PARAMETER(Vec3<double>, Kd_ori);
  DECLARE_PARAMETER(Vec3<double>, Kp_foot);
  DECLARE_PARAMETER(Vec3<double>, Kd_foot);
  DECLARE_PARAMETER(Vec3<double>, Kp_joint);
  DECLARE_PARAMETER(Vec3<double>, Kd_joint);
  DECLARE_PARAMETER(double, stance_legs);
  // Gait Scheduler
  DECLARE_PARAMETER(double, gait_type);
  DECLARE_PARAMETER(double, gait_period_time);
  DECLARE_PARAMETER(double, gait_switching_phase);
  DECLARE_PARAMETER(double, gait_override);
  DECLARE_PARAMETER(double, xmaxfs);
  DECLARE_PARAMETER(double, ymaxfs);
  DECLARE_PARAMETER(double, yawmaxfs);
  DECLARE_PARAMETER(double, quasig);
  DECLARE_PARAMETER(Vec3<double>, Qwei2);
  DECLARE_PARAMETER(Vec3<double>, Qwei3);
  DECLARE_PARAMETER(Vec3<double>, Qwei4);
  DECLARE_PARAMETER(Vec3<double>, Qwei5);
  DECLARE_PARAMETER(Vec3<double>, Qweiu);
  DECLARE_PARAMETER(double, Npar1);
  DECLARE_PARAMETER(double, Npar2);
  DECLARE_PARAMETER(double, Npar3);
  DECLARE_PARAMETER(double, bo_hei_jum);
  DECLARE_PARAMETER(double, bo_height);
  DECLARE_PARAMETER(double, rollmaxfs);
  DECLARE_PARAMETER(double, pitchdesin);
  DECLARE_PARAMETER(double, coef_fil);
  DECLARE_PARAMETER(double, muest);
  DECLARE_PARAMETER(double, shiftleg);
  DECLARE_PARAMETER(double, slidshift);
  DECLARE_PARAMETER(double, det_terrn);
  DECLARE_PARAMETER(double, ori_pla_body);
  DECLARE_PARAMETER(double, ifleap);
  DECLARE_PARAMETER(double, dxjump);
  DECLARE_PARAMETER(double, dyjump);
  DECLARE_PARAMETER(double, noncopjump);
  DECLARE_PARAMETER(double, erlylnd);
  DECLARE_PARAMETER(double, ltelnd);
  DECLARE_PARAMETER(double, apla2align);
  DECLARE_PARAMETER(double, jumppg);
  DECLARE_PARAMETER(double, ywrtjmp);
  DECLARE_PARAMETER(double, savdat);
  DECLARE_PARAMETER(double, bckflstab);
  DECLARE_PARAMETER(double, sideflstab);
};

#endif //PROJECT_MITUSERPARAMETERS_H
