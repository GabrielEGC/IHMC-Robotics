#ifndef LINK_SLID_TASK
#define LINK_SLID_TASK

// (X, Y, Z)
#include <WBC/Task.hpp>

template <typename T>
class FloatingBaseModel;

template <typename T>
class LinkSlidTask : public Task<T> {
 public:
  LinkSlidTask(const FloatingBaseModel<T>*, int link_idx,
              bool virtual_depend = true);
  virtual ~LinkSlidTask();
  void setRpla(Eigen::Matrix<T,3,3> Rpla) {_Rpla=Rpla;};
  DVec<T> _Kp, _Kd, _Kp_kin;
  Eigen::Matrix<T,3,3> _Rpla;

 protected:
  // Update op_cmd_
  virtual bool _UpdateCommand(const void* pos_des, const DVec<T>& vel_des,
                              const DVec<T>& acc_des);
  // Update Jt_
  virtual bool _UpdateTaskJacobian();
  // Update JtDotQdot_
  virtual bool _UpdateTaskJDotQdot();
  virtual bool _AdditionalUpdate() { return true; }

  const FloatingBaseModel<T>* robot_sys_;
  int link_idx_;
  bool virtual_depend_;
};

#endif
