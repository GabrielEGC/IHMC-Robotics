#ifndef Cheetah_SLIDING_CONTACT
#define Cheetah_SLIDING_CONTACT

#include <Dynamics/FloatingBaseModel.h>
#include <Dynamics/Quadruped.h>
#include <WBC/ContactSpec.hpp>

template <typename T>
class SlidingContact : public ContactSpec<T> {
 public:
  SlidingContact(const FloatingBaseModel<T>* robot, int contact_pt);
  virtual ~SlidingContact();

  void setMaxFz(T max_fz) { _max_Fz = max_fz; }
  void setMu(T muv) { _muv = muv; }
  void setRpla(Eigen::Matrix<T,3,3> Rpla) {_Rpla=Rpla;};
  void setRFDesired(const DVec<T>& Fr_des) { ContactSpec<T>::Fr_des_(0) = Fr_des[2]; }//+15

 protected:
  T _max_Fz;
  T _muv;
  Eigen::Matrix<T,3,3> _Rpla;
  int _contact_pt;
  int _dim_U;

  virtual bool _UpdateJc();
  virtual bool _UpdateJc_FBD();
  virtual bool _UpdateJcDotQdot();
  virtual bool _UpdateUf();
  virtual bool _UpdateInequalityVector();

  const FloatingBaseModel<T>* robot_sys_;
};

#endif
