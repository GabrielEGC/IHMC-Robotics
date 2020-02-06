#ifndef Cheetah_SINGLE_CONTACT
#define Cheetah_SINGLE_CONTACT

#include <Dynamics/FloatingBaseModel.h>
#include <Dynamics/Quadruped.h>
#include <WBC/ContactSpec.hpp>

template <typename T>
class SingleContact : public ContactSpec<T> {
 public:
  SingleContact(const FloatingBaseModel<T>* robot, int contact_pt);
  virtual ~SingleContact();

  void setMaxFz(T max_fz) { _max_Fz = max_fz; }
  void setMu(T muv) { _muv = muv; }
  void setRpla(Eigen::Matrix<T,3,3> Rpla) {_Rpla=Rpla;};

 protected:
  T _max_Fz;
  T _muv;
  Eigen::Matrix<T,3,3> _Rpla;
  int _contact_pt;
  int _dim_U;

  virtual bool _UpdateJc();
  virtual bool _UpdateJcDotQdot();
  virtual bool _UpdateUf();
  virtual bool _UpdateInequalityVector();

  const FloatingBaseModel<T>* robot_sys_;
};

#endif
