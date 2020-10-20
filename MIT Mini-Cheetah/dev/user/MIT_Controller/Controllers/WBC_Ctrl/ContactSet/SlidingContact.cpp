#include "SlidingContact.hpp"
#include <Utilities/Utilities_print.h>

// [ Fx, Fy, Fz ]
template <typename T>
SlidingContact<T>::SlidingContact(const FloatingBaseModel<T>* robot, int pt)
    : ContactSpec<T>(1), _max_Fz(1500.), _muv(0.4), _contact_pt(pt), _dim_U(2) {
  _Rpla.setIdentity();
  Contact::idx_Fz_ = 0;
  robot_sys_ = robot;
  Contact::Jc_ = DMat<T>(Contact::dim_contact_, cheetah::dim_config);
  Contact::JcDotQdot_ = DVec<T>::Zero(Contact::dim_contact_);
  Contact::Uf_ = DMat<T>::Zero(_dim_U, Contact::dim_contact_);

  Contact::Uf_(0, 0) = 1.;
  Contact::Uf_(1, 0) = -1.;// Upper bound of normal force
}

template <typename T>
SlidingContact<T>::~SlidingContact() {}

template <typename T>
bool SlidingContact<T>::_UpdateJc() {
  Contact::Jc_ = _Rpla.block(0,2,3,1).transpose()*robot_sys_->_Jc[_contact_pt]; ///VERIFYYYYYYY
  // Quat<T> quat = robot_sys_->_state.bodyOrientation;
  // Mat3<T> Rot = ori::quaternionToRotationMatrix(quat);
  // Contact::Jc_.block(0,3, 3,3) = Rot*Contact::Jc_.block(0,3,3,3);
  // Contact::Jc_.block(0,0, 3,3) = Rot.transpose()*Contact::Jc_.block(0,0,3,3);
  // pretty_print(Rot, std::cout, "body ori");
  // pretty_print(Contact::Jc_, std::cout, "Jc");
  return true;
}

template <typename T>
bool SlidingContact<T>::_UpdateJc_FBD() {// GO BACK CHANGE FOR TANH VERIFYYYYYYYYYYYYYYYYYYYY
  Vec2<T> vslctc=robot_sys_->_vGC[_contact_pt].block(0,0,2,1);// VERIFYYYYYYYYYYYYYYYYYYYY ADD RPla
  //std::cout<<"vslctc"<<robot_sys_->_vGC[_contact_pt]<<std::endl;
  T _muv4FBD;
  if (vslctc.norm()>0.04){
    _muv4FBD=_muv;
  }else{
    _muv4FBD=(T)0;
  }
  vslctc.normalize();
  Contact::JcFBD_ = Vec3<T>(-_muv4FBD*vslctc(0),-_muv4FBD*vslctc(1),1).transpose()*robot_sys_->_Jc[_contact_pt]; ///VERIFYYYYYYY
  return true;
}

template <typename T>
bool SlidingContact<T>::_UpdateJcDotQdot() {
  Contact::JcDotQdot_ = _Rpla.block(0,2,3,1).transpose()*robot_sys_->_Jcdqd[_contact_pt];///VERIFYYYYYYY
  // pretty_print(Contact::JcDotQdot_, std::cout, "JcDotQdot");
  return true;
}

template <typename T>
bool SlidingContact<T>::_UpdateUf() {
  Contact::Uf_ = DMat<T>::Zero(_dim_U, Contact::dim_contact_);
  Contact::Uf_(0, 0) = 1.;
  Contact::Uf_(1, 0) = -1.;// Upper bound of normal force
  Contact::Uf_=Contact::Uf_;///VERIFYYYYYYY//*_Rpla.block(0,2,3,1).transpose()
  /*std::cout << std::endl;
  std::cout << std::endl;
  std::cout << "Contact::Uf_" << Contact::Uf_ << std::endl;
  std::cout << "_Rpla" <<_Rpla << std::endl;*/
  return true;
}

template <typename T>
bool SlidingContact<T>::_UpdateInequalityVector() {
  Contact::ieq_vec_ = DVec<T>::Zero(_dim_U);
  Contact::ieq_vec_[1] = -_max_Fz;
  return true;
}

template class SlidingContact<double>;
template class SlidingContact<float>;
