#include <eigen3/Eigen/Dense>
#ifndef _convexmpc_interface
#define _convexmpc_interface
#define K_MAX_GAIT_SEGMENTS 36

//#include "common_types.h"

#ifdef __cplusplus
#define EXTERNC extern "C"
#else
#define EXTERNC
#endif

struct problem_setup
{
  float dt;
  float mu[4];
  float f_max;
  int horizon;
  Eigen::Matrix<float,3,3> Rpla[4];
};

struct update_data_t
{
  float p[3];
  float v[3];
  float q[4];
  float w[3];
  float r[12];
  double weights[13];
  float traj[13*K_MAX_GAIT_SEGMENTS];
  float alpha;
  unsigned char gait[K_MAX_GAIT_SEGMENTS];
  unsigned char hack_pad[1000];
  int max_iterations;
  double rho, sigma, solver_alpha, terminate;
  int use_jcqp;
  float x_drag;
};

EXTERNC void setup_problem(double dt, int horizon, float* mu, double f_max, int sga = -1, Eigen::Matrix<float,3,3>* Rpla = nullptr);
EXTERNC double get_solution(int index);
EXTERNC void update_solver_settings(int max_iter, double rho, double sigma, double solver_alpha, double terminate, double use_jcqp);
EXTERNC void update_problem_data_floats(float* p, float* v, float* q, float* w,
                                        float* r, double* weights,
                                        float* state_trajectory, float alpha, int* gait);
EXTERNC void update_traj(float *trajSome);
#endif