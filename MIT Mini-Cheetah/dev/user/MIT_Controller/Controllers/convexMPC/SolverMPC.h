#ifndef _solver_mpc
#define _solver_mpc
#include <qpOASES.hpp>
#include "common_types.h"
#include "convexMPC_interface.h"
#include <iostream>
#include "../Printl/printheader.h"
#include <stdio.h>

using Eigen::Matrix;
using Eigen::Dynamic;
using Eigen::Quaternionf;
using Eigen::Quaterniond;

template <typename T>
T square(T a) {
  return a * a;
}

Matrix<fpt,4,3> QMATquat(Quaternionf q);
void solve_mpc(update_data_t* update, problem_setup* setup);
//void ct_ss_mats(Matrix<fpt,3,3> I_world, fpt m, Matrix<fpt,3,4> r_feet, Quaternionf q, Matrix<fpt,13,13>& A, Matrix<fpt,13,12>& B, Matrix<fpt,13,1>& C);
void resize_qp_mats(s16 horizon, int sga);
void getvNHacH(unsigned char* gait, int horiz);
void c2qp(Matrix<fpt,12,13> ACc, Matrix<fpt,12,12> Bc,problem_setup* setup);
void real_t_to_matrix(Matrix<fpt,Dynamic,1> &dst, qpOASES::real_t* src, s32 n_items);

extern void BDrdsim(float m, float rGx, float rGy, float rGz, float rfo1x, float
                    rfo1y, float rfo1z, Matrix<fpt,12,3>& BDrd);//[36]
extern void ACrdsim(float EA1, float EA2, float EA3, float I1, float I2, float
                    I3, float Ld1, float Ld2, float Ld3, float dr1, float dr2,
                    float dr3, float m, float r1, float r2, float r3, float rGx,
                    float rGy, float rGz, Matrix<fpt,12,13>& ACc);
/*extern void ADrdsim(float EA1, float EA2, float EA3, float I1, float I2, float
                    I3, float Ld1, float Ld2, float Ld3, float dr1, float dr2,
                    float dr3, float m, float r1, float r2, float r3, float rGx,
                    float rGy, float rGz, Matrix<fpt,12,12>& ADrd);//[144]
extern void CDrdsim(float EA1, float EA2, float EA3, float I1, float I2, float
                    I3, float Ld1, float Ld2, float Ld3, float dr1, float dr2,
                    float dr3, float m, float r1, float r2, float r3, float rGx,
                    float rGy, Matrix<fpt,12,1>& CDrd);//[12]*/
//void matrixLogRot(const Eigen::Matrix<float,3,3> & R, Eigen::Matrix<float,3,1> & omega);
void Quat2LogMap(const Eigen::Quaternionf & qpreR, Eigen::Matrix<float,3,1> & eta0);
Eigen::Matrix<float,3,1> quatToRPYsMPC(const Eigen::Quaternionf& q);
mfp* get_q_soln();
Matrix<fpt,Dynamic,1> get_Xsolin();
#endif