#ifndef _solver_mpc
#define _solver_mpc

#include <qpOASES.hpp>
#include <eigen3/Eigen/Dense>
#include "common_types.h"
#include "convexMPC_interface.h"
#include <iostream>
#include "../Printl/printheader.h"
#include <stdio.h>

using Eigen::Matrix;
using Eigen::Dynamic;
using Eigen::Quaternionf;
using Eigen::Quaterniond;

template <class T>
T t_min(T a, T b)
{
    if(a<b) return a;
    return b;
}

template <class T>
T sq(T a)
{
    return a*a;
}

Matrix<fpt,4,3> QMATquat(Quaternionf q);
void solve_mpc(update_data_t* update, problem_setup* setup);
void ct_ss_mats(Matrix<fpt,3,3> I_world, fpt m, Matrix<fpt,3,4> r_feet, Quaternionf q, Matrix<fpt,13,13>& A, Matrix<fpt,13,12>& B, Matrix<fpt,13,1>& C);
void resize_qp_mats(s16 horizon, int sga);
void getvNHacH(unsigned char* gait, int horiz);
void c2qp(Matrix<fpt,13,13> Ac, Matrix<fpt,13,12> Bc, Matrix<fpt,13,1> Cc,problem_setup* setup);
void real_t_to_matrix(Matrix<fpt,Dynamic,1> &dst, qpOASES::real_t* src, s32 n_items);

mfp* get_q_soln();
Matrix<fpt,Dynamic,1> get_Xsolin();
#endif
