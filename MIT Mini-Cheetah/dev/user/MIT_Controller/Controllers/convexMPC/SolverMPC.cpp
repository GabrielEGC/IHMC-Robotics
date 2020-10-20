#include "SolverMPC.h"
#include <cmath>
#include <eigen3/unsupported/Eigen/MatrixFunctions>
#include <sys/time.h>
#include <Utilities/Timer.h>
//#include <Utilities/Utilities_print.h>
//#include <fstream>
//#define K_PRINT_EVERYTHING

using namespace std;
using std::cout;
using std::endl;
using Eigen::Array;
//qpOASES::real_t a;
Matrix<fpt,Dynamic,12> A_qp;
Matrix<fpt,Dynamic,Dynamic> B_qp;
Matrix<fpt,Dynamic,1> C_qp;
Matrix<fpt,Dynamic,1> Xsolin;
Matrix<fpt,12,12> Bdt;
Matrix<fpt,12,12> Adt;
Matrix<fpt,12,1> Cdt;
Matrix<fpt,25,25> ABc,expmm;
Matrix<fpt,Dynamic,Dynamic> S,SU;
Matrix<fpt,Dynamic,1> X_d;
Matrix<fpt,Dynamic,Dynamic> fmat;
Matrix<fpt,Dynamic,Dynamic> qH;
Matrix<fpt,Dynamic,1> qg;
Matrix<fpt,Dynamic,1> U_sol,U_d;
qpOASES::real_t* H_qpoases;
qpOASES::real_t* g_qpoases;
qpOASES::real_t* A_qpoases;
qpOASES::real_t* lb_qpoases;
qpOASES::real_t* q_soln;
u8 real_allocated = 0;

Matrix<int,Dynamic,1> vNH;
Array<Matrix<int,Dynamic,1>,Dynamic,1> acH;
int NModes;
//int* vNH;
//int** acH;
Matrix<fpt,Dynamic,1> get_Xsolin()
{return Xsolin;}
mfp* get_q_soln()
{return q_soln;}

void real_t_to_matrix(Matrix<fpt,Dynamic,1> &dst, qpOASES::real_t* src, s32 n_items){
  for(s32 i = 0; i < n_items; i++)
    dst(i) = src[i];
}

void matrix_to_real(qpOASES::real_t* dst, Matrix<fpt,Dynamic,Dynamic> src, s16 rows, s16 cols){
  s32 a = 0;
  for(s16 r = 0; r < rows; r++){
    for(s16 c = 0; c < cols; c++){
      dst[a] = src(r,c);
      a++;
}}}

void getvNHacH(unsigned char* gait, int horiz){
  vNH.resize(horiz, Eigen::NoChange);
  s16 HC = 0;
  for(int i = 1; i < horiz; i++){
    for(s16 j = 0; j < 4; j++){
      if (gait[4*i+j]!=gait[4*(i-1)+j]){
        vNH(HC)=i-1;
        HC++;
        break;
  }}}

  NModes=HC+1;
  vNH(HC)=horiz-1;
  vNH.conservativeResize(NModes, Eigen::NoChange);
  acH.conservativeResize(NModes, Eigen::NoChange);
  for(int HCi = 0; HCi<HC+1; HCi++){
    acH(HCi).conservativeResize(4,Eigen::NoChange);
    s16 lki = 0;
    for(int ki = 0;ki < 4; ki++){
      if ((int)gait[4*vNH(HCi)+ki]==1){
        acH(HCi)(lki)=ki;
        lki++;
    }}
    acH(HCi).conservativeResize(lki,Eigen::NoChange);
  }
}

void c2qp(Matrix<fpt,12,13> ACc, Matrix<fpt,12,12> Bc, problem_setup* setup){
  fpt dt=setup->dt;
  s16 horizon=setup->horizon;
  Matrix<fpt,3,3> Rplac2qp[4];
  for(s16 i=0;i<4;i++){
    Rplac2qp[i]=setup->Rpla[i];}

  ABc.setZero();
  ABc.block(0,0,12,12) = ACc.block(0,0,12,12);
  ABc.block(0,12,12,12) = Bc;
  ABc.block(0,24,12,1) = ACc.block(0,12,12,1);
  ABc = dt*ABc;
  expmm = ABc.exp();
  Adt = expmm.block(0,0,12,12);
  Bdt = expmm.block(0,12,12,12);
  Cdt = expmm.block(0,24,12,1);
  #ifdef K_PRINT_EVERYTHING
    cout<<"Adt: \n"<<Adt<<"\nBdt:\n"<<Bdt<<"\nCdt:\n"<<Cdt<<endl;
  #endif
  if(horizon > (K_MAX_GAIT_SEGMENTS-1)) {
    throw std::runtime_error("horizon is too long!");
  }
  Matrix<fpt,12,12> powerMats[30];
  Matrix<fpt,12,1> sumpowerMats[30];
  powerMats[0].setIdentity();
  sumpowerMats[0]=Cdt;
  for(int i = 1; i < horizon+1; i++) {
    powerMats[i] = Adt * powerMats[i-1];
    sumpowerMats[i]=Adt*sumpowerMats[i-1]+Cdt;
  }
  //Eigen::Array<Matrix<fpt,12,Dynamic>,NModes,1> BdtH;//wondering if initialization is needed
  Eigen::Array<Matrix<fpt,12,Dynamic>,Dynamic,1> BdtH;
  //float* BdtH[NModes];
  BdtH.resize(NModes,1);
  fpt mu = 1.f/setup->mu[0];
  Matrix<fpt,6,3> f_block;
  f_block <<  mu, 0,  1.f,
    -mu, 0,  1.f,
    0,  mu, 1.f,
    0, -mu, 1.f,
    0,   0, 1.f,
    0,   0, -1.f;
    /*f_block <<  1.f, 0, setup->mu[0],
    -1.f, 0,  setup->mu[0],
    0,  1.f, setup->mu[0],
    0, -1.f, setup->mu[0],
    0,   0, 1.f,
    0,   0, -1.f;*/
  /*std::cout << "Rplac2qp0: "<< Rplac2qp[0] <<std::endl;
  std::cout << "Rplac2qp1: "<< Rplac2qp[1] <<std::endl;
  std::cout << "Rplac2qp2: "<< Rplac2qp[2] <<std::endl;
  std::cout << "Rplac2qp3: "<< Rplac2qp[3] <<std::endl;
  std::cout << "f_block: "<< f_block <<std::endl;*/
  u32 ifm=0;
    for (u16 iN = 0; iN < vNH(0)+1;iN++){
      u16 nacHN=acH(0).rows();
      fpt U_dz=88.29f/nacHN;
      for (u16 N3 = 0; N3 < nacHN;N3++){
        fmat.block(ifm*6,ifm*3,6,3) = f_block*(Rplac2qp[acH(0)(N3)].transpose());
        U_d(ifm*3+2) = U_dz;//U_d.block(ifm*3,1,3,1) = (0,0,44.145f);
        ifm++;
    }}
  for (u16 N = 1; N < NModes;N++){
    for (u16 iN = vNH(N-1)+1; iN < vNH(N)+1;iN++){  
      u16 nacHN=acH(N).rows();
      fpt U_dz=88.29f/nacHN;
      for (u16 N3 = 0; N3 < nacHN;N3++){
        fmat.block(ifm*6,ifm*3,6,3) = f_block*(Rplac2qp[acH(N)(N3)].transpose());
        U_d(ifm*3+2) = U_dz;//U_d.block(ifm*3,1,3,1) = (0,0,44.145f);
        ifm++;
  }}}
  //for(s16 i = 0; i < fmat.cols(); i++)
  //  U_d(i) = (i%3==2)*44.145f;//redefineeeee
  /*cout << "\nfmat.rows()/6): " << fmat.rows()/6 << endl;
  cout << "\nifm: " << ifm << endl;*/
  for (u16 N = 0; N < NModes;N++){
    BdtH(N).resize(12,acH(N).rows()*3);
    for (u16 N3 = 0; N3 < acH(N).rows();N3++){
      BdtH(N).block(0,3*N3,12,3)=Bdt.block(0,(acH(N)(N3))*3,12,3);
  }}
    A_qp.block(0,0,12,12) = Adt;
    C_qp.block(0,0,12,1) = Cdt;
    B_qp.block(0,0,12,3*acH(0).rows()) = BdtH(0);
    int aB_qp =0;
    aB_qp+=3*acH(0).rows();
  for(s16 r = 1; r < vNH(0)+1; r++){
    A_qp.block(12*r,0,12,12) = Adt*A_qp.block(12*(r-1),0,12,12);
    C_qp.block(12*r,0,12,1) = Adt*C_qp.block(12*(r-1),0,12,1)+Cdt;  
    B_qp.block(12*r,0,12,aB_qp) = Adt*B_qp.block(12*(r-1),0,12,aB_qp);
    B_qp.block(12*r,aB_qp,12,3*acH(0).rows()) = BdtH(0);
    aB_qp+=3*acH(0).rows();
  }
  for (u16 N = 1; N < NModes;N++){
    for(s16 r = vNH(N-1)+1; r < vNH(N)+1; r++){
      A_qp.block(12*r,0,12,12) = Adt*A_qp.block(12*(r-1),0,12,12);
      C_qp.block(12*r,0,12,1) = Adt*C_qp.block(12*(r-1),0,12,1)+Cdt;
      B_qp.block(12*r,0,12,aB_qp) = Adt*B_qp.block(12*(r-1),0,12,aB_qp);
      B_qp.block(12*r,aB_qp,12,3*acH(N).rows()) = BdtH(N);
      aB_qp+=3*acH(N).rows();
  }}//cout<<"AQP_err:\n"<<(A_qp-A_qp2).sum()<<"\nBQP_err:\n"<<(B_qp-B_qp2).sum()<<"\nCQP_err:\n"<<(C_qp-C_qp2).sum()<<endl;
#ifdef K_PRINT_EVERYTHING
  cout<<"AQP:\n"<<A_qp<<"\nBQP:\n"<<B_qp<<"\nCQP:\n"<<C_qp<<endl;
#endif
}

void resize_qp_mats(s16 horizon, int sga){
  int mcount = 0;
  int h2 = horizon*horizon;
  int sghor = sga*horizon;
  int sga2 = sga*sga;

  A_qp.resize(12*horizon, Eigen::NoChange);
  mcount += 12*horizon;//??
  U_sol.resize(3*sga, Eigen::NoChange);
  mcount += 3*sga;//??
  U_d.resize(3*sga, Eigen::NoChange);
  mcount += 3*sga;
  Xsolin.resize(12*horizon, Eigen::NoChange);
  mcount += 12*horizon;
  B_qp.resize(12*horizon, 3*sga);
  mcount += 12*3*sghor;
  C_qp.resize(12*horizon, Eigen::NoChange);
  mcount += 12*horizon;
  S.resize(12*horizon, 12*horizon);
  mcount += 12*12*h2;
  SU.resize(3*sga, 3*sga);
  mcount += 3*3*sga2;
  X_d.resize(12*horizon, Eigen::NoChange);
  mcount += 12*horizon;
  fmat.resize(6*sga, 3*sga);
  mcount += 6*3*sga2;
  qH.resize(3*sga, 3*sga);
  mcount += 3*3*sga2;
  qg.resize(3*sga, Eigen::NoChange);
  mcount += 3*sga;
  //printf("realloc'd %d floating point numbers.\n",mcount);
  mcount = 0;
  A_qp.setZero();
  B_qp.setZero();
  C_qp.setZero();
  S.setZero();
  SU.setZero();
  X_d.setZero();
  fmat.setZero();
  qH.setZero();
  U_sol.setZero();
  U_d.setZero();
  Xsolin.setZero();
  //TODO: use realloc instead of free/malloc on size changes
  if(real_allocated){
    free(H_qpoases);
    free(g_qpoases);
    free(A_qpoases);
    free(lb_qpoases);
    free(q_soln);
  }

  H_qpoases = (qpOASES::real_t*)malloc(3*3*sga2*sizeof(qpOASES::real_t));
  mcount += 3*3*sga2;
  g_qpoases = (qpOASES::real_t*)malloc(3*sga*sizeof(qpOASES::real_t));
  mcount += 3*sga;
  A_qpoases = (qpOASES::real_t*)malloc(3*6*sga2*sizeof(qpOASES::real_t));
  mcount += 3*6*sga2;
  lb_qpoases = (qpOASES::real_t*)malloc(6*sga*sizeof(qpOASES::real_t));
  mcount += 6*sga;
  q_soln = (qpOASES::real_t*)malloc(3*sga*sizeof(qpOASES::real_t));
  mcount += 3*sga;
  real_allocated = 1;
  //printf("malloc'd %d floating point numbers.\n",mcount);
#ifdef K_DEBUG
  printf("RESIZED MATRICES FOR HORIZON: %d\n",horizon);
#endif
}

Matrix<fpt,12,13> Ac_ct;
Matrix<fpt,12,12> B_ct_r;

const Matrix<fpt,3,1> I_D(.07f*1.4, 0.26f, 0.242f);//Id << 0.3f, 2.1f, 2.1f; // DH
/*const fpt I1 = I_D(0);
const fpt I2 = I_D(1);
const fpt I3 = I_D(2);*/
const Matrix<fpt,3,3> I_bd((Eigen::Matrix3f()<<I_D(0), 0, 0, 0, I_D(1), 0, 0, 0, I_D(2)).finished());
const fpt Mss=9.0f; //fpt Mss = 50.236; //DH

Quaternionf rsq;
Matrix<fpt,3,1> rsp,rsv,rsL0;
Matrix<fpt,3,4> rsr_feet;
Matrix<fpt,3,3> rsRrot;


void solve_mpc(update_data_t* update, problem_setup* setup){
  const fpt I1 = I_D(0);
  const fpt I2 = I_D(1);
  const fpt I3 = I_D(2);
  fpt EA1;
  fpt EA2;
  fpt EA3;
  fpt r1;
  fpt r2;
  fpt r3;
  fpt Ld1;
  fpt Ld2;
  fpt Ld3;
  fpt dr1;
  fpt dr2;
  fpt dr3;
  fpt rGx;
  fpt rGy;
  fpt rGz;
  //rs.set(update->p, update->v, update->q, update->w, update->r);
  rsq.w() = update->q[0];
  rsq.x() = update->q[1];
  rsq.y() = update->q[2];
  rsq.z() = update->q[3];
  rsRrot = rsq.toRotationMatrix();
  rsL0 = rsRrot * I_bd * rsRrot.transpose()*Matrix<fpt,3,1>(update->w[0],update->w[1],update->w[2]);//rsL0
  //x_0traj<<rsq.w(), rsq.x(), rsq.y(), rsq.z(), rs.p, rsL0, rs.v;
  for(u8 i = 0; i < 3; i++){
    rsp(i) = update->p[i];
    rsv(i) = update->v[i];
  }
  for(u8 rs = 0; rs < 3; rs++)
    for(u8 c = 0; c < 4; c++)
      rsr_feet(rs,c) = update->r[rs*4 + c];
  
  Eigen::Matrix<float,3,1> EAalt=quatToRPYsMPC(rsq);//cout<<"EAalt: \n"<<EAalt<<endl; 
  EA1 = EAalt(2);
  EA2 = EAalt(1);
  EA3 = EAalt(0);
  r1=rsp(0);//update->p[0];
  r2=rsp(1);//update->p[1];
  r3=rsp(2);//update->p[2];
  Ld1=rsL0(0)-update->seLcomp[0];
  Ld2=rsL0(1)-update->seLcomp[1];
  Ld3=rsL0(2)-update->seLcomp[2];
  dr1=rsv(0);//update->v[0];
  dr2=rsv(1);//update->v[1];
  dr3=rsv(2);//update->v[2];
  rGx=r1;
  rGy=r2;
  rGz=r3;
  Matrix<fpt,12,1> x_0;//NOT ITS PLACE BUT OK
  x_0<<0,0,0, rsp, rsL0+Mss*(rsp-Matrix<float,3,1>(rGx,rGy,rGz)).cross(rsv), rsv;//VERIFY 0 ,0,0 FROM BEGGINING
  Matrix<fpt,13,1> x_0traj;
  x_0traj<<rsq.w(), rsq.x(), rsq.y(), rsq.z(), rsp, rsL0, rsv;
  printf("   PROBLEM DATA  \n");
  print_problem_setup(setup);
  print_named_arraytr("Gait",update->gait,4,setup->horizon); 
  
  /*print_named_array("\nGait",update->gait,setup->horizon,4); */
  /////printf("\n   PROBLEM DATA  \n");
  /////print_problem_setup(setup);
  /////printf("    ROBOT DATA   \n");
  /////cout<<"Robot State:"<<endl<<"Position: "<<rsp<<"\nVelocity: "<<rsv<<"\nAngular Momentum: "<<rsL0<<"\nRotation : "<<rsRrot<<"\nFoot Locations\n"<<rsr_feet<<"\nI_bdy\n"<<I_bd<<endl;//<<"\nOrientation\n"<<rsq
  /////cout<<"Initial state (MPC vars): \n"<<x_0.transpose()<<endl;
  /////cout<<"Initial state (CMCL vars): \n"<<x_0traj.transpose()<<endl;
  //print_update_data(update,setup->horizon);
#ifdef K_PRINT_EVERYTHING
    printf("-----------------\n");
    printf("   PROBLEM DATA  \n");
    printf("-----------------\n");
    print_problem_setup(setup);
    printf("-----------------\n");
    printf("    ROBOT DATA   \n");
    printf("-----------------\n");
    cout<<"Robot State:"<<endl<<"Position: "<<rsp<<"\nVelocity: "<<rsv<<"\nAngular Momentum: "<<rsL0<<"\nRotation : "<<rsRrot<<"\nOrientation\n"<<rsq<<"\nFoot Locations\n"<<rsr_feet<<"\nI_bdy\n"<<I_bd<<endl;
    print_update_data(update,setup->horizon);
#endif
  ACrdsim(EA1, EA2, EA3, I1, I2, I3, Ld1, Ld2, Ld3, dr1, dr2, dr3, Mss, r1, r2, r3, rGx, rGy, rGz, Ac_ct);
  for(s16 b = 0; b < 4; b++){
    Matrix<fpt,12,3> BDrd;
    BDrdsim(Mss, rGx, rGy, rGz, rsr_feet(0,b)+rsp(0), rsr_feet(1,b)+rsp(1), rsr_feet(2,b)+rsp(2), BDrd);
    B_ct_r.block(0,b*3,12,3) = BDrd;
  }
#ifdef K_PRINT_EVERYTHING
    //cout<<"Initial state: \n"<<x_0<<endl;
    //cout<<"World Inertia: \n"<<I_world<<endl;
    cout<<"A CT: \n"<<A_ct<<endl;
    cout<<"B CT (simplified): \n"<<B_ct_r<<endl;
#endif
  Timer TimegetvNHacH;
  getvNHacH(update->gait, setup->horizon);
  printf("TimegetvNHacH: %.3f ms\n", TimegetvNHacH.getMs());
  Timer Timc2qp;
  c2qp(Ac_ct,B_ct_r,setup);//****needed??
  printf("Timc2qp: %.3f ms\n", Timc2qp.getMs());
  Matrix<fpt,12,1> full_weight;
  Matrix<fpt,3,1> full_uweight;
  for(u8 i = 0; i < 12; i++)
    full_weight(i) = (float)update->weights[i];
  S.diagonal() = full_weight.replicate(setup->horizon,1);
  for(u8 i = 0; i < 3; i++)
    full_uweight(i) = (float)update->uweights[i];
  SU.diagonal() = full_uweight.replicate(SU.cols()/3,1);
  for(s16 i = 0; i < setup->horizon; i++){
    Matrix<float,3,1> drdtraj(update->traj[13*i+10],update->traj[13*i+11],update->traj[13*i+12]);
    Matrix<float,3,1> rdtraj(update->traj[13*i+4],update->traj[13*i+5],update->traj[13*i+6]);
    for(s16 j = 9; j < 12; j++)
      X_d(12*i+j,0) = drdtraj(j-9);
    for(s16 j = 3; j < 6; j++)
      X_d(12*i+j,0) = rdtraj(j-3);
    Matrix<float,3,1> Lwcorr(Mss*(rdtraj-Matrix<float,3,1>(rGx,rGy,rGz)).cross(drdtraj));
    for(s16 j = 6; j < 9; j++)
      X_d(12*i+j,0) = update->traj[13*i+j+1]+Lwcorr(j-6);
    Quaternionf qdt(update->traj[13*i+0],update->traj[13*i+1],update->traj[13*i+2],update->traj[13*i+3]);
    /*Matrix<float,3,3> Rqdt=rs.R.transpose()*qdt.toRotationMatrix();Matrix<float,3,1> omega0dt=;*/
    Quaternionf qpreR=rsq.inverse()*qdt;
    qpreR.normalize();
    Matrix<float,3,1> eta0dt;
    Quat2LogMap(qpreR,eta0dt);//matrixLogRot(Rqdt, eta0dt);//Matrix<float,3,3> Rqdt=qpreR.normalized().toRotationMatrix();//Matrix<float,3,1> eta0dt;
    for(s16 j = 0; j < 3; j++)
      X_d(12*i+j,0) = eta0dt(j);//Eigen::Matrix<float,3,1> EA(rs.R.eulerAngles(0,1,2));
  }
  //cout<<"XD:\n"<<X_d.block(0,0,12,1).transpose()<<endl;//cout<<"UD:\n"<<U_d.transpose()<<endl;//cout << "\nRot trans in MPC: " << setup->Rpla.transpose() << endl;
  qH = 2*(B_qp.transpose()*S*B_qp+ SU);//wtf//update->alpha*eye_12h*0.001//alpha = 4e-5
  qg = 2*B_qp.transpose()*S*(A_qp*x_0 - (X_d-C_qp))-2*SU*U_d;
  matrix_to_real(H_qpoases,qH,qH.rows(), qH.cols());
  matrix_to_real(g_qpoases,qg,qg.rows(), 1);
  matrix_to_real(A_qpoases,fmat,fmat.rows(), fmat.cols());
  Timer itlbtime;
  float fmaxQP = - (setup->f_max);
  for(s16 i = 0; i < fmat.rows(); i++)
    lb_qpoases[i] = (i%6==5)*fmaxQP;
  //printf("Timeitlbtime: %.4f ms\n", itlbtime.getMs());
  s16 num_constraints = fmat.rows();
  s16 num_variables = fmat.cols();
  qpOASES::int_t nWSR = 100;
  qpOASES::QProblem problem_red (num_variables, num_constraints);
  qpOASES::Options op;
  op.setToMPC();
  op.printLevel = qpOASES::PL_NONE;
  problem_red.setOptions(op);//int_t nWSR = 50000;
  Timer solve_timer;
  int rval = problem_red.init(H_qpoases, g_qpoases, A_qpoases, NULL, NULL, lb_qpoases, NULL, nWSR);(void)rval;//int rval = problem_red.init(&H_SDM, g_qpoases,&A_SDM, NULL, NULL, lb_qpoases, NULL, nWSR, 0);//q_soln_guess
  int rval2 = problem_red.getPrimalSolution(q_soln);
  if(rval2 != qpOASES::SUCCESSFUL_RETURN)
    printf("failed to solve!\n");
  //printf("solve time: %.3f ms, size %d, %d\n", solve_timer.getMs(), num_variables, num_constraints);
  Timer timer_rtm;
  real_t_to_matrix(U_sol,q_soln,U_sol.rows());
  //cout<<"U_sol.rows: \n"<<U_sol.rows()<<endl;  //cout<<"U_sol: "<<U_sol.transpose()<<endl;  //cout<<"q_soln: "<<q_soln[0]<<", "<<q_soln[1]<<", "<<q_soln[2]<<", "<<q_soln[3]<<endl;
  if(U_sol.rows()>0){
    Xsolin=A_qp*x_0+B_qp*U_sol+C_qp;
  }else{
    Xsolin=A_qp*x_0+C_qp;
  }
  //printf("CompTimer2m: %.3f ms.\n", timer_rtm.getMs());
  for(s16 i = 0; i < (setup->horizon); i++){//SHIFTED 12 and (setup->horizon-1)
    Matrix<float,3,1> drpred(Xsolin(12*i+9),Xsolin(12*i+10),Xsolin(12*i+11));
    Matrix<float,3,1> rpred(Xsolin(12*i+3),Xsolin(12*i+4),Xsolin(12*i+5));
    for(s16 j = 9; j < 12; j++)
      update->traj[13*i+j+1] = drpred(j-9);
    for(s16 j = 3; j < 6; j++)
      update->traj[13*i+j+1] = rpred(j-3);
    Matrix<float,3,1> Lwpred(Mss*(rpred-Matrix<float,3,1>(rGx,rGy,rGz)).cross(drpred));
    for(s16 j = 6; j < 9; j++)
      update->traj[13*i+j+1] = Xsolin(12*i+j)-Lwpred(j-6);
    Matrix<float,3,1> etapred(Xsolin(12*i+0),Xsolin(12*i+1),Xsolin(12*i+2));
    float thpr=etapred.norm();
    Quaternionf qpreR;
    if (thpr>10e-05){
      qpreR.w()=cos(thpr/2);
      qpreR.vec()=etapred*sin(thpr/2)/thpr;
    }else{
      qpreR.setIdentity();
    }
    Quaternionf qpredf=rsq*qpreR;
    update->traj[13*i] = qpredf.w();
    update->traj[13*i+1] = qpredf.x();
    update->traj[13*i+2] = qpredf.y();
    update->traj[13*i+3] = qpredf.z();
  }
  /*for(s16 itra = 0; itra < 13; itra++)
    update->traj[13*(setup->horizon-1)+itra] = update->traj[13*(setup->horizon-2)+itra];
  */
  /////print_named_array("predicted trajectory: \n",update->traj,setup->horizon,13);
  //pretty_print(Xsolin, std::cout, "Xsolin");
  /*ofstream DatXsol;
  DatXsol.open("/home/ggarcia/DXsolin.txt",std::ios_base::app);
  DatXsol << "];iX=iX+1;X{iX}=[: " << setup->horizon << endl;
  DatXsol << Xsolin << endl;
  DatXsol.close();
  ofstream DatX_d;
  DatX_d.open("/home/ggarcia/DX_d.txt",std::ios_base::app);
  DatX_d << "];iX=iX+1;X{iX}=[: " << setup->horizon << endl;
  DatX_d << X_d << endl;
  DatX_d.close(); */
  //cout<<"\n"<<Xsolin<<endl;
  //cout<<"\n"<<X_d<<endl;
  #ifdef K_PRINT_EVERYTHING
    //cout<<"fmat:\n"<<fmat<<endl;
  #endif
}

void Quat2LogMap(const Eigen::Quaternionf & qpreR, Eigen::Matrix<float,3,1> & eta0) {//theta = acos( (Trace(R) - 1)/2 );
   eta0=qpreR.vec();
   float th=2*atan2(eta0.norm(),qpreR.w());
   if (th > 10e-5) {
      eta0*=th/(sin(th/2));
   }else {
      eta0*=2;
   } 
}

Eigen::Matrix<float,3,1> quatToRPYsMPC(const Eigen::Quaternionf& q) {
  Eigen::Matrix<float,3,1> rpy;
  float as = std::min(-2. * (q.x() * q.z() - q.w() * q.y()), .99999);
  rpy(2) =std::atan2(2 * (q.x() * q.y() + q.w() * q.z()),square(q.w()) + square(q.x()) - square(q.y()) - square(q.z()));
  rpy(1) = std::asin(as);
  rpy(0) =std::atan2(2 * (q.y() * q.z() + q.w() * q.x()),square(q.w()) - square(q.x()) - square(q.y()) + square(q.z()));
  return rpy;
}

void BDrdsim(float m, float rGx, float rGy, float rGz, float rfo1x, float rfo1y, float rfo1z, Matrix<fpt,12,3>& BDrd){//(36)
  float t2;
  t2 = 1.0F / m;
  BDrd(0) = 0.0F;
  BDrd(1) = 0.0F;
  BDrd(2) = 0.0F;
  BDrd(3) = 0.0F;
  BDrd(4) = 0.0F;
  BDrd(5) = 0.0F;
  BDrd(6) = 0.0F;
  BDrd(7) = -rGz + rfo1z;
  BDrd(8) = rGy - rfo1y;
  BDrd(9) = t2;
  BDrd(10) = 0.0F;
  BDrd(11) = 0.0F;
  BDrd(12) = 0.0F;
  BDrd(13) = 0.0F;
  BDrd(14) = 0.0F;
  BDrd(15) = 0.0F;
  BDrd(16) = 0.0F;
  BDrd(17) = 0.0F;
  BDrd(18) = rGz - rfo1z;
  BDrd(19) = 0.0F;
  BDrd(20) = -rGx + rfo1x;
  BDrd(21) = 0.0F;
  BDrd(22) = t2;
  BDrd(23) = 0.0F;
  BDrd(24) = 0.0F;
  BDrd(25) = 0.0F;
  BDrd(26) = 0.0F;
  BDrd(27) = 0.0F;
  BDrd(28) = 0.0F;
  BDrd(29) = 0.0F;
  BDrd(30) = -rGy + rfo1y;
  BDrd(31) = rGx - rfo1x;
  BDrd(32) = 0.0F;
  BDrd(33) = 0.0F;
  BDrd(34) = 0.0F;
  BDrd(35) = t2;
}

void ACrdsim(float EA1, float EA2, float EA3, float I1, float I2, float I3, float Ld1, float Ld2, float Ld3, float dr1, float dr2, float dr3, float m, float r1, float r2, float r3, float rGx, float rGy, float rGz, Matrix<fpt,12,13>& ACc){
  float t2;
  float t3;
  float t4;
  float t5;
  float t6;
  float t7;
  float t17;
  float t18;
  float t19;
  float t60;
  float t26;
  float t27;
  float t28;
  float t29;
  float t33_tmp;
  float t39;
  float t40;
  float t41;
  float t45;
  float t46;
  float t47;
  float t48;
  float t67;
  float t68;
  float t72;
  float t64;
  float t65;
  float t70;
  float t71;
  float t73;
  float out;
  float b_out;
  float c_out;
  float d_out;
  float e_out;
  float f_out;
  float g_out;
  float h_out;
  float i_out;
  float j_out;
  float k_out;
  float l_out;
  float m_out;
  float n_out;
  float o_out;
  t2 = static_cast<float>(cos(static_cast<double>(EA1)));
  t3 = static_cast<float>(cos(static_cast<double>(EA2)));
  t4 = static_cast<float>(cos(static_cast<double>(EA3)));
  t5 = static_cast<float>(sin(static_cast<double>(EA1)));
  t6 = static_cast<float>(sin(static_cast<double>(EA2)));
  t7 = static_cast<float>(sin(static_cast<double>(EA3)));
  t17 = 1.0F / I1;
  t18 = 1.0F / I2;
  t19 = 1.0F / I3;
  t60 = m * 9.81F;
  t26 = t2 * t4;
  t27 = t2 * t7;
  t28 = t4 * t5;
  t29 = t5 * t7;
  t33_tmp = Ld3 * t3;
  t39 = r1 + -rGx;
  t40 = r2 + -rGy;
  t41 = r3 + -rGz;
  t45 = t6 * t26;
  t46 = t6 * t27;
  t47 = t6 * t28;
  t48 = t6 * t29;
  t67 = t26 + t48;
  t68 = t29 + t45;
  t72 = (Ld2 * t3 * t5 + -(Ld3 * t6)) + Ld1 * t2 * t3;
  t64 = Ld3 + m * (dr1 * r2 + -(dr2 * r1));
  t65 = Ld1 + m * (dr2 * r3 + -(dr3 * r2));
  t70 = t27 + -t47;
  t71 = t28 + -t46;
  t73 = (((Ld1 * t29 + t33_tmp * t4) + -(Ld2 * t27)) + Ld1 * t45) + Ld2 * t47;
  t27 = (((t33_tmp * t7 + Ld2 * t26) + -(Ld1 * t28)) + Ld1 * t46) + Ld2 * t48;
  t29 = Ld2 + -(m * (dr1 * r3 + -(dr3 * r1)));
  t45 = -t18 * t19 * (I2 + -I3);
  t47 = t17 * t19 * (I1 + -I3);
  t26 = -t17 * t18 * (I1 + -I2);
  t28 = dr2 * m;
  t33_tmp = t28 * t3;
  t46 = dr3 * m * t18;
  t48 = -dr3 * m * t19;
  out = dr1 * m;
  b_out = out * t3;
  c_out = m * t3;
  d_out = -m * t18 * t41;
  e_out = c_out * t7 * t18;
  f_out = m * t19;
  g_out = f_out * t41;
  h_out = c_out * t4 * t19;
  i_out = c_out * t5 * t17;
  j_out = m * t18;
  k_out = t2 * t3 * t17;
  l_out = t3 * t5 * t17;
  m_out = t18 * t67;
  n_out = t3 * t7 * t18;
  o_out = t3 * t4 * t19;
    ACc(0) = 0.0F;
  ACc(1) = t45 * t73;
  ACc(2) = t45 * t27;
  ACc(3) = 0.0F;
  ACc(4) = 0.0F;
  ACc(5) = 0.0F;
  ACc(6) = 0.0F;
  ACc(7) = 0.0F;
  ACc(8) = 0.0F;
  ACc(9) = 0.0F;
  ACc(10) = 0.0F;
  ACc(11) = 0.0F;
  ACc(12) = t47 * t73;
  ACc(13) = 0.0F;
  ACc(14) = t47 * t72;
  ACc(15) = 0.0F;
  ACc(16) = 0.0F;
  ACc(17) = 0.0F;
  ACc(18) = 0.0F;
  ACc(19) = 0.0F;
  ACc(20) = 0.0F;
  ACc(21) = 0.0F;
  ACc(22) = 0.0F;
  ACc(23) = 0.0F;
  ACc(24) = t26 * t27;
  ACc(25) = t26 * t72;
  ACc(26) = 0.0F;
  ACc(27) = 0.0F;
  ACc(28) = 0.0F;
  ACc(29) = 0.0F;
  ACc(30) = 0.0F;
  ACc(31) = 0.0F;
  ACc(32) = 0.0F;
  ACc(33) = 0.0F;
  ACc(34) = 0.0F;
  ACc(35) = 0.0F;
  ACc(36) = m * t17 * (dr2 * t6 + dr3 * t3 * t5);
  ACc(37) = t46 * t67 - t33_tmp * t7 * t18;
  ACc(38) = t48 * t70 - t33_tmp * t4 * t19;
  ACc(39) = 0.0F;
  ACc(40) = 0.0F;
  ACc(41) = 0.0F;
  ACc(42) = 0.0F;
  ACc(43) = t60;
  ACc(44) = 0.0F;
  ACc(45) = 0.0F;
  ACc(46) = 0.0F;
  ACc(47) = 0.0F;
  ACc(48) = -m * t17 * (dr1 * t6 + dr3 * t2 * t3);
  ACc(49) = t46 * t71 + b_out * t7 * t18;
  ACc(50) = t48 * t68 + b_out * t4 * t19;
  ACc(51) = 0.0F;
  ACc(52) = 0.0F;
  ACc(53) = 0.0F;
  ACc(54) = -t60;
  ACc(55) = 0.0F;
  ACc(56) = 0.0F;
  ACc(57) = 0.0F;
  ACc(58) = 0.0F;
  ACc(59) = 0.0F;
  ACc(60) = c_out * t17 * (dr2 * t2 - dr1 * t5);
  ACc(61) = -dr1 * m * t18 * t67 - t28 * t18 * t71;
  ACc(62) = t28 * t19 * t68 + out * t19 * t70;
  ACc(63) = 0.0F;
  ACc(64) = 0.0F;
  ACc(65) = 0.0F;
  ACc(66) = 0.0F;
  ACc(67) = 0.0F;
  ACc(68) = 0.0F;
  ACc(69) = 0.0F;
  ACc(70) = 0.0F;
  ACc(71) = 0.0F;
  ACc(72) = k_out;
  ACc(73) = -t18 * t71;
  ACc(74) = t19 * t68;
  ACc(75) = 0.0F;
  ACc(76) = 0.0F;
  ACc(77) = 0.0F;
  ACc(78) = 0.0F;
  ACc(79) = 0.0F;
  ACc(80) = 0.0F;
  ACc(81) = 0.0F;
  ACc(82) = 0.0F;
  ACc(83) = 0.0F;
  ACc(84) = l_out;
  ACc(85) = m_out;
  ACc(86) = -t19 * t70;
  ACc(87) = 0.0F;
  ACc(88) = 0.0F;
  ACc(89) = 0.0F;
  ACc(90) = 0.0F;
  ACc(91) = 0.0F;
  ACc(92) = 0.0F;
  ACc(93) = 0.0F;
  ACc(94) = 0.0F;
  ACc(95) = 0.0F;
  ACc(96) = -t6 * t17;
  ACc(97) = n_out;
  ACc(98) = o_out;
  ACc(99) = 0.0F;
  ACc(100) = 0.0F;
  ACc(101) = 0.0F;
  ACc(102) = 0.0F;
  ACc(103) = 0.0F;
  ACc(104) = 0.0F;
  ACc(105) = 0.0F;
  ACc(106) = 0.0F;
  ACc(107) = 0.0F;
  ACc(108) = -m * t6 * t17 * t40 - i_out * t41;
  ACc(109) = d_out * t67 + e_out * t40;
  ACc(110) = g_out * t70 + h_out * t40;
  ACc(111) = 1.0F;
  ACc(112) = 0.0F;
  ACc(113) = 0.0F;
  ACc(114) = 0.0F;
  ACc(115) = 0.0F;
  ACc(116) = 0.0F;
  ACc(117) = 0.0F;
  ACc(118) = 0.0F;
  ACc(119) = 0.0F;
  ACc(120) = m * t6 * t17 * t39 + m * t2 * t3 * t17 * t41;
  ACc(121) = d_out * t71 - e_out * t39;
  ACc(122) = g_out * t68 - h_out * t39;
  ACc(123) = 0.0F;
  ACc(124) = 1.0F;
  ACc(125) = 0.0F;
  ACc(126) = 0.0F;
  ACc(127) = 0.0F;
  ACc(128) = 0.0F;
  ACc(129) = 0.0F;
  ACc(130) = 0.0F;
  ACc(131) = 0.0F;
  ACc(132) = -m * t2 * t3 * t17 * t40 + i_out * t39;
  ACc(133) = j_out * t39 * t67 + j_out * t40 * t71;
  ACc(134) = -m * t19 * t40 * t68 - f_out * t39 * t70;
  ACc(135) = 0.0F;
  ACc(136) = 0.0F;
  ACc(137) = 1.0F;
  ACc(138) = 0.0F;
  ACc(139) = 0.0F;
  ACc(140) = 0.0F;
  ACc(141) = 0.0F;
  ACc(142) = 0.0F;
  ACc(143) = 0.0F;
  ACc(144) = (t6 * t17 * t64 - k_out * t65) - l_out * t29;
  ACc(145) = (t18 * t65 * t71 - m_out * t29) - n_out * t64;
  ACc(146) = (-t19 * t65 * t68 + t19 * t29 * t70) - o_out * t64;
  ACc(147) = 0.0F;
  ACc(148) = 0.0F;
  ACc(149) = 0.0F;
  ACc(150) = rGy * t60;
  ACc(151) = m * rGx * -9.81F;
  ACc(152) = 0.0F;
  ACc(153) = 0.0F;
  ACc(154) = 0.0F;
  ACc(155) = -9.81F;
}