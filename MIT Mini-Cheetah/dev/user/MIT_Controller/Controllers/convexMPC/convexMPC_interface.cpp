#include "SolverMPC.h"
#include <pthread.h>
#include <string.h>

#define K_NUM_LEGS 4

problem_setup problem_configuration;
pthread_mutex_t problem_cfg_mt;
pthread_mutex_t update_mt;
update_data_t update;
pthread_t solve_thread;

u8 first_run = 1;

void initialize_mpc(){
  //printf("Initializing MPC!\n");
  if(pthread_mutex_init(&problem_cfg_mt,NULL)!=0)
    printf("[MPC ERROR] Failed to initialize problem configuration mutex.\n");

  if(pthread_mutex_init(&update_mt,NULL)!=0)
    printf("[MPC ERROR] Failed to initialize update data mutex.\n");

#ifdef K_DEBUG
  printf("[MPC] Debugging enabled.\n");
    printf("[MPC] Size of problem setup struct: %ld bytes.\n", sizeof(problem_setup));
    printf("      Size of problem update struct: %ld bytes.\n",sizeof(update_data_t));
    printf("      Size of MATLAB floating point type: %ld bytes.\n",sizeof(mfp));
    printf("      Size of flt: %ld bytes.\n",sizeof(flt));
#else
  //printf("[MPC] Debugging disabled.\n");
#endif
}

void setup_problem(double dt, int horizon, float* mu, double f_max, int sga, Eigen::Matrix<float,3,3>* Rpla, double savdatatxt){
  //mu = 0.6;
  if(first_run){
    first_run = false;
    initialize_mpc();
  }
#ifdef K_DEBUG
  printf("[MPC] Got new problem configuration!\n");
    printf("[MPC] Prediction horizon length: %d\n      Force limit: %.3f, friction %.3f\n      dt: %.3f\n",
            horizon,f_max,mu,dt);
#endif
  //pthread_mutex_lock(&problem_cfg_mt);
  problem_configuration.savdatatxt = savdatatxt;
  problem_configuration.horizon = horizon;
  problem_configuration.f_max = f_max;
  //problem_configuration.mu = mu;
  memcpy((void*)problem_configuration.mu,(void*)mu,sizeof(float)*4);
  problem_configuration.dt = dt;
  for(s16 i=0;i<4;i++){
  problem_configuration.Rpla[i] = Rpla[i];}
  if (sga ==-1){sga = 4*horizon;}
  //pthread_mutex_unlock(&problem_cfg_mt);
  resize_qp_mats(horizon, sga);
}

//inline to motivate gcc to unroll the loop in here.
inline void mfp_to_flt(flt* dst, mfp* src, s32 n_items){
  for(s32 i = 0; i < n_items; i++)
    *dst++ = *src++;
}

inline void mint_to_u8(u8* dst, mint* src, s32 n_items){
  for(s32 i = 0; i < n_items; i++)
    *dst++ = *src++;
}

int has_solved = 0;

//void *call_solve(void* ptr){
//  solve_mpc(&update, &problem_configuration);
//}
//safely copies problem data and starts the solve

void update_problem_data_floats(float* p, float* v, float* q, float* w,
                                float* r, double* weights,
                                float* state_trajectory, double* uweights, int* gait, float* seLcomp){
  mint_to_u8(update.gait,gait,4*problem_configuration.horizon);
  memcpy((void*)update.seLcomp,(void*)seLcomp,sizeof(float)*3);
  memcpy((void*)update.p,(void*)p,sizeof(float)*3);
  memcpy((void*)update.v,(void*)v,sizeof(float)*3);
  memcpy((void*)update.q,(void*)q,sizeof(float)*4);
  memcpy((void*)update.w,(void*)w,sizeof(float)*3);
  memcpy((void*)update.r,(void*)r,sizeof(float)*12);
  memcpy((void*)update.weights,(void*)weights,sizeof(double)*12);
  memcpy((void*)update.uweights,(void*)uweights,sizeof(double)*3);
  memcpy((void*)update.traj,(void*)state_trajectory, sizeof(float) * 13 * problem_configuration.horizon);
  solve_mpc(&update, &problem_configuration);
  has_solved = 1;
}

void update_traj(float *trajSome) {
  for(s16 itra = 0; itra < (13*(problem_configuration.horizon)); itra++)
    trajSome[itra]=update.traj[itra];
}

void update_trajh(float *trajSome, int horiupt) {
  for(s16 itra = 0; itra < (13*horiupt); itra++)
    trajSome[itra]=update.traj[itra];
}

double get_solution(int index){
  if(!has_solved){
    printf("\n[Warning] Not solved yet!\n");
    return 0.f;
  }
  mfp* qs = get_q_soln();
  if (update.gait[index/3]==0){return 0.f;}
  else{
    int sumindc=0;
    for (int indc=0; indc<(index/3)+1; indc++){
      sumindc+=update.gait[indc];
    }
    return qs[3*(sumindc-1)+(index % 3)];
  }
}
/*
Matrix<fpt,Dynamic,1> get_Xsolin_int()
{
    Matrix<fpt,Dynamic,1> Xsolin_int;
    Xsolin_int=get_Xsolin();
    return Xsolin_int;
}
*/