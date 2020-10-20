#ifndef PROJECT_GAIT_H
#define PROJECT_GAIT_H

#include <string>
#include <queue>

#include "cppTypes.h"

template <typename T>
using MatDx1 = typename Eigen::Matrix<T, Eigen::Dynamic, 1>;

using Eigen::Matrix;
using Eigen::Dynamic;
using Eigen::Array;
using Eigen::Array4f;
using Eigen::Array4i;

class Gait {
public:
  virtual ~Gait()= default;//
  
  virtual Vec4<float> getContactState() = 0;
  virtual Vec4<float> getSwingState() = 0;
  virtual Vec4<float> getis_slidcntct() = 0;
  virtual int* getMpcTable() = 0;
  virtual void setIterations(int iterationsBetweenMPC, int currentIteration) = 0;
  virtual void setItCorrected(int iterationsBetweenMPC, int currentIteration) = 0;
  virtual void setItCorrectedStart(int iterationsBetweenMPC, int currentIteration) = 0;
  virtual float getCurrentStanceTime(float dtMPC, int leg) = 0;
  virtual float getCurrentSwingTime(float dtMPC, int leg) = 0;
  virtual int getCurrentGaitPhase() = 0;
  virtual int getCurrentHybridMode() = 0;
  virtual int getItNxtHybMode() = 0;
  virtual void debugPrint() {};
  virtual int getHoriz() = 0;
  virtual bool getFlightState() = 0;
  virtual void setdurations(int dur){(void) dur;};
  virtual void set_nIterations(int _nIt){(void) _nIt;};
  virtual void recomputeVals(int jumpConf, int _Np1, int _Np2, int _Np3) = 0;
  virtual float get_phaseLO() = 0;
  virtual float get_phasela() = 0;
  virtual int getFliMacNH(int NumFlightM) = 0;
  virtual void setIterationssb(int _iteration1, float _phase1){(void) _iteration1;(void) _phase1;};
  int gaitphascorr = 0; 
  int countflystage = 0; 
  int jumpConfGait = 0;
  bool ischained = false;

protected:
  std::string _name;
};

class ParkourGait : public Gait {
public:
  ParkourGait(int jumpConf, int _Np1, int _Np2, int _Np3, const std::string &name);
  ~ParkourGait();
  Vec4<float> getContactState();
  Vec4<float> getSwingState();
  Vec4<float> getis_slidcntct();
  int* getMpcTable();
  void getMpcTablebase();
  void setIterations(int iterationsBetweenMPC, int currentIteration);
  void setItCorrectedStart(int iterationsBetweenMPC, int currentIteration);
  void setIterationssb(int _iteration1, float _phase1);
  void setvNHacH(int jumpConf, int _Np1, int _Np2, int _Np3);//void getvNHacH();
  void set_nIterations(int _nIt);
  void setdurations(int dur);
  float getCurrentStanceTime(float dtMPC, int leg);
  float getCurrentSwingTime(float dtMPC, int leg);
  int getCurrentGaitPhase();
  int getCurrentHybridMode();
  int getItNxtHybMode();
  void setItCorrected(int iterationsBetweenMPC, int currentIteration);
  void debugPrint();
  int getHoriz();
  void recomputeVals(int jumpConf, int _Np1, int _Np2, int _Np3);
  bool getFlightState();
  float get_phaseLO();
  float get_phasela();
  int getFliMacNH(int NumFlightM);
  MatDx1<int> FliMacNH;

private:
  MatDx1<int> vNH;
  Array<MatDx1<int>,Dynamic,1> acH;
  int NModes;
  int* _mpc_table;
  int* _mpc_table_base;
  int _iteration;
  int _nIterations;
  int _CurrHybridMode;
  float _phase;
};

class OffsetDurationGait : public Gait {
public:
  OffsetDurationGait(int nSegment, Vec4<int> offset, Vec4<int> durations, const std::string& name);
  ~OffsetDurationGait();
  Vec4<float> getContactState();
  Vec4<float> getSwingState();
  Vec4<float> getis_slidcntct();
  int* getMpcTable();
  void setIterations(int iterationsBetweenMPC, int currentIteration);
  void setItCorrectedStart(int iterationsBetweenMPC, int currentIteration);
  void setIterationssb(int _iteration1, float _phase1);
  void setdurations(int dur);
  void set_nIterations(int _nIt);
  float getCurrentStanceTime(float dtMPC, int leg);
  float getCurrentSwingTime(float dtMPC, int leg);
  int getCurrentGaitPhase();
  int getCurrentHybridMode();
  int getItNxtHybMode();
  void setItCorrected(int iterationsBetweenMPC, int currentIteration);
  void debugPrint();
  int getHoriz();
  void recomputeVals(int jumpConf, int _Np1, int _Np2, int _Np3);
  bool getFlightState();
  float get_phaseLO();
  float get_phasela();
  int getFliMacNH(int NumFlightM);

private:
  int* _mpc_table;
  Array4i _offsets; // offset in mpc segments
  Array4i _durations; // duration of step in mpc segments
  Array4f _offsetsFloat; // offsets in phase (0 to 1)
  Array4f _durationsFloat; // durations in phase (0 to 1)
  int _stance;
  int _swing;
  int _iteration;
  int _nIterations;
  float _phase;
};

#endif //PROJECT_GAIT_H