/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov
------------------------------------------------------------------------- */

#ifdef FIX_CLASS

FixStyle(robot_learning,FixRobotLearning)

#else

#ifndef LMP_FIX_ROBOT_LEARNING_H
#define LMP_FIX_ROBOT_LEARNING_H

#include "fix.h"

namespace LAMMPS_NS {

class FixRobotLearning : public Fix {
 public:
  FixRobotLearning(class LAMMPS *, int, char **);
  ~FixRobotLearning();
  int setmask();
  void setup(int);
  void init();
  void post_force(int);

 private:
 class RanMars *random;
 int seed;
 int Nn;
 double alphaq;
 double Dp;
 char* idregion0;
 class Region *region0;
 char* idregion1;
 class Region *region1;
 double alpha;
 double comm_radius2;
 double dt;

};

}

#endif
#endif
