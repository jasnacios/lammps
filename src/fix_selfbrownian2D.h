/* ----------------------------------------------------------------------
 LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
 http://lammps.sandia.gov, Sandia National Laboratories
 Steve Plimpton, sjplimp@sandia.gov
------------------------------------------------------------------------- */
#ifdef FIX_CLASS
// clang-format off
FixStyle(SelfBrownian2D,FixSelfBrownian2D);
// clang-format on
#else

#ifndef LMP_FIX_SELF_BROWNIAN2D_H
#define LMP_FIX_SELF_BROWNIAN2D_H

#include "fix.h"

namespace LAMMPS_NS {

class FixSelfBrownian2D : public Fix {
 public:
  FixSelfBrownian2D(class LAMMPS *, int, char **);
  ~FixSelfBrownian2D();
  int setmask();
  void init();
  void setup(int);
  void initial_integrate(int);
  void final_integrate();
  

 private:
  double frtt, frtr, Dt;
  int seed;
  double dt;
  class RanMars *random;
};

}

#endif
#endif
