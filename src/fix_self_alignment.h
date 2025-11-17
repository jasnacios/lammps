/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov
------------------------------------------------------------------------- */

#ifdef FIX_CLASS

FixStyle(self_alignment,FixSelfAlignment)

#else

#ifndef LMP_FIX__SELF_ALIGNMENT_H
#define LMP_FIX__SELF_ALIGNMENT_H

#include "fix.h"

namespace LAMMPS_NS {

class FixSelfAlignment : public Fix {
 public:
  FixSelfAlignment(class LAMMPS *, int, char **);
  ~FixSelfAlignment();
  int setmask();
  void setup(int);
  void post_force(int);

 private:
};

}

#endif
#endif
