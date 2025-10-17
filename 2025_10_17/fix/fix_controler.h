/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov
------------------------------------------------------------------------- */

#ifdef FIX_CLASS

FixStyle(controler,FixControler)

#else

#ifndef LMP_FIX_CONTROLER_H
#define LMP_FIX_CONTROLER_H

#include "fix.h"

namespace LAMMPS_NS {

class FixControler : public Fix {
 public:
  FixControler(class LAMMPS *, int, char **);
  ~FixControler();
  int setmask();
  void setup(int);
  void post_force(int);

 private:
 double nc;
};

}

#endif
#endif
