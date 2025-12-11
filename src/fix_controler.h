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
 int nc;
 int Nn;
 int h;
 inline double sigmoid(double x) {
  return 1.0 / (1.0 + std::exp(-x));
}
 inline double f_slope(double S){ return std::pow(4.0*S, 4.0); }
 inline double th_trun(double x){return 0.5*(1 + std::tanh(11.5*(x-0.39)));}
 inline double Leaky_Relu(double x){
   if (x<0){ return 0.1*x;}
   else {return x;}
 }
};

}

#endif
#endif
