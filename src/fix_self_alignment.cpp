/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov
------------------------------------------------------------------------- */

#include "fix_self_alignment.h"
#include <cmath>
#include "atom.h"
#include "update.h"
#include "error.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixSelfAlignment::FixSelfAlignment(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg < 3) error->all(FLERR,"Illegal fix /active/force command");

}

/* ---------------------------------------------------------------------- */

FixSelfAlignment::~FixSelfAlignment()
{
} 
/* ---------------------------------------------------------------------- */

int FixSelfAlignment::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixSelfAlignment::setup(int vflag)
{
  post_force(vflag);
} 

/* ---------------------------------------------------------------------- */

void FixSelfAlignment::post_force(int vflag)
{
  double **f = atom->f;
  int *mask = atom->mask;    
  int nlocal = atom->nlocal;
  double *ang2D = atom->ang2D;

  double *ztorque = atom->ztorque;      
  double **v = atom->v;

  double *Fa = atom->Fa;
  double *zeta = atom->zeta;
  int step = update->ntimestep;

  if (step <= 1) {
    for (int i = 0; i < nlocal; i++){       
        zeta[i] = 0.0;
       }
  }

  if (step > 1) {
    for (int i = 0; i < nlocal; i++) {
        if (mask[i] & groupbit) {
           zeta[i] = 0.0;
          double fx = Fa[i] * cos(ang2D[i]);
          double fy = Fa[i] * sin(ang2D[i]);
          double Cz = zeta[i] *  (cos(ang2D[i]) * v[i][1] - sin(ang2D[i]) * v[i][0]);
          
          // Appliquer la force
          f[i][0] += fx;
          f[i][1] += fy;
          ztorque[i] += Cz;
        }
    }
  }
}
