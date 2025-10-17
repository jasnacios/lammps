/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov
------------------------------------------------------------------------- */

#include "fix_controler.h"
#include <cmath>
#include "atom.h"
#include "update.h"
#include "error.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixControler::FixControler(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg < 4) error->all(FLERR,"Illegal fix /active/force command");

  nc = utils::numeric(FLERR,arg[3],false,lmp);

}

/* ---------------------------------------------------------------------- */

FixControler::~FixControler()
{
} 
/* ---------------------------------------------------------------------- */

int FixControler::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixControler::setup(int vflag)
{
  post_force(vflag);
} 

/* ---------------------------------------------------------------------- */

void FixControler::post_force(int vflag)
{
  double *lightintensity = atom->lightintensity;
  double *Fa = atom->Fa;
  double *zeta = atom->zeta;
  double **poidsnn = atom->poidsnn;


  int *mask = atom->mask;    
  int nlocal = atom->nlocal;

  
  for (int i = 0; i < nlocal; i++) {
    //Fa[i] = poidsnn[i][0]; //poidsnn[i][1];

    if (mask[i] & groupbit) {
        if (lightintensity[i] >= poidsnn[i][0]) {
          Fa[i] = 0.1;
        } else {
          Fa[i] = 1.0;
        }
   }
  }
}
