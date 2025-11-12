/* ----------------------------------------------------------------------
 LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
 http://lammps.sandia.gov, Sandia National Laboratories
 Steve Plimpton, sjplimp@sandia.gov
------------------------------------------------------------------------- */
#include "fix_selfbrownian2D.h"
#include <cmath>
#include "atom.h"
#include "update.h"
#include "error.h"
#include "comm.h"
#include "random_mars.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixSelfBrownian2D::FixSelfBrownian2D(LAMMPS *lmp, int narg, char **arg) :
Fix(lmp, narg, arg)
{
  if (narg < 7) error->all(FLERR,"Illegal fix brownian/2d command");
  
  frtt = utils::numeric(FLERR,arg[3],false,lmp);
  frtr = utils::numeric(FLERR,arg[4],false,lmp);
  Dt = utils::numeric(FLERR,arg[5],false,lmp);
  seed = utils::numeric(FLERR,arg[6],false,lmp);
  
  random = new RanMars(lmp, seed + comm->me);
  
  // Indiquer que ce fix réalise une intégration temporelle 
  time_integrate = 1;
}

/* ---------------------------------------------------------------------- */

FixSelfBrownian2D::~FixSelfBrownian2D()
{
  delete random;
}

/* ---------------------------------------------------------------------- */

void FixSelfBrownian2D::init()
{
  dt = update->dt;
  
}

/* ---------------------------------------------------------------------- */

int FixSelfBrownian2D::setmask()
{
  int mask = 0;
  mask |= INITIAL_INTEGRATE;
  mask |= FINAL_INTEGRATE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixSelfBrownian2D::setup(int vflag)
{
  initial_integrate(vflag);
}

/* ---------------------------------------------------------------------- */

void FixSelfBrownian2D::initial_integrate(int /* vflag */)
{
  double **f = atom->f;
  double **x = atom->x;
  double **v = atom->v;
  double *ztorque = atom->ztorque;
  double *Dr = atom->Dr;

  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  double *ang2D = atom->ang2D;
  int step = update->ntimestep;
  
  if (step <= 1) {
    for (int i = 0; i < nlocal; i++){       
  
      ang2D[i] = random->uniform() * 2.0 * M_PI;
      Dr[i] = 0.02;
  
       }
  }
  // Active force respecting the 2D orientation
  if (step > 1) {
    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) {
        double dx = f[i][0]/frtt*dt + sqrt(2.0*Dt*dt/(frtt*frtt))*random->gaussian();
        double dy = f[i][1]/frtt*dt + sqrt(2.0*Dt*dt/(frtt*frtt))*random->gaussian();
        double dtheta = dt*ztorque[i]/frtr + sqrt(2.0*Dr[i]*dt/(frtr*frtr))*random->gaussian();
      
        // Appliquer la force
        x[i][0] += dx;
        x[i][1] += dy;
        v[i][0] = dx/dt;
        v[i][1] = dy/dt;
        
        // Appliquer le couple
        ang2D[i] += dtheta;
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixSelfBrownian2D::final_integrate()
{
  // Rien à faire ici, mais nécessaire pour l'intégration temporelle
}

/* ---------------------------------------------------------------------- */

