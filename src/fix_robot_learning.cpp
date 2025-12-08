/* ----------------------------------------------------------------------
 LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
 http://lammps.sandia.gov, Sandia National Laboratories
 Steve Plimpton, sjplimp@sandia.gov
------------------------------------------------------------------------- */
#include "fix_robot_learning.h"
#include <cmath>
#include "atom.h"
#include "update.h"
#include "error.h"
#include "comm.h"
#include "random_mars.h"
#include "region.h"
#include "domain.h"
#include "neigh_list.h"
#include "neighbor.h"
#include "memory.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixRobotLearning::FixRobotLearning(LAMMPS *lmp, int narg, char **arg) :
Fix(lmp, narg, arg)
{
  if (narg < 9) error->all(FLERR,"Illegal fix brownian/2d command");
  
 
  alphaq = utils::numeric(FLERR,arg[3],false,lmp);
  Dp = utils::numeric(FLERR,arg[4],false,lmp);
  idregion0 = utils::strdup(arg[5]);
  region0 = domain->get_region_by_id(arg[5]);
  idregion1 = utils::strdup(arg[6]);
  region1 = domain->get_region_by_id(arg[6]);
  comm_radius2 = utils::numeric(FLERR,arg[7],false,lmp);
  alpha = utils::numeric(FLERR,arg[8],false,lmp);
  Nn = utils::numeric(FLERR,arg[9],false,lmp);
  seed = utils::numeric(FLERR,arg[10],false,lmp);
  
  random = new RanMars(lmp, seed + comm->me);

}

/* ---------------------------------------------------------------------- */

FixRobotLearning::~FixRobotLearning()
{
  delete random;
}

/* ---------------------------------------------------------------------- */

void FixRobotLearning::init()
{
  dt = update->dt;
  
}

/* ---------------------------------------------------------------------- */

int FixRobotLearning::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixRobotLearning::setup(int vflag)
{
  post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixRobotLearning::post_force(int vflag)
{
  int i, j, ii, jj, inum, jnum;
  double xtmp, ytmp, ztmp, delx, dely, delz;
  double rsq;
  int *ilist, *jlist, *numneigh, **firstneigh;

  double **poidsnn = atom->poidsnn;
  double *ztorque = atom->ztorque;
  double *qreward = atom->qreward;
  double *lightintensity = atom->lightintensity;

  double **x = atom->x;
  int step = update->ntimestep;

  NeighList *list = neighbor->lists[0];
  
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

 
  
  if (step < 1) {
    for (int i = 0; i < nlocal; i++) {
    qreward[i] = 0.0;
    // ztorque[i] = 0.0;
    for (int k = 0; k < Nn; k++){
      poidsnn[i][k] = random->uniform();}}
  }
  if (step >= 1) {
    for (ii = 0; ii < inum; ii++) {
      i = ilist[ii];
      xtmp = x[i][0];
      ytmp = x[i][1];
      ztmp = x[i][2];
        
      jlist = firstneigh[i];
      jnum = numneigh[i];

      for (jj = 0; jj < jnum; jj++) {
        j = jlist[jj];
        j &= NEIGHMASK;
    
        delx = xtmp - x[j][0];
        dely = ytmp - x[j][1];
        delz = ztmp - x[j][2];
        rsq = delx * delx + dely * dely + delz * delz;
          
        if(rsq == 0) continue;
        if (rsq < comm_radius2) {
          if (qreward[i] > qreward[j] ) {
              qreward[j] += alpha*(qreward[i] - qreward[j])*dt;
              for (int k = 0; k<Nn; k++) {
                poidsnn[j][k] += alpha*(poidsnn[i][k] - poidsnn[j][k])*dt;
              }
            }
          if (qreward[i] < qreward[j]) {
              qreward[i] += alpha*(qreward[j] - qreward[i])*dt;
              for (int k = 0; k<Nn; k++) {
                poidsnn[i][k] += alpha*(poidsnn[j][k] - poidsnn[i][k])*dt;
              }
            }
          }
        }
      }
    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) {
        for(int k=0;k<Nn;k++) {
          poidsnn[i][k] += random->gaussian() * sqrt(2*dt*Dp);
        }
        for (int k=0;k<Nn;k++) {
          if(poidsnn[i][k] > 1.0) poidsnn[i][k] = 2 - poidsnn[i][k];
          if(poidsnn[i][k] < 0.0) poidsnn[i][k] = - poidsnn[i][k];
        }
      }
    }
    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) {
        if(region0->match(x[i][0], x[i][1], x[i][2])|| region1->match(x[i][0], x[i][1], x[i][2])) {
          lightintensity[i] = 0.9;
        } else {
          lightintensity[i] = 0.1;
        } 
        // Update reward
        qreward[i] += alphaq*(lightintensity[i] - qreward[i]) *dt;
      }
    }  
  }
}

