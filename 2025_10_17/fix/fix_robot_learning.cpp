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
  comm_radius = utils::numeric(FLERR,arg[7],false,lmp);
  alpha = utils::numeric(FLERR,arg[8],false,lmp);

  seed = utils::numeric(FLERR,arg[9],false,lmp);
  
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
    for (int i = 0; i < nlocal; i++){   

  
      qreward[i] = 0.0;
      // poidsnn[i][0] = 0.0;

      //if (i%2 == 0) {
        //poidsnn[i][0] = 1.0 ;
      //} else {
        //poidsnn[i][0] = 0.0;
     //}

      poidsnn[i][0] = random->uniform();
      poidsnn[i][1] = random->uniform();
      poidsnn[i][2] = random->uniform(); 
      poidsnn[i][3] = random->uniform();
      poidsnn[i][4] = random->uniform();
      poidsnn[i][5] = random->uniform();
      poidsnn[i][6] = random->uniform();
      poidsnn[i][7] = random->uniform();

       }
  }
  if (step >= 1) {
    
      // loop over neighbors of my atoms
    
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
        if (rsq < comm_radius) {
          if (qreward[i] > qreward[j]) {
              qreward[j] += alpha*(qreward[i] - qreward[j])*dt;
              poidsnn[j][0] += alpha*(poidsnn[i][0] - poidsnn[j][0])*dt;
              poidsnn[j][1] += alpha*(poidsnn[i][1] - poidsnn[j][1])*dt;
              poidsnn[j][2] += alpha*(poidsnn[i][2] - poidsnn[j][2])*dt;  
              poidsnn[j][3] += alpha*(poidsnn[i][3] - poidsnn[j][3])*dt;
              poidsnn[j][4] += alpha*(poidsnn[i][4] - poidsnn[j][4])*dt;
              poidsnn[j][5] += alpha*(poidsnn[i][5] - poidsnn[j][5])*dt;
              poidsnn[j][6] += alpha*(poidsnn[i][6] - poidsnn[j][6])*dt;
              poidsnn[j][7] += alpha*(poidsnn[i][7] - poidsnn[j][7])*dt;
            }
          if (qreward[i] < qreward[j]) {
              qreward[i] += alpha*(qreward[j] - qreward[i])*dt;
              poidsnn[i][0] += alpha*(poidsnn[j][0] - poidsnn[i][0])*dt;
              poidsnn[i][1] += alpha*(poidsnn[j][1] - poidsnn[i][1])*dt;
              poidsnn[i][2] += alpha*(poidsnn[j][2] - poidsnn[i][2])*dt;  
              poidsnn[i][3] += alpha*(poidsnn[j][3] - poidsnn[i][3])*dt;
              poidsnn[i][4] += alpha*(poidsnn[j][4] - poidsnn[i][4])*dt;
              poidsnn[i][5] += alpha*(poidsnn[j][5] - poidsnn[i][5])*dt;
              poidsnn[i][6] += alpha*(poidsnn[j][6] - poidsnn[i][6])*dt;
              poidsnn[i][7] += alpha*(poidsnn[j][7] - poidsnn[i][7])*dt;
            }
           
          }
        }
      }
    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) {
        
        poidsnn[i][0] += random->gaussian() * sqrt(2*dt*Dp); 
        poidsnn[i][1] += random->gaussian() * sqrt(2*dt*Dp); 
        poidsnn[i][2] += random->gaussian() * sqrt(2*dt*Dp);
        poidsnn[i][3] += random->gaussian() * sqrt(2*dt*Dp);
        poidsnn[i][4] += random->gaussian() * sqrt(2*dt*Dp);
        poidsnn[i][5] += random->gaussian() * sqrt(2*dt*Dp);
        poidsnn[i][6] += random->gaussian() * sqrt(2*dt*Dp);
        poidsnn[i][7] += random->gaussian() * sqrt(2*dt*Dp);

        // Clamp values to 0 and 1
        if(poidsnn[i][0] > 1.0) poidsnn[i][0] = 2 - poidsnn[i][0];
        if(poidsnn[i][0] < 0.0) poidsnn[i][0] = - poidsnn[i][0];
      
        if(poidsnn[i][1] > 1.0) poidsnn[i][1] = 2 - poidsnn[i][1];
        if(poidsnn[i][1] < 0.0) poidsnn[i][1] = - poidsnn[i][1];
      
        if(poidsnn[i][2] > 1.0) poidsnn[i][2] = 2 - poidsnn[i][2];
        if(poidsnn[i][2] < 0.0) poidsnn[i][2] = - poidsnn[i][2];

        if(poidsnn[i][3] > 1.0) poidsnn[i][3] = 2 - poidsnn[i][3];
        if(poidsnn[i][3] < 0.0) poidsnn[i][3] = - poidsnn[i][3];

        if(poidsnn[i][4] > 1.0) poidsnn[i][4] = 2 - poidsnn[i][4];
        if(poidsnn[i][4] < 0.0) poidsnn[i][4] = - poidsnn[i][4];

        if(poidsnn[i][5] > 1.0) poidsnn[i][5] = 2 - poidsnn[i][5];
        if(poidsnn[i][5] < 0.0) poidsnn[i][5] = - poidsnn[i][5];
        
        if(poidsnn[i][6] > 1.0) poidsnn[i][6] = 2 - poidsnn[i][6]; 
        if(poidsnn[i][6] < 0.0) poidsnn[i][6] = - poidsnn[i][6];

        if(poidsnn[i][7] > 1.0) poidsnn[i][7] = 2 - poidsnn[i][7];
        if(poidsnn[i][7] < 0.0) poidsnn[i][7] = - poidsnn[i][7];  
      }
    }
    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) {
        if(region0->match(x[i][0], x[i][1], x[i][2])|| region1->match(x[i][0], x[i][1], x[i][2])) {
          lightintensity[i] = 0.66;
        } else {
          lightintensity[i] = 0.33;
        } 
        // Update reward
        qreward[i] += alphaq*(lightintensity[i] - qreward[i]) *dt;
      }
    }  
  }
}

