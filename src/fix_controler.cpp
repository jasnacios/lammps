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
#include <span>

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixControler::FixControler(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg < 4) error->all(FLERR,"Illegal fix /active/force command");

  nc = utils::numeric(FLERR,arg[3],false,lmp);
  Nn = utils::numeric(FLERR,arg[4],false,lmp);
  h = utils::numeric(FLERR,arg[5],false,lmp);


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
  double *ztorque = atom->ztorque;
  double *zeta = atom->zeta;
  double **poidsnn = atom->poidsnn;



  int *mask = atom->mask;    
  int nlocal = atom->nlocal;

  
  for (int i = 0; i < nlocal; i++) {
    double I = lightintensity[i];
    double *w = poidsnn[i];

    double slope;
    double sigm;

    double slopeg;
    double sigmg;
    double sloped;
    double sigmd;


    if (nc ==1) {
    if (mask[i] & groupbit) 
        if (I >= w[0]) {
          Fa[i] = 0.1;
        } else {
          Fa[i] = 1.0;
        }
      }
    if (nc==2){
      if (mask[i] & groupbit){
        if (w[1]> w[2]){std::swap(w[1], w[2]); std::swap(w[0], w[3]);}
        if (w[1] == w[2]){Fa[i] = 0.5*(w[0]+w[3]);}

        if (I<= w[1]){
          Fa[i] = w[0];
        }
        else if (I >= w[2]) {
          Fa[i] = w[3];
        }
        else {
          Fa[i] = w[0] + (w[3]-w[0])*(I-w[1])/(w[2]-w[1]);
        }
      }
    }
    if (nc==3){
      if (mask[i]& groupbit){
        slope = f_slope(w[1]);
        sigm = sigmoid(slope*(I - w[2]));
        Fa[i] = w[0] + (w[3] - w[0])*sigm;
      }
    }
    if (nc==4){
      if (mask[i]& groupbit){
        slopeg = f_slope(w[1]);
        sigmg = sigmoid(slopeg*(I - w[2]));
        double Fg = w[0] + (w[3] - w[0])*sigmg;
        sloped = f_slope(w[5]);
        sigmd = sigmoid(sloped*(I - w[6]));
        double Fd =  w[4] + (w[7] - w[4])*sigmd;
        double PWM = (Fg+Fd)/2;
        Fa[i] = th_trun(PWM);
        ztorque[i] = 0.5*(Fd - Fg);
      } 
    }
    if (nc ==5){
      int index = 0;
      double layer_1[h];
      double layer_2[h];
      double layer_out[2];
      if (mask[i]& groupbit){
        for (int l = 0; l<h; l++){
          double Wl1 = w[index++]; // les poids paires correspondent aux w
          double b1l = w[index++]; // Les poids impaires correspondent aux seuil 
          double x1l =  Wl1*I +b1l; // on récupère les valeurs des premiers neurones en appliquant W*I  +b (W dans M_{h,1}, b dans R^{h))
          layer_1[l] = Leaky_Relu(x1l); // on applique la fonction d'activation
        }

        for (int l = 0; l<h; l++){
          double x2l = 0;
          for (int c=0; c<h; c++){
            double W2lc = w[index++];
            x2l += W2lc*layer_1[c]; 
          }
          x2l += w[index++];
          layer_2[l] = Leaky_Relu(x2l);
        }
        for (int l=0; l<2; l++){
          double x_out_l = 0;
          for (int c=0; c<h;c++){
            double W_out_lc = w[index++];
            x_out_l += W_out_lc*layer_2[c];
          }
          x_out_l += w[index++];
          if (l==0){
            layer_out[l] = sigmoid(x_out_l);
          }
          if (l==1){
            layer_out[l] = std::tanh(x_out_l);
          }
          
        }
        double f_d= layer_out[0];
        double f_g = layer_out[1];

        Fa[i] = (f_d+f_g)/2;
        ztorque[i] = (f_d - f_g)/2;
        
      }
    }
  }
}
