/* ----------------------------------------------------------------------
   ComputePoidsNNAtom : expose les poids du NN stockés dans atom->poidsnn
------------------------------------------------------------------------- */

#ifdef COMPUTE_CLASS
ComputeStyle(poidsnn/atom,ComputeLearning)
#else

#ifndef LMP_COMPUTE_learning_H
#define LMP_COMPUTE_learning_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputeLearning : public Compute {
 public:
  ComputeLearning(class LAMMPS *, int, char **);
  ~ComputeLearning() override = default;

  void init() override;
  void compute_peratom() override;
  double memory_usage() override;

 private:
  int ncols;  // nombre de paramètres NN par atome
  int nmax_alloc; // nombre max d'atomes
};

} // namespace LAMMPS_NS

#endif
#endif
