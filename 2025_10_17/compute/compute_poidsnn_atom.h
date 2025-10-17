/* ----------------------------------------------------------------------
   ComputePoidsNNAtom : expose les poids du NN stockés dans atom->poidsnn
------------------------------------------------------------------------- */

#ifdef COMPUTE_CLASS
ComputeStyle(poidsnn/atom,ComputePoidsNNAtom)
#else

#ifndef LMP_COMPUTE_POIDSNN_ATOM_H
#define LMP_COMPUTE_POIDSNN_ATOM_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputePoidsNNAtom : public Compute {
 public:
  ComputePoidsNNAtom(class LAMMPS *, int, char **);
  ~ComputePoidsNNAtom() override = default;

  void init() override;
  void compute_peratom() override;
  double memory_usage() override;

 private:
  int ncols;  // nombre de paramètres NN par atome
};

} // namespace LAMMPS_NS

#endif
#endif
