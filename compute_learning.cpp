/* ----------------------------------------------------------------------
   ComputePoidsNNAtom : expose les poids du NN stockés dans atom->poidsnn
------------------------------------------------------------------------- */

#include "compute_learning.h"
#include "atom.h"
#include "update.h"
#include "error.h"
#include "utils.h"
#include "memory.h"

using namespace LAMMPS_NS;

ComputeLearning::ComputeLearning(LAMMPS *lmp, int narg, char **arg)
  : Compute(lmp, narg, arg)
{
  if (narg < 4)
    error->all(FLERR,"Illegal compute learning command");

  ncols = utils::inumeric(FLERR, arg[3], false, lmp);
  if (ncols <= 0)
    error->all(FLERR,"compute learning: N must be > 0");

  peratom_flag = 1;
  size_peratom_cols = ncols+1;

  nmax_alloc = 0;
  array_atom = nullptr;
}

/* ---------------------------------------------------------------------- */

void ComputeLearning::init()
{
  if (!atom->poidsnn)
    error->all(FLERR,"compute learning: atom->poidsnn not defined");
}

/* ---------------------------------------------------------------------- */

void ComputeLearning::compute_peratom()
{
  invoked_peratom = update ->ntimestep;
  if (atom->nmax > nmax_alloc) {
    nmax_alloc = atom->nmax;
    if (array_atom) memory -> destroy(array_atom);
    memory->create(array_atom, atom->nmax, size_peratom_cols, "learning:array");
  }

  double *q = atom-> qreward;
  double **w = atom-> poidsnn;
  int nlocal = atom-> nlocal;

  for (int i = 0; i<nlocal; i++) {
    array_atom[i][0] = q[i];
    for (int j=0; j<ncols; j++){
      array_atom[i][j+1] = w[i][j];
    }
  }

  // Zéro-copie : on pointe directement vers atom->poidsnn
  // array_atom = atom->poidsnn;*/
}

/* ---------------------------------------------------------------------- */

double ComputeLearning::memory_usage()
{
  return atom->nmax * ncols * sizeof(double);
}