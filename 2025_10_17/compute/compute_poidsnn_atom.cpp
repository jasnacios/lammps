/* ----------------------------------------------------------------------
   ComputePoidsNNAtom : expose les poids du NN stockés dans atom->poidsnn
------------------------------------------------------------------------- */

#include "compute_poidsnn_atom.h"
#include "atom.h"
#include "update.h"
#include "error.h"
#include "utils.h"

using namespace LAMMPS_NS;

ComputePoidsNNAtom::ComputePoidsNNAtom(LAMMPS *lmp, int narg, char **arg)
  : Compute(lmp, narg, arg)
{
  if (narg < 4)
    error->all(FLERR,"Illegal compute poidsnn/atom command");

  ncols = utils::inumeric(FLERR, arg[3], false, lmp);
  if (ncols <= 0)
    error->all(FLERR,"compute poidsnn/atom: N must be > 0");

  peratom_flag = 1;
  size_peratom_cols = ncols;
}

/* ---------------------------------------------------------------------- */

void ComputePoidsNNAtom::init()
{
  if (!atom->poidsnn)
    error->all(FLERR,"compute poidsnn/atom: atom->poidsnn not defined");
}

/* ---------------------------------------------------------------------- */

void ComputePoidsNNAtom::compute_peratom()
{
  invoked_peratom = update->ntimestep;

  // Zéro-copie : on pointe directement vers atom->poidsnn
  array_atom = atom->poidsnn;
}

/* ---------------------------------------------------------------------- */

double ComputePoidsNNAtom::memory_usage()
{
  return atom->nmax * ncols * sizeof(double);
}
