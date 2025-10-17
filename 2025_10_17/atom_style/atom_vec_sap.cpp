/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "atom_vec_sap.h"

#include "atom.h"
#include "domain.h"
#include "memory.h"

#include <cmath>
 
using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

AtomVecSap::AtomVecSap(LAMMPS *lmp) : AtomVec(lmp)
{
  molecular = Atom::ATOMIC;
  mass_type = PER_TYPE;

  atom->Fa_flag = 1;
  atom->Dr_flag = 1;
  atom->zeta_flag = 1;
  

  // strings with peratom variables to include in each AtomVec method
  // strings cannot contain fields in corresponding AtomVec default strings
  // order of fields in a string does not matter
  // except: fields_data_atom & fields_data_vel must match data file

  fields_grow = { "Fa", "Dr", "zeta"};
  fields_copy = { "Fa", "Dr", "zeta"};
  fields_border = { "Fa", "Dr","zeta"};
  fields_exchange = { "Fa", "Dr", "zeta"};
  fields_restart = { "Fa", "Dr", "zeta"};
  fields_create = { "Fa", "Dr", "zeta"};
  fields_data_atom = {"id", "type", "x"};
  fields_data_vel = {"id", "v"};

  setup_fields();
}

/* ----------------------------------------------------------------------
   set local copies of all grow ptrs used by this class, except defaults
   needed in replicate when 2 atom classes exist and it calls pack_restart()
------------------------------------------------------------------------- */

void AtomVecSap::grow_pointers()
{
  Fa = atom->Fa;
  Dr = atom->Dr;
  zeta = atom->zeta;
}

/* ----------------------------------------------------------------------
   modify what AtomVec::data_atom() just unpacked
   or initialize other atom quantities
------------------------------------------------------------------------- */

void AtomVecSap::data_atom_post(int ilocal)
{
    
}

