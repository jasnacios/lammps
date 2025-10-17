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

#include "atom_vec_polar2D.h"

#include "atom.h"
#include "domain.h"
#include "memory.h"

#include <cmath>

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

AtomVecPolar2D::AtomVecPolar2D(LAMMPS *lmp) : AtomVec(lmp)
{
  molecular = Atom::ATOMIC;
  mass_type = PER_TYPE;

  atom->ang2D_flag = atom->ztorque_flag = 1;



  // strings with peratom variables to include in each AtomVec method
  // strings cannot contain fields in corresponding AtomVec default strings
  // order of fields in a string does not matter
  // except: fields_data_atom & fields_data_vel must match data file

  fields_grow = {"ang2D", "ztorque"};
  fields_copy = {"ang2D", "ztorque"};
  fields_comm = {"ang2D", "ztorque"};
  fields_border = {"ang2D", "ztorque"};
  fields_exchange = {"ang2D", "ztorque"};
  fields_create = {"ang2D", "ztorque"};
  fields_data_atom = {"id", "type", "ang2D", "x", "ztorque"};
  fields_data_vel = {"id", "v"};

  setup_fields();
}

/* ----------------------------------------------------------------------
   set local copies of all grow ptrs used by this class, except defaults
   needed in replicate when 2 atom classes exist and it calls pack_restart()
------------------------------------------------------------------------- */

void AtomVecPolar2D::grow_pointers()
{
  ang2D = atom->ang2D;
  ztorque = atom->ztorque;
}
void AtomVecPolar2D::data_atom_post(int ilocal)
{
}
