import numpy as np
from math import sqrt, pi
from ase.optimize import BFGS
from ase.lattice import bulk
from ase.lattice.cubic import BodyCenteredCubic
from ase.constraints import FixAtoms
import ase.units as units
from qlab import set_fortran_indexing, view
from quippy.atoms import Atoms
from quippy.structures import transform, supercell, MillerPlane, MillerDirection
from quippy.potential import Potential, Minim
from quippy.elasticity import youngs_modulus, poisson_ratio
from quippy.io import write
from atomeye import AtomEyeViewer
from quippy.crack import (print_crack_system,
                          G_to_strain,
                          thin_strip_displacement_y,
                          find_crack_tip_coordination,
                          find_crack_tip_stress_field,
                          make_crack_advance_map)

# additional requirements for the QM/MM simulation:
from quippy.potential import ForceMixingPotential
from quippy.lotf import LOTFDynamics, update_hysteretic_qm_region
#from gpaw import GPAW
from ase.calculators.lammps import LAMMPS
from quippy.io import AtomsWriter
# ******* Start of parameters ***********

# There are three possible crack systems, choose one and uncomment it

# System 1. (111)[0-11]
crack_direction = (-2, 1, 0)      # Miller index of x-axis
cleavage_plane = (1, 2, 0)        # Miller index of y-axis
crack_front = (0, 0, 1)          # Miller index of z-axis

# System 2. (110)[001]
# crack_direction = (1,-1,0)
# cleavage_plane = (1,1,0)
# crack_front = (0,0,1)

# System 3. (110)[1-10]
# crack_direction = (0,0,-1)
# cleavage_plane = (1,1,0)
# crack_front = (1,-1,0)

# Maximum force criteria for relaxation
relax_fmax = 0.001 * units.eV / units.Ang
# File to which structure will be written
output_file = 'crack-fe.xyz'
traj_file = 'traj-crack-fe.nc'
# ******* End of parameters *************

# set_fortran_indexing(False)

# ********** Build unit cell ************
ase_atoms = bulk('Fe', 'bcc', a=2.8553, cubic=True)

# ********** LAMMPS Setting ************
el1 = "Fe"
el2 = "Fe"
lp = 2.8553
species = [el1]
import potpara
print potpara.parameters
lammps_cal = LAMMPS(parameters=potpara.parameters, specorder=["Fe", "H"],
                    tmp_dir='lammps_tmp', keep_tmp_files=True)
# ********** Trasfer ASE to QUIP ************
# convert from LAMMPS calculator to QUIP potential
mm_pot = Potential(calculator=lammps_cal)
# mm_pot.set_default_quantities(['stresses'])
Fe_bulk = Atoms(ase_atoms)  # convert from ase.Atoms to quippy.Atoms
Fe_bulk.set_calculator(mm_pot)
# ***** Find eqm. lattice constant ******

# find the equilibrium lattice constant by minimising atoms wrt virial
# tensor
print('Minimising bulk unit cell...')
minim = Minim(Fe_bulk, relax_positions=True, relax_cell=True)
minim.run(fmax=relax_fmax)
print Fe_bulk.cell
a0 = Fe_bulk.cell[0, 0]
print('Lattice constant %.4f A\n' % a0)
# make a new bulk cell with correct a0 (so that off-diagonal lattice
# values are exactly zero)
Fe_bulk = bulk('Fe', 'bcc', a=a0, cubic=True)
Fe_bulk = Atoms(Fe_bulk)
Fe_bulk.set_calculator(mm_pot)
# ******* Find elastic constants *******
# Get 6x6 matrix of elastic constants C_ij
c = mm_pot.get_elastic_constants(Fe_bulk)
print('Elastic constants (GPa):')
print((c / units.GPa).round(2))
print('')
E = youngs_modulus(c, cleavage_plane)
print('Young\'s modulus %.1f GPa' % (E / units.GPa))
nu = poisson_ratio(c, cleavage_plane, crack_direction)
print('Poisson ratio % .3f\n' % nu)
# System 1. (111)[0-11]
GB_up_x = (0, 2, 1)      # Miller index of x-axis
GB_up_y = (0, -1, 2)        # Miller index of y-axis
GB_up_z = (1, 0, 0)          # Miller index of z-axis
GB_dn_x = (0, 2, -1)      # Miller index of x-axis
GB_dn_y = (0, 1, 2)        # Miller index of y-axis
GB_dn_z = (1, 0, 0)          # Miller index of z-axis
GB_up = BodyCenteredCubic(directions=[GB_up_x,
                                      GB_up_y,
                                      GB_up_z],
                          size=(1, 1, 1),
                          symbol='Fe',
                          pbc=True,
                          latticeconstant=a0)
GB_dn = BodyCenteredCubic(directions=[GB_dn_x,
                                      GB_dn_y,
                                      GB_dn_z],
                          size=(1, 1, 1),
                          symbol='Fe',
                          pbc=True,
                          latticeconstant=a0)
GB = GB_up + GB_dn
nl = NeighborList([2.3, 1.7])
nl.update(GB)
indices, offsets = nl.get_neighbors(0)
# view(GB)
