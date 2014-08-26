import numpy as np
import shutil
import os
import glob
from math import sqrt, pi, ceil
from ase.optimize.sciopt import SciPyFminCG
from ase.lattice.cubic import BodyCenteredCubic
from ase.lattice.surface import add_adsorbate
from ase.constraints import FixAtoms
import ase.units as units
from ase.calculators.lammps import LAMMPS, write_lammps_data
from qlab import Atoms
from qlab import view
import ase.units as units
from quippy.io import write
from quippy.potential import Potential, Minim
from quippy.farray import frange, convert_farray_to_ndarray, farray
from quippy.structures import bcc, MillerIndex, supercell
import mpi4py
import time
start = time.time()
a0 = 2.8557332
x0 = 1
y0 = 5
grid_x = 4
grid_z = 4
height = 60.0 * units.Ang
vacuum = 20.0 * units.Ang
eng_tol = '1e-15'
force_tol = '1e-15'
relax_fmax = 0.00001 * units.eV / units.Ang
GB_o_cor = [[0, 1, 0], [0, 0, 1], [1, 0, 0]]

# Functions #########################################

# in this function the last para, 0 means up, 1 means dn atoms


def find_atoms_within_cutoff(atoms, cutoff, d=0):
    tmp_atoms = Atoms(atoms, fortran_indexing=False)
    tmp_atoms.calc_connect()
    over = []
    rij = farray(0.0)
    for i in frange(len(tmp_atoms)):
        for n in frange(tmp_atoms.n_neighbours(i)):
            j = tmp_atoms.neighbour(i, n, distance=rij)
            # print j
            if rij < cutoff and i > j:
                if d == 0:
                    over.append(j - 1)
                else:
                    over.append(i - 1)

    # for i in range(len(atoms)): #python method, may be slower than fortran index?
    #     indices, offsets = tmp_atoms.neighbours.get_neighbors(i)
    #     for j, offset in zip(indices, offsets):
    #         diff=tmp_atoms.get_distance(i,j)
    # if diff < cutoff and i > j : #i>j is important for returning only one atom id
    #             over.append(j)
    return over


def create_unit_GB(cor, a):
    unit_GB = BodyCenteredCubic(directions=[cor[0], cor[1], cor[2]],
                                size=(1, 1, 1),
                                symbol='Fe',
                                pbc=(1, 1, 1),
                                latticeconstant=a)
    return unit_GB


def cal_surface_energy(atoms):
    surface = atoms.copy()
    surface.set_calculator(mm_pot)
    E_per_atom_bulk = surface.get_potential_energy() / len(surface)
    surface.center(vacuum, axis=1)
    # potpara.parameters["minimize"] = "0 1.0e-6 10000 10000"
    E_surf = surface.get_potential_energy()
    area = surface.get_volume() / surface.cell[1, 1]
    gamma = ((E_surf - E_per_atom_bulk * len(surface)) /
             (2.0 * area))
    return gamma / (units.J / units.m ** 2)


def build_GB(unit_GB_up, unit_GB_down, tx, tz):
    ny = int(ceil(height / unit_GB_up.cell[1, 1]))
    # print ny
    ny2 = ny * 2.
    GB_up = supercell(unit_GB_up, 2, ny2, 6)
    GB_up = GB_up.select(GB_up.positions[:, 1] >= GB_up.positions[:, 1].mean())
    print "created %d atoms in up grain" % (len(GB_up))
    GB_down = supercell(unit_GB_down, 2, ny2, 6)
    GB_down = GB_down.select(
        GB_down.positions[:, 1] <= ny * unit_GB_down.cell[1, 1])
    print "created %d atoms in down grain" % (len(GB_down))
    # GB_up.positions[:,1]+= GB_up.positions[:, 1].max()
    GB_down.positions[:, 0] += tx
    GB_down.positions[:, 2] += tz
    GB = GB_down + GB_up
    # GB.set_cell([GB_down.cell[0],GB_down.cell[1]*2,GB_down.cell[2]])
    return GB


def lammps_atoms(GB):
    GB = lammps_cal.atoms
    GB = Atoms(GB)
    GB.set_calculator(mm_pot)
    return GB


def euler(GB_o):
    trans_matrix_up = [[y0, -x0, 0], [x0, y0, 0], [0, 0, 1]]
    trans_matrix_down = [[y0, x0, 0], [-x0, y0, 0], [0, 0, 1]]
    GB_cor_up = np.dot(trans_matrix_up, GB_o)
    GB_cor_down = np.dot(trans_matrix_down, GB_o)
    return GB_cor_up, GB_cor_down
# variables############################
GB_cor_up, GB_cor_down = euler(GB_o_cor)
count = 0
nm_GB = 0

tmp_dir = 'lammps_tmp'
rangle = np.degrees(np.arctan2(x0, y0))
print rangle

print('up grain is like :\n %s ' % (GB_cor_up))
print('down grain is like :\n %s ' % (GB_cor_down))

GB_plane = (GB_cor_up[1, 0], GB_cor_up[1, 1], GB_cor_up[1, 2])
Sigma = GB_cor_up[1, 0] ** 2 + GB_cor_up[1, 1] ** 2 + GB_cor_up[1, 2] ** 2
while Sigma % 2 == 0:
    Sigma = Sigma / 2
DSC = (GB_cor_up[1, 0] ** 2 + GB_cor_up[1, 1] ** 2 + GB_cor_up[1, 2] ** 2)

xinc = 1. / grid_x / DSC
zinc = 1. / grid_z
overlapdist = np.arange(0.275, 0.705, 0.005) * a0
lable = 'S%d_%d%d%d_' % (Sigma, GB_plane[0], GB_plane[1], GB_plane[2])
traj_name = 'Fe_GB_' + lable
el1 = "Fe"
el2 = "H"
species = [el1, el2]
import potpara
print potpara.parameters
#, keep_tmp_files=False)
lammps_cal = LAMMPS(
    parameters=potpara.parameters, specorder=species, tmp_dir=tmp_dir)
mm_pot = Potential(calculator=lammps_cal)
f = open('GB_result_S%d_%d%d%d_.txt' %
         (Sigma, GB_plane[0], GB_plane[1], GB_plane[2]), 'w')
f.write('%s %s %s %s %s %s %s %s %s %s\n' %
       ('Angle', 'Count', 'del_cri', 'RBT_X', 'RBT_Z', 'Bulk_pe_GB', 'N_atoms', 'Surface_e', 'GB_e', 'GB_coh_e'))
f.close()
# create unit cell for GB#################
unit_GB_up = create_unit_GB(GB_cor_up, a0)
unit_GB_down = create_unit_GB(GB_cor_down, a0)
# view(unit_GB)
unit_GB_up = Atoms(unit_GB_up, fortran_indexing=False)
unit_GB_up.set_calculator(mm_pot)
print('Minimising GB_up unit cell...')
potpara.parameters["fix"] = "relax all box/relax y 0 vmax 0.001"
potpara.parameters["minimize"] = "%s %s 10000 10000" % (eng_tol, force_tol)
E_per_atom_bulk = unit_GB_up.get_potential_energy() / len(unit_GB_up)
print "Get atom positions......"
unit_GB_up = lammps_atoms(unit_GB_up)
tx_inc = xinc * unit_GB_up.get_cell()[0, 0]
tz_inc = zinc * unit_GB_up.get_cell()[2, 2] / 2
print tx_inc, tz_inc
print('Energy per atom of up GB is: %f eV' % (E_per_atom_bulk))
surf_eng = cal_surface_energy(unit_GB_up)
print('Surface energy of %s surface is %.4f J/m^2\n' %
     (GB_plane, surf_eng))

unit_GB_down = Atoms(unit_GB_down, fortran_indexing=False)
unit_GB_down.set_calculator(mm_pot)
print('Minimising bulk unit cell...')
E_per_atom_bulk_down = unit_GB_down.get_potential_energy() / len(unit_GB_down)
print "Get atom positions......"
unit_GB_down = lammps_atoms(unit_GB_down)
print('Energy per atom of down GB is: %f eV' % (E_per_atom_bulk_down))
del potpara.parameters["minimize"]
del potpara.parameters["fix"]
# copy GB unit to form bi-crystal ######
for i in range(0, grid_x + 1):
    tx = i * tx_inc
    for j in range(0, grid_z + 1):
        tz = j * tz_inc
        print "x move %.4f" % (tx)
        print "z move %.4f" % (tz)
        for cut in overlapdist:
            if cut == 0.275 * a0:
                natom_pre = 1
            for d in range(1):
                traj_name_new = traj_name + \
                    '_x%d_z%d_cutoff%.4f_%d' % (i, j, cut, d)
                print(traj_name_new + '......')
                GB = build_GB(unit_GB_up, unit_GB_down, tx, tz)
# delete overlaped atoms ###############
                if d == 0:
                    print "Eating all overlaped atoms (UP)....."
                else:
                    print "Eating all overlaped atoms (DOWN)....."
                over = find_atoms_within_cutoff(GB, cut, d)
                print "Eat %d overlaped atoms: " % (len(over)), over
                if over != []:
                    del GB[over]
                count += 1
                print "%d We eat some atoms, before: %d; now: %d " % (
                    count, natom_pre, len(GB))
                if len(GB) == natom_pre:
                    print "We must eat more...\n"
                    continue
                nm_GB += 1
                print "We found %d GB(s)! His number is %d! Send it to troop!" % (
                    nm_GB, count)
                natom_pre0 = len(GB)
                natom_pre = natom_pre0
# minimization #########################
                # GB.pbc = True
                GB = Atoms(GB, fortran_indexing=False)
                # print GB.get_pbc()
                GB.set_calculator(mm_pot)
                print "First Minimization....."
                potpara.parameters["minimize"] = "%s %s 10000 10000" % (
                    eng_tol, force_tol)
                bulk_GB_pe0 = GB.get_potential_energy()
                print "Get atom positions......"
                GB = lammps_atoms(GB)
                del potpara.parameters["minimize"]
                print "Second Minimization....."
                potpara.parameters[
                    "fix"] = "relax all box/relax y 0.0 vmax 0.001"
                potpara.parameters["minimize"] = "%s %s 10000 10000" % (
                    eng_tol, force_tol)
                bulk_GB_pe = GB.get_potential_energy()
                print "Get atom positions......"
                GB = lammps_atoms(GB)
                tol_atoms = len(GB)
                GB_eng = (bulk_GB_pe - E_per_atom_bulk * tol_atoms) / (
                    2.0 * GB.get_volume() / GB.cell[1, 1]) / (
                        units.J / units.m ** 2)
                print('Total energy of GB is %.6f' % (bulk_GB_pe))
                print('GB energy of S%d(%d%d%d) is %.6f J/m^2\n' %
                     (Sigma, GB_plane[0], GB_plane[1], GB_plane[2], GB_eng))

                # sort Traj with GB energy mj/m^2
                traj_sort = int(GB_eng * 1000)
                traj_name_new = '%d_' % (
                    traj_sort) + '%d_' % (count) + traj_name_new

                GB_cohesive_eng = 2 * surf_eng - GB_eng
                f = open('GB_result_S%d(%d%d%d).txt' %
                        (Sigma, GB_plane[0], GB_plane[1], GB_plane[2]), 'a')
                f.write('%.6f %d %.4f %.6f %.6f %.6f %d %.6f %.6f %.6f\n' %
                       (rangle, count, cut, tx, tz, bulk_GB_pe,
                        tol_atoms, surf_eng, GB_eng, GB_cohesive_eng))
                f.close()
                del potpara.parameters["minimize"]
                del potpara.parameters["fix"]
                GB.add_property('ack', lammps_cal.get_ack(GB))
# writing-results#######################
                print "Writing traj file....."
                if not os.path.exists(traj_name):
                    os.mkdir(traj_name)
                cwd = os.getcwd()
                os.chdir(traj_name)
                # write_lammps_data('out.lammps',GB,['Fe'])
                write(traj_name_new + '.xyz', GB)
                lammps_cal.write_lammps_traj(
                    traj_name=traj_name_new + '.lammpstrj')
                lammps_cal.write_lammps_data(
                    lammps_data=traj_name_new + '.lammpsdata')

                os.chdir('../' + tmp_dir)
                filelist = glob.glob("*lammps*")
                for f in filelist:
                    os.remove(f)
                os.chdir(cwd)
                # view(GB)

end = time.time()
mi = (end - start) / 60
print "All done! Takes %f min." % (mi)
