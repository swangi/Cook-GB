import sys
import os
import errno
import shutil
import glob
import time
import numpy as np
from lammps import lammps
# from ctypes import *
# import GBindex
# import matplotlib.pyplot as plt


def mkdir_p(path):
    try:
        os.makedirs(path)
    except os.error as e:
        if e.errno != 17:
            raise
        else:
            pass


def rename_p(old, new):
    try:
        os.rename(old, new)
    except os.error as e:
        if e.errno != 2:
            raise
        else:
            pass
# +-+-+-+ +-+-+-+-+ +-+-+-+-+-+
# |r|u|n| |w|i|t|h| |p|y|p|a|r|
# +-+-+-+ +-+-+-+-+ +-+-+-+-+-+
# import pypar
# myid = pypar.rank()
# nprocs = pypar.size()
# node = pypar.get_processor_name()
from mpi4py import MPI
comm = MPI.COMM_WORLD
myid = comm.Get_rank()
nprocs = comm.Get_size()
node = MPI.Get_processor_name()
# # +-+-+-+-+-+
# |P|a|t|h|s|
# +-+-+-+-+-+
lammps_scr_path = "/home/wang/cal_scripts/lammps_scripts/"
pot_path = "/home/wang/potentials/eam/Ram_song_potA.fs"
# lammps_scr_path = "/home/swang/share/lammps_scripts/"
# pot_path = "/home/swang/share/potentials/eam/Ram_song_potA.fs"
# lammps_scr_path = "/Users/swang/Simulation/LAMMPS/Py/lammps_scripts/"
# pot_path = "/Users/swang/Simulation/potentials/EAM/Ram_song_potA.fs"
# +-+-+-+-+-+-+-+-+-+-+ +-+-+-+ +-+-+-+-+-+-+
# |p|a|r|a|m|e|t|e|r|s| |f|o|r| |l|a|m|m|p|s|
# +-+-+-+-+-+-+-+-+-+-+ +-+-+-+ +-+-+-+-+-+-+
n_A = 10
T = 300
T2 = T * 2
dump_id = "mydump"
dump_step = 1
dump_name = "mydump"
h_coh = -2.36800333126474
fe_bulk = -4.155410449999999755
maxiter = 1000000
maxeval = 1000000
ftol = "1.0e-5"
etol = "0.0"
counter = 0
lmp = lammps("", ["-screen", "none"])
pair = "eam/fs"
potfile = "%s Fe H" % (pot_path)


def read_lmp_line(name):
    # print lammps_scr_path + name
    lines = open(lammps_scr_path + name, 'r').readlines()
    for line in lines:
        lmp.command(line)


def lammps_initial(data):
    global num_fe
    ritual = ["clear",
              "atom_modify sort 0 0.0",
              "variable pair string %s" % (pair),
              'variable potfile string "%s"' % (potfile),
              "variable maxiter equal %d" % (maxiter),
              "variable maxeval equal %d" % (maxeval),
              "variable ftol equal %s" % (ftol),
              "variable etol equal %s" % (etol),
              "variable T equal %s" % (T),
              "variable T2 equal %s" % (T2),
              "variable dump_name string %s" % (dump_name),
              "variable dump_id string %s" % (dump_id),
              "variable dump_step equal %d" % (dump_step),
              "units metal",
              "boundary p p p",
              "atom_style atomic"
              ]
    for r in ritual:
        lmp.command(r)
    lmp.command("read_data %s" % (data))
    lmp.command("replicate %d %d %d" % (nx, ny, nz))
    read_lmp_line("pot_thermo.lammps")
    read_lmp_line("compute.lammps")
    lmp.command("group B type 1")
    num_fe = lmp.get_natoms()
#
# GB energy without H calculation ###################
#


def GB_eng():
    # dump_name = "%s_GB" % (GB_type)
    # lmp.command("variable dump_name string %s" % (dump_name))
    # read_lmp_line("dump.lammps")
    read_lmp_line("min_GB.lammps")
    pote_gb = lmp.extract_variable("pot", "all", 0)
    gbe_0 = (pote_gb - (fe_bulk * num_fe)
             ) / lmp.extract_variable("area", "all", 0) / 2
    gbemJm2_0 = gbe_0 * 16021.7733
    area_gb = lmp.extract_variable("area", "all", 0)
    lmp.command("write_data %s" % (GB_traj))
    return gbe_0, gbemJm2_0, area_gb, pote_gb


def FS_eng(gbe):
    # lmp.command("undump %s" % (dump_id))
    # dump_name = "%s_FS" % (GB_type)
    # lmp.command("variable dump_name string %s" % (dump_name))
    # read_lmp_line("dump.lammps")
    read_lmp_line("FS_create.lammps")
    read_lmp_line("min.lammps")
    pote_fs = lmp.extract_variable("pot", "all", 0)
    fse_0 = (pote_fs - (fe_bulk * num_fe)
             ) / lmp.extract_variable("area", "all", 0) - gbe
    fsemJm2_0 = fse_0 * 16021.7733
    area_fs = lmp.extract_variable("area", "all", 0)
    lmp.command("write_data %s" % (FS_traj))
    return fse_0, fsemJm2_0, area_fs, pote_fs


def interface_H_min(gamma_int, gamma_prefix):
    # global h_pos_face
    h_pos_face = [0, 0, 0]
    data_pre = gamma_prefix
    lmp.command("create_atoms 2 single %5f %5f %5f units box" %
                (pos1[0], pos1[1], pos1[2]))
    # print pos1
    lmp.command("group A type 2")
    # comment following line if you want to relax hydrogen
    if freeze == 1:
        lmp.command("fix freeze A setforce 0.0 0.0 0.0")
    lmp.command("delete_atoms overlap 0.1 A B")
    num_tol = lmp.get_natoms()
    num_h = lmp.get_natoms() - num_fe
    if myid == 0:
        print "number of H is %d" % (num_h)
    if num_h == 0:
        return 0, 0., 0., pos1
    read_lmp_line("min.lammps")
    pote_h_face = lmp.extract_variable("pot", "all", 0)
    ed = pote_h_face - gamma_int - h_coh * num_h
    area_face_H = lmp.extract_variable("area", "all", 0)
    f = "temp.lammpstrj"
    lmp.command("write_data %s" %
               (f))
    if myid == 0:
        with open(f, 'r') as f:
                    lines = f.readlines()
                    # print lines
        # lines = pypar.broadcast(lines, 0)
        for i in lines[16:num_tol + 16]:
            if int(i.split()[1]) == 2:
                h_pos_face = map(float, i.split()[2:5])
        data_name_p = "%s_%.2fx_%.2fy_%.2fz.lammpstrj" % (
            data_pre, h_pos_face[0], h_pos_face[1], h_pos_face[2])
        rename_p("temp.lammpstrj", '%s/%s' % (Result_dir, data_name_p))
    #     for i in range(1, pypar.size()):
    #         pypar.send(h_pos_face, i)
    # h_pos_face = pypar.receive(0)
    return num_h, ed, area_face_H, h_pos_face


def make_traj(Btrj, log):
    file_B = Btrj  # solvent
    file_A = log  # solute
    with open(file_B, 'r') as f:
        linesb = f.readlines()
    with open(file_A, 'r') as f:
        linesa = f.readlines()
    num_B = int(linesb[2].split()[0])
    # print num_B
    data_A = []
    data_B = []
    for i in linesa[1:]:
        data_A.append(map(float, i.split()[0:8]))
    for i in linesb[16:num_B + 16]:
        data_B.append(map(float, i.split()[0:8]))
    num_A = len(data_A)
    num_all = num_A + num_B
    output_file = file_B + '_H_linescan.lammpstrj'
    # int(data_B[:2])
    with open(output_file, 'w') as f:
        f.write(
            "ITEM: TIMESTEP\n0\nITEM: NUMBER OF ATOMS\n%d\nITEM: BOX BOUNDS pp pp pp\n%.5f %.5f\n%.5f %.5f\n%.5f %.5f\nITEM: ATOMS id type x y z Ed_H area layer_n\n" %
            (num_all,
                float(linesb[5].split()[0]), float(linesb[5].split()[1]),
                float(linesb[6].split()[0]), float(linesb[6].split()[1]),
                float(linesb[7].split()[0]), float(linesb[7].split()[1])))
        # print data_B[0:20]
        for i in data_B:
            f.write('\t'.join(['%d' % (n) for n in i[0:2]]))
            f.write('\t')
            f.write('\t'.join(['%f' % (n) for n in i[2:]]))
            f.write("\n")
        for i in data_A:
            f.write('\t'.join(['%d' % (n) for n in i[0:2]]))
            f.write('\t')
            f.write('\t'.join(['%f' % (n) for n in i[2:]]))
            f.write("\n")
    rename_p(output_file, '%s/%s' % (Result_dir, output_file))


def main():
    lammps_initial(initial_trj)
    gbe_0, gbemJm2_0, area_gb, pote_gb = GB_eng()
    fse_0, fsemJm2_0, area_fs, pote_fs = FS_eng(gbe_0)
    N_GB = 0
    N_FS = 0
    global H_pos_GB
    for k in range(0, zresol):
        for j in range(-yresol, yresol + 1):
            for i in range(0, xresol + 1):
                dxx = dx * i
                dyy = dy * j
                dzz = dz * k
                pos1[0] = x0 + dxx
                pos1[1] = y0 + dyy
                pos1[2] = z0 + dzz
                # print pos1[0]
                lammps_initial(GB_traj)
                num_h, ed_h_gb, h_area_gb, h_pos_gb = interface_H_min(
                    pote_gb, "GB")
                if num_h == 0:
                    continue
                if myid == 0:
                    if H_pos_GB == []:
                        H_pos_GB.append(h_pos_gb)
                        nn_tol = [1]
                    else:
                        nn_tol = [np.linalg.norm(
                            np.array(h_pos_gb) - np.array(pos)
                        ) for pos in H_pos_GB]
                        # print nn_tol
                    if all(np.array(nn_tol) > 0.05):
                        Ed_H_GB.append(ed_h_gb)
                        # print H_pos_GB
                        H_pos_GB.append(h_pos_gb)
                        # H_area_GB.append(h_area_gb)
                        N_GB += 1
                        id_H_gb = num_fe + N_GB
                        print "NO. ~%d~  H position:[%.5f, %.5f, %.5f] at %sGB (%.5f mJm2); Ed_H^GB : %.5f eV" % (
                            N_GB, h_pos_gb[0], h_pos_gb[1], h_pos_gb[2],
                            GB_type, gbemJm2_0, ed_h_gb)
                        with open(GB_log, 'a') as f:
                            f.write(
                                '%d\t2\t%.5f\t%.5f\t%.5f\t%.5f\t%.5f\t%d\n' %
                                (id_H_gb, h_pos_gb[
                                    0], h_pos_gb[1], h_pos_gb[2],
                                    ed_h_gb, h_area_gb, k + 1))
    for k in range(0, zresol):
        for j in range(-yresol, 1):
            for i in range(0, xresol + 1):
                dxx = dx * i
                dyy = dy * j
                dzz = dz * k
                pos1[0] = x0 + dxx
                pos1[1] = y0 + dyy
                pos1[2] = z0 + dzz
                lammps_initial(FS_traj)
                num_h, ed_h_fs, h_area_fs, h_pos_fs = interface_H_min(
                    pote_fs, "FS")
                if num_h == 0:
                    continue
                if myid == 0:
                    if H_pos_FS == []:
                        H_pos_FS.append(h_pos_fs)
                        nn_tol = [1]
                    else:
                        nn_tol = [np.linalg.norm(np.array(h_pos_fs)
                                                 - np.array(pos)
                                                 ) for pos in H_pos_FS]
                        # print nn_tol
                    if all(np.array(nn_tol) > 0.05):
                        Ed_H_FS.append(ed_h_fs)
                        H_pos_FS.append(h_pos_fs)
                        # H_area_FS.append(h_area_fs)
                        N_FS += 1
                        id_H_fs = num_fe + N_FS
                        print "NO. ~%d~ H position:[%.5f, %.5f, %.5f] at %sFS (%.5f mJm2); Ed_H^FS : %.5f eV" % (
                            N_FS, h_pos_fs[0], h_pos_fs[1], h_pos_fs[2],
                            GB_type, fsemJm2_0, ed_h_fs)
                        with open(FS_log, 'a') as f:
                            f.write(
                                '%d\t2\t%.5f\t%.5f\t%.5f\t%.5f\t%.5f\t%.d\n' %
                                (id_H_fs, h_pos_fs[
                                    0], h_pos_fs[1], h_pos_fs[2],
                                    ed_h_fs, h_area_fs, k + 1))
    if myid == 0:
        mu_h_poor_GB = min(Ed_H_GB) + h_coh
        mu_h_poor_FS = min(Ed_H_FS) + h_coh
        print "GB type %s finished, H poor case GB energy is %.5f eV" % (
            GB_type, mu_h_poor_GB)
        print "FS type %s finished, H poor case FS energy is %.5f eV" % (
            GB_type, mu_h_poor_FS)
if __name__ == '__main__':
    # +-+-+-+-+-+-+ +-+-+-+-+-+-+ +-+-+-+ +-+-+-+-+-+-+-+-+-+ +-+-+
    # |A|d|j|u|s|t| |b|e|l|o|w|s| |f|o|r| |d|i|f|f|e|r|e|n|t| |G|B|
    # +-+-+-+-+-+-+ +-+-+-+-+-+-+ +-+-+-+ +-+-+-+-+-+-+-+-+-+ +-+-+
    file_names = [os.path.basename(x)
                  for x in glob.glob('input_dumps/*.lammpstrj')]
    type_list = [i.split('.')[0] for i in file_names]
    typeGB = [i.split('_') for i in type_list]
    # h_pos_face = [0, 0, 0]
    for t in typeGB:
        GB_type = "_".join(t)
        GB_index = map(int, t)
        m = GB_index[0]
        n = GB_index[2]  # need to adjust for dif tilt axis
        oe = 1
        # if myid == 0:
        #     if (m + n) % 2 != 0:
        #         print "NOTE: m+n is odd number!"
        #         oe = 2
        initial_trj = "input_dumps/%s.lammpstrj" % (GB_type)
        FS_traj = "%s_FS.lammpstrj" % (GB_type)
        GB_traj = "%s_GB.lammpstrj" % (GB_type)
        a0 = 2.855730057
        lz = a0 * np.sqrt(2.0)  # need to adjust for dif tilt axis
        # need to adjust for dif tilt axis
        lx = a0 * np.sqrt(n ** 2 + n ** 2 + (2 * m) ** 2) * 3  # / 2.0
        # need to adjust for dif tilt axis
        dy = a0 / np.sqrt(m ** 2 + m ** 2 + n ** 2)  # / oe
        y0 = 0.000
        z0 = 0.000
        x0 = 0.000
        pos1 = [x0, y0, z0]
        freeze = 0

        # +-+-+-+-+ +-+-+-+-+-+-+-+ +-+-+-+-+-+-+
        # |S|o|m|e| |i|n|i|t|i|a|l| |v|a|l|u|e|s|
        # +-+-+-+-+ +-+-+-+-+-+-+-+ +-+-+-+-+-+-+
        zresol = 4
        xresol = 200
        yresol = 30
        dx = lx / xresol
        dz = lz / zresol
        nx = 1
        ny = 1
        nz = 1
        Ed_H_GB = []
        Ed_H_FS = []
        H_pos_GB = []
        H_pos_FS = []
        H_area_GB = []
        H_area_FS = []
        surfelist = []
        msd = []
        GB_log = '%s_Dissolution_H_GB.txt' % (GB_type)
        FS_log = '%s_Dissolution_H_FS.txt' % (GB_type)
        Result_dir = "%s_Results" % (GB_type)
        if myid == 0:
            with open(GB_log, 'w') as f:
                f.write('#id\ttype\tx\ty\tz\tEdH\tarea\tlayer_n\n')
            with open(FS_log, 'w') as f:
                f.write('#id\ttype\tx\ty\tz\tEdH\tarea\tlayer_n\n')
            mkdir_p(Result_dir)
        main()
        if myid == 0:
            make_traj(GB_traj, GB_log)
            make_traj(FS_traj, FS_log)
            rename_p(GB_traj, '%s/%s' % (Result_dir, GB_traj))
            rename_p(FS_traj, '%s/%s' % (Result_dir, FS_traj))
            rename_p(GB_log, '%s/%s' % (Result_dir, GB_log))
            rename_p(FS_log, '%s/%s' % (Result_dir, FS_log))
    # print "Proc %d out of %d procs has" % (myid, nprocs), lmp
    # pypar.finalize()
