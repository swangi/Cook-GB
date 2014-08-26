import sys
import os
import errno
import shutil
import glob
import time
import numpy as np
from lammps import lammps
from ctypes import *
# import GBindex
# import matplotlib.pyplot as plt

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


def mkdirp(directory):
    if not os.path.isdir(directory):
        os.makedirs(directory)
# +-+-+-+-+-+
# |P|a|t|h|s|
# +-+-+-+-+-+
# +-+-+-+-+-+-+-+-+-+-+-+-+-+
# |T|s|u|j|i|-|C|l|u|s|t|e|r|
# +-+-+-+-+-+-+-+-+-+-+-+-+-+
# lammps_scr_path = "/home/wang/cal_scripts/lammps_scripts/"
# pot_path = "/home/wang/potentials/eam/Ram_song_potA.fs"
# +-+-+-+-+-+-+-+-+-+-+-+-+
# |A|l|e|x|-|c|l|u|s|t|e|r|
# +-+-+-+-+-+-+-+-+-+-+-+-+
lammps_scr_path = "/home/swang/share/lammps_scripts/"
pot_path = "/home/swang/share/potentials/eam/Ram_song_potA.fs"
# +-+-+-+-+-+-+
# |M|y|-|m|a|c|
# +-+-+-+-+-+-+
# lammps_scr_path = "/Users/swang/Simulation/LAMMPS/Py/lammps_scripts/"
# pot_path = "/Users/swang/Simulation/potentials/EAM/Ram_song_potA.fs"
# +-+-+-+-+-+-+-+-+-+-+ +-+-+-+ +-+-+-+-+-+-+
# |p|a|r|a|m|e|t|e|r|s| |f|o|r| |l|a|m|m|p|s|
# +-+-+-+-+-+-+-+-+-+-+ +-+-+-+ +-+-+-+-+-+-+
input_dir = 'input_dumps_60A'
n_A = 10
T = 300
T2 = T * 2
dump_id = "mydump"
dump_step = 1000
dump_name = "mydump"
h_coh = -2.36800333126474
fe_bulk = -4.155410449999999755
maxiter = 1000000
maxeval = 1000000
fixup = "172.766"
fixlow = "62.8917"
ftol = "1.0e-3"
etol = "0.0"
srate = "1.0e7"

displace_initial = 0.08
displace_final = 0.10
dis_target = 6
displace_delta = (displace_final - displace_initial) / (dis_target - 1.)
ten_step = 5000000
eq_step = 2000000
counter = 0
lmp = lammps("")  # , ["-screen", "none"])
pair = "eam/fs"
potfile = "%s Fe H" % (pot_path)
(nx, ny, nz) = (15, 1, 5)
cut_pos = "-0.637773"
n_frame = 1000  # size of dumped file
H_pos = [[12.53121, -0.96147, 2.37933],
        [18.91537, -0.96142, 2.37919],
        [6.14576, -0.96152, 2.37928]]


def read_lmp_line(name):
    # print lammps_scr_path + name
    lines = open(lammps_scr_path + name, 'r').readlines()
    for line in lines:
        lmp.command(line)


def lammps_initial(data):
    global num_fe, lx0
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
              "variable fixup equal %s" % (fixup),
              "variable fixlow equal %s" % (fixlow),
              "variable dump_name string %s" % (dump_name),
              "variable dump_id string %s" % (dump_id),
              "variable dump_step equal %d" % (dump_step),
              "variable srate equal %s" % (srate),
              "variable nx equal %s" % (nx),
              "variable lx0 equal lx",
              "units metal",
              "boundary p p p",
              "atom_style atomic"
              ]
    for r in ritual:
        lmp.command(r)
    lmp.command("read_data %s" % (data))
    lx0 = lmp.extract_variable("lx0", "all", 0)
    lmp.command("replicate %d %d %d" % (nx, ny, nz))
    read_lmp_line("pot_thermo.lammps")
    read_lmp_line("compute.lammps")
    lmp.command("group B type 1")
    num_fe = lmp.get_natoms()


def nvt(n_step, strain):
    read_lmp_line("nvt.lammps")
    dump_name = strain + "_eq_nvt"
    dump_step = n_step / 10
    lmp.command("variable dump_step equal %d" % (dump_step))
    lmp.command("variable dump_name string %s" % (dump_name))
    read_lmp_line("dump.lammps")
    lmp.command("run %d" % (n_step))
    lmp.command("undump %s" % (dump_id))
    lmp.command("unfix nvt")
    lmp.command("write_data %s" % (strain + "_eq.lammpstrj"))


def displace(init_sr):
    lmp.command("variable init_sr equal %.5f" % (init_sr))
    read_lmp_line("displace.lammps")
    # freeze 0.1 ly layer (up and low)
    # read_lmp_line("group.lammps")
    read_lmp_line("freeze.lammps")
    read_lmp_line("min.lammps")
    lmp.command("write_data %s" % ("min_crack.lammpstrj"))
    lmp.command("write_restart restart.min_crack")
    lmp.command("unfix freeze")


def cut(cut_pos):
    lmp.command("variable cut_pos equal %s" % (cut_pos))
    read_lmp_line("lxlylz.lammps")
    read_lmp_line("group.lammps")
    read_lmp_line("cut_potential.lammps")  # cut lx/nx distance for crack
    lmp.command("write_data %s" %
                ("initial_before_min.lammpstrj"))
    read_lmp_line("min.lammps")
    lmp.command("write_data %s" % ("min_initial.lammpstrj"))


def initial_npt(n_step):  # initial pressure equlibrilium
    lmp.command("reset_timestep 0")
    v_seed = int(time.time())
    lmp.command("variable v_seed equal %d" % (v_seed))
    read_lmp_line("initial_v.lammps")
    read_lmp_line("npt.lammps")
    dump_name = "eq_npt"
    dump_step = n_step / 50
    lmp.command("variable dump_step equal %d" % (dump_step))
    lmp.command("variable dump_name string %s" % (dump_name))
    read_lmp_line("dump.lammps")
    lmp.command("run %d" % (n_step))
    lmp.command("undump %s" % (dump_id))
    lmp.command("unfix npt")
    lmp.command("write_data %s" % ("eq_crack.lammpstrj"))
    lmp.command("write_restart restart.eq_crack")
    lmp.command("reset_timestep 0")


def main():
    if myid == 0:
        result_dirs = {'results_dir': "Results",
                       'traj_dir': "Results/Trajs",
                       'log_dir': "Results/Logs",
                       'restart_dir': "Results/Restarts"
                       }
        for k, v in result_dirs.items():
            mkdirp(v)
    temp_strain = displace_initial
    # START OF INITIAL
    lammps_initial(initial_trj)
    if 'H_pos' in globals():
        for p in H_pos:
            for n in range(nx + 1):
                lmp.command("create_atoms 2 single %5f %5f %5f units box" %
                            (n * lx0 + p[0], p[1], p[2]))
        read_lmp_line("msd.lammps")
    cut(cut_pos)
    lmp.command("log log.initial_npt append")
    initial_npt(eq_step)
    # INITIAL DONE
    # START OF DISPLACE
    # displace to a initial value to save time
    lmp.command("log log.displace_min append")
    displace(displace_initial)
    lmp.command("log log.displace_nvt append")
    nvt(eq_step, '{0:.2f}'.format(temp_strain))
    for i in range(0, dis_target):
        temp_strain += displace_delta
        lmp.command("log log.displace_min append")
        displace(displace_delta)
        lmp.command("log log.displace_nvt append")
        nvt(eq_step, '{0:.3f}'.format(temp_strain))

if __name__ == '__main__':
    # +-+-+-+-+-+-+ +-+-+-+-+-+-+ +-+-+-+ +-+-+-+-+-+-+-+-+-+ +-+-+
    # |A|d|j|u|s|t| |b|e|l|o|w|s| |f|o|r| |d|i|f|f|e|r|e|n|t| |G|B|
    # +-+-+-+-+-+-+ +-+-+-+-+-+-+ +-+-+-+ +-+-+-+-+-+-+-+-+-+ +-+-+
    file_names = [os.path.basename(x)
                  for x in glob.glob(input_dir + '/*.lammpstrj')]
    type_list = [i.split('.')[0] for i in file_names]
    typeGB = [i.split('_') for i in type_list]
    # h_pos_face = [0, 0, 0]
    for t in typeGB:
        GB_type = "_".join(t)
        GB_index = map(int, t)
        m = GB_index[0]
        n = GB_index[1]
        oe = 1
        if myid == 0:
            if (m + n) % 2 != 0:
                print "NOTE: m+n is odd number!"
                oe = 2
        initial_trj = input_dir + "/%s.lammpstrj" % (GB_type)
        #initial_trj = "initial_crack_0.02.lammpstrj"
        main()
        if myid == 0:
            traj_file_names = [glob.glob('*.lammpstrj')
                               ]
            log_file_names = [glob.glob('log.*')
                              ]
            restart_file_names = [glob.glob('restart.*')
                                  ]
            for traj_file in traj_file_names:
                os.rename(traj_file, result_dirs['traj_dir'] + "/" + traj_file)
            for log_file in log_file_names:
                os.rename(log_file, result_dirs['log_dir'] + "/" + log_file)
            for restart_file in restart_file_names:
                os.rename(traj_file, result_dirs[
                          'restart_dir'] + "/" + restart_file)
    print('====================================Fin====================================')

        # +-+-+-+-+ +-+-+-+-+-+-+-+ +-+-+-+-+-+-+
        # |S|o|m|e| |i|n|i|t|i|a|l| |v|a|l|u|e|s|
        # +-+-+-+-+ +-+-+-+-+-+-+-+ +-+-+-+-+-+-+
    # print "Proc %d out of %d procs has" % (myid, nprocs), lmp
    # pypar.finalize()
