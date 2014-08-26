#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (c) 2014 - Swang <swangi@outlook.com>
# Filename:
# +-+-+-+-+-+-+ +-+-+-+-+-+-+-+-+-+
# |I|m|p|o|r|t| |s|o|m|e|t|h|i|n|g|
# +-+-+-+-+-+-+ +-+-+-+-+-+-+-+-+-+
#
from pylab import *
import os
import glob
import time
import errno
import shutil
import sys
import re


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
#
# +-+-+-+-+-+ +-+-+ +-+-+-+-+
# |S|t|a|r|t| |o|f| |m|a|i|n|
# +-+-+-+-+-+ +-+-+ +-+-+-+-+
#
from lammps import lammps
lmp = lammps("", ["-screen", "none"])
lammps_scr_path = "/Users/swang/Simulation/LAMMPS/Py/lammps_scripts/"
pot_path = "/Users/swang/Simulation/potentials/EAM/Ram_song_potA.fs"
pair = "eam/fs"
potfile = "%s Fe H" % (pot_path)


def read_lmp_line(name):
    # print lammps_scr_path + name
    lines = open(lammps_scr_path + name, 'r').readlines()
    for line in lines:
        lmp.command(line)


def dump2data(input_data, output_data):
    global num_fe
    ritual = ["clear",
              "atom_modify sort 0 0.0",
              "units metal",
              "boundary p p p",
              "atom_style atomic",
              'lattice hcp 2.8553',
              'region whole block 0  10 0 10 0 10',
              'create_box 2 whole',
              "variable pair string %s" % (pair),
              'variable potfile string "%s"' % (potfile),
              ]
    for r in ritual:
        lmp.command(r)
    lmp.command("read_dump %s 0 x y z box yes add yes" % input_data)
    lmp.command('mass 1 55.847')
    lmp.command('mass 2 1.008')
    # lmp.command("replicate %d %d %d" % (nx, ny, nz))
    # read_lmp_line("pot_thermo.lammps")
    # read_lmp_line("compute.lammps")
    # lmp.command('shell cd %s' % output_dir)
    lmp.command("write_data %s" % output_data)
    # lmp.command('shell cd ..')
    # num_fe = lmp.get_natoms()


def read_type(input_files):
    file_names = [os.path.basename(x) for x in glob.glob(input_files)]
    name = [re.sub('[(){}<>]', '', f) for f in file_names]
    typeGB = [i.split('_') for i in name]
    return file_names, typeGB

if __name__ == "__main__":
    output_dir = 'input_dumps/'
    mkdir_p(output_dir)
    input_dir = 'results/'
    dump_file, data_pre = read_type(input_dir + '*.lammpstrj')
    for i, t in enumerate(data_pre):
        input_file = input_dir + dump_file[i]
        m = t[2]
        n = t[4]
        data_file = m+"_"+m+"_"+n
        output_data = output_dir + data_file + '.lammpstrj'
        dump2data(input_file, output_data)
