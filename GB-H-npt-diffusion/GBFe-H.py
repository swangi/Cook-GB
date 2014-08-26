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
import pypar
myid = pypar.rank()
nprocs = pypar.size()
node = pypar.get_processor_name()
# +-+-+-+-+-+-+-+-+-+-+ +-+-+-+ +-+-+-+-+-+-+
# |p|a|r|a|m|e|t|e|r|s| |f|o|r| |l|a|m|m|p|s|
# +-+-+-+-+-+-+-+-+-+-+ +-+-+-+ +-+-+-+-+-+-+
initial_trj = "initial.lammpstrj"
vac = 2
n_A = 10
T = 300
T2 = T*2
maxiter = 1000000
maxeval = 1000000
ftol = "1.0e-15"
etol = "1.0e-15"
surfelist = []

counter = 0
# if not os.path.exists("Results"): os.makedirs("Results")
if myid == 0:
    mkdir_p("Results")

# +-+-+-+-+-+-+-+-+-+-+
# |P|o|t|e|n|t|i|a|l|s|
# +-+-+-+-+-+-+-+-+-+-+
# -------------------Input EAM pot---------------------------------
pair = "eam/fs"
potfile = "/home/wang/potentials/eam/PotentialA3410-song.fs Fe H"
# potfile="/home/wang/potentials/eam/PotentialB3410.fs Fe H"
# ------------------Input MEAM pot--------------------------------
# pair="meam"
# potfile="/home/wang/potentials/meam/lib-fe-p.meam Fe /home/wang/potentials/meam/Fe-meam-p.meam Fe"


# +-+-+-+-+-+ +-+-+-+-+-+-+-+-+-+-+-+-+-+
# |S|t|a|r|t| |l|a|m|m|p|s|r|u|n|n|i|n|g|
# +-+-+-+-+-+ +-+-+-+-+-+-+-+-+-+-+-+-+-+
# screen is not good look, rev me when needed
lmp = lammps("") #, ["-screen", "none"])

lmp.command("variable pair string %s" % (pair))
lmp.command('variable potfile string "%s"' % (potfile))
lmp.command("variable maxiter equal %d" % (maxiter))
lmp.command("variable maxeval equal %d" % (maxeval))
lmp.command("variable ftol equal %s" % (ftol))
lmp.command("variable T equal %s" % (T))
lmp.command("variable T2 equal %s" % (T2))
lmp.command("variable etol equal %s" % (etol))
lmp.command("read_data %s" % (initial_trj))
lmp.command("region B block EDGE EDGE INF INF INF INF")
lmp.command("group B region B")
# lmp.command("change_box all x delta -%d %d boundary p p p remap" %
#             (vac, vac))
lmp.command("region A1 block -%d 0 INF INF INF INF" %
            (vac))
# lmp.command("delete_atoms region A1")
# lmp.command("region A2 block -%d %d INF INF INF INF" %
            # (vac))
#lmp.command("create_atoms 2 random %d 121434 A1" % (n_A))
#lmp.command("group A type 2")
lines = open("lammps_scripts/npt.lammps", 'r').readlines()
for line in lines:
    lmp.command(line)
lmp.command("run 1000")
# uncomment if running in parallel via Pypar
print "Proc %d out of %d procs has" % (myid, nprocs), lmp
pypar.finalize()






# infile = "in.trivial"
#
#
#
# from lammps import lammps
# lmp = lammps()
#
# run infile all at once
#
# lmp.file(infile)
#
# run infile one line at a time
#
# lines = open(infile,'r').readlines()
# for line in lines: lmp.command(line)
#
