import sys
import os
import errno
import shutil
import glob
import time
import numpy as np
from lammps import lammps
from ctypes import *
import GBindex
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
# uncomment if running in parallel via Pypar
import pypar
myid = pypar.rank()
nprocs = pypar.size()
node = pypar.get_processor_name()
mlist = GBindex.mlist_100
nlist = GBindex.nlist_100
sigmalist = GBindex.sigmalist_100
anglelist = GBindex.angle_100
num = len(mlist)
ydis = 60
xcell = 3
zcell = 1
grid = 6
overlapinc = 86
overlapboth = 1
natoms_pre = 1

maxiter = 1000000
maxeval = 1000000
ftol = "1.0e-15"
etol = "1.0e-15"
surfelist = []

counter = 0
# if not os.path.exists("Results"): os.makedirs("Results")
if myid == 0:
    mkdir_p("Results")
# -------------------Input EAM pot---------------------------------
pair = "eam/fs"
potfile = "/home/wang/potentials/eam/PotentialA3410-song.fs Fe H"
# potfile="/home/wang/potentials/eam/PotentialB3410.fs Fe H"
# ------------------Input MEAM pot--------------------------------
# pair="meam"
# potfile="/home/wang/potentials/meam/lib-fe-p.meam Fe /home/wang/potentials/meam/Fe-meam-p.meam Fe"
# --------------------Para for lammpsrunning------------------------------
# screen is not good look, rev me when needed
lmp = lammps("", ["-screen", "none"])
lmp.command("variable pair string %s" % (pair))
lmp.command('variable potfile string "%s"' % (potfile))
lmp.command("variable maxiter equal %d" % (maxiter))
lmp.command("variable maxeval equal %d" % (maxeval))
lmp.command("variable ftol equal %s" % (ftol))
lmp.command("variable etol equal %s" % (etol))
# ------------------Minimization of unitslab--------------------------------
lines = open("lammps_scripts/in.unitslab.lammps", 'r').readlines()
for line in lines:
    lmp.command(line)
unite = lmp.extract_variable("unite", 1, 0)
latparam = lmp.extract_variable("latp", 1, 0)
overlapdist = np.arange(0.275, 0.705, 0.005) * latparam
if myid == 0:
    f = open('all-GBs.txt', 'w')
    f.write('#Angle\tGB_energy\tSurf_energy\tCoh_energy\tCount\tX\tZ\n')
    f.close
for i in np.arange(num):
    gbelist = []
    cohelist = []
    counterlist = []
    tolnum = 0
    tolnumlist = []
    # ---------------folder for results---------------------------------
    file_name = "Fe_S%d_(%d_%d_%d)_%dA_STGB" % (
        sigmalist[i], mlist[i], nlist[i], 0, ydis)
    # if not os.path.exists("Results/"+file_name): os.makedirs("Results/"+file_name)
    if myid == 0:
        mkdir_p("Results/" + file_name)
    vari = [
        # ----------------variables defined by python------------------------
        "variable ydis equal %d" % (ydis),
        "variable m equal %d" % (mlist[i]),
        "variable n equal %d" % (nlist[i]),
        "variable sigma equal %d" % (sigmalist[i]),
        "variable angle equal %f" % (anglelist[i]),
        "variable gbname string %s" % (file_name),
        "log log.${gbname}.lammps",
        "variable xcell equal %d" % (xcell),
        "variable zcell equal %d" % (zcell),
        "variable grid equal %d" % (grid),
        # "variable overlapinc equal %d"%(overlapinc),
        "variable overlapboth equal %d" % (overlapboth),
        "variable minimumenergy equal %f" % (unite),
        "variable latparam equal %f" % (latparam),
        # ---------------orientation!----------------------------
        "variable x1 equal ${n}",
        "variable x2 equal -${m}",
        "variable x3 equal 0",
        "variable y1 equal ${m}",
        "variable y2 equal ${n}",
        "variable y3 equal 0",
        "variable x1p equal ${n}",
        "variable x2p equal ${m}",
        "variable x3p equal 0",
        "variable y1p equal -${m}",
        "variable y2p equal ${n}",
        "variable y3p equal 0",
        "variable z1 equal 0",
        "variable z2 equal 0",
        "variable z3 equal 1",
        "variable z1p equal 0",
        "variable z2p equal 0",
        "variable z3p equal 1",
        # -----------------other important paras-----------------------
        'variable ysize equal "sqrt(v_y1^2 + v_y2^2 + v_y3^2)"',
        'variable xsize1 equal "sqrt(v_x1^2 + v_x2^2 + v_x3^2)"',
        'variable zsize1 equal "sqrt(v_z1^2 + v_z2^2 + v_z3^2)"',
        'variable xsize2 equal "sqrt(v_x1^2 + v_x2^2 + v_x3^2)"',
        'variable zsize2 equal "sqrt(v_z1^2 + v_z2^2 + v_z3^2)"',
        'variable xlen equal "v_xsize * v_latparam"',
        'variable ylen equal "v_ysize * v_latparam"',
        'variable zlen equal "v_zsize * v_latparam"',
        'variable ycell equal "ceil(v_ydis/v_ylen)"',
        'if "${xsize1} <= ${xsize2}" then "variable xsize equal ${xsize1}" else "variable xsize equal ${xsize2}"',
        'if "${zsize1} <= ${zsize2}" then "variable zsize equal ${zsize1}" else "variable zsize equal ${zsize2}" '
    ]
    for v in vari:
        lmp.command(v)
    # --------------Surface calculation deck-------------------
    lines = open("lammps_scripts/in.surface.lammps", 'r').readlines()
    for line in lines:
        lmp.command(line)
    surfe = lmp.extract_variable("surfemJm2", "all", 0)
    lmp.command('variable surfemJm2 equal %f' % (surfe))
    surfelist.append(surfe)
    # --------------GB calculation deck--------------------------
    # --------------define grid displacement unit vector-------------
    xlen = lmp.extract_variable("xlen", "all", 0)
    zlen = lmp.extract_variable("zlen", "all", 0)
    xsize = lmp.extract_variable("xsize", "all", 0)
    zsize = lmp.extract_variable("zsize", "all", 0)
    # print xsize
    inc = latparam / grid
    xinc = np.floor(xlen / inc)
    # print xinc
    zinc = np.floor(zlen / inc)
    # print zinc
    if myid == 0:
        f = open('temp_results.txt', 'w')
        f.write('#Angle\tGB_energy\tSurf_energy\tCoh_energy\tCount\tX\tZ\n')
        f.close

    for a in np.arange(xinc):
        # print a
        tx = a / xinc * xsize
        lmp.command('variable tx equal %f' % (tx))
        for b in np.arange(zinc):
            tz = a / zinc * zsize
            lmp.command('variable tz equal %f' % (tz))
            lmp.command(
                'variable d equal %f' %
                (overlapboth))  # not for non-sym
            for c in overlapdist:
                tolnum += 1
                lmp.command('variable overlapdist equal %f' % (c))
                if c == 0.275 * latparam:
                    natoms_pre = 1
                lines = open(
                    "lammps_scripts/in.GB_create.lammps",
                    'r').readlines()
                for line in lines:
                    lmp.command(line)
                # print "read"
                natoms = lmp.get_natoms()
                if natoms != natoms_pre:
                    counter += 1
                    tolnumlist.append(tolnum)
                    natoms_pre = natoms
                    lmp.command('variable counter equal %d' % (counter))
                    lines = open(
                        "lammps_scripts/in.GB_min.lammps",
                        'r').readlines()
                    for line in lines:
                        lmp.command(line)
                        #lmp.command('variable tx delete')
                        #lmp.command('variable tz delete')
                    gbemJm2 = lmp.extract_variable("gbemJm2", "all", 0)
                    gbelist.append(gbemJm2)
                    cohe = lmp.extract_variable("cohe", "all", 0)
                    cohelist.append(cohe)
                    trjname = "%.2f_x%d_z%d_%d.lammpstrj" % (
                        gbemJm2, a, b, tolnum)
                    if myid == 0:
                        print "Caught number %d:" % (counter)
                        print file_name
                        print trjname
                        print "................"
                        rename_p(
                            "temp.lammpstrj",
                            "Results/" +
                            file_name +
                            "/" +
                            trjname)
                        f = open('temp_results.txt', 'a')
                        f.write(
                            '%f\t%f\t%f\t%f\t%d\t%f\t%f\n' %
                            (anglelist[i], gbemJm2, surfe, cohe, tolnum, tx, tz))
                        f.close
    # ---------------Sort results-------------------------------
    if myid == 0:
        print "Total number of model is %d, we took %d" % (tolnum, counter)
        f = open('temp_results.txt', 'r')
        data = np.genfromtxt('temp_results.txt', names=True)
        i = np.argmin(np.abs(data['GB_energy']))
        min_angle = data['Angle'][i]
        min_gb = data['GB_energy'][i]
        min_surf = data['Surf_energy'][i]
        min_coh = data['Coh_energy'][i]
        min_cont = data['Count'][i]
        min_x = data['X'][i]
        min_z = data['Z'][i]
        f.close
        print 'The ID of optimized GB is %d\nLowest GB energy is %f, Cohensive energy is %f' % (min_cont, min_gb, min_coh)
        f = open('all-GBs.txt', 'a')
        f.write(
            '%f\t%f\t%f\t%f\t%d\t%f\t%f\n' %
            (min_angle,
             min_gb,
             min_surf,
             min_coh,
             min_cont,
             min_x,
             min_z))
        f.close
        mintrjname = glob.glob(
            "Results/" +
            file_name +
            "/" +
            "*_%d.lammpstrj" %
            (min_cont))
                # print mintrjname
        shutil.copy2(mintrjname[0], "Results/" + file_name + ".lammpstrj")
        newlog_GB = "Results/" + file_name + "/log." + file_name + ".lammps"
        newsum_GB = "Results/" + file_name + "/sum." + file_name + ".txt"
        rename_p("log." + file_name + ".lammps", newlog_GB)
        rename_p("temp_results.txt", newsum_GB)
        #         f = open('all-GBs.txt', 'r')
        #         data = np.genfromtxt('all-GBs.txt',names = True)
        #         fig=plt.figure()
        #         plt.ylabel('$\mathrm{Grain\ boundary\ energy\ (J/m^2)}$')
        #         plt.xlabel('$\mathrm{Angle\ (degrees)}$')
        # plt.axis([0,90,0,1500])
        #         plt.grid()
        #         plt.plot(data['Angle'],data['GB_energy'],'--ro')
        # plt.show()
        #         fig.savefig('allGB.png', transparent=True, bbox_inches='tight', pad_inches=0, dpi=300)
    counter = 0
    tolnum = 0

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
