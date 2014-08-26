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
from lammps import lammps
#
# +-+-+-+-+-+ +-+-+ +-+-+-+-+
# |S|t|a|r|t| |o|f| |m|a|i|n|
# +-+-+-+-+-+ +-+-+ +-+-+-+-+
#


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


def GB_eng():
    global gbe_0
    global gbemJm2_0
    global area_gb
    global pote_gb
    read_lmp_line("min_GB.lammps")
    num_fe = lmp.get_natoms()
    pote_gb = lmp.extract_variable("pot", "all", 0)
    area_gb = lmp.extract_variable("area", "all", 0)
    gbe_0 = (pote_gb - (fe_bulk * num_fe)
             ) / area_gb / 2
    gbemJm2_0 = gbe_0 * 16021.7733
    # lmp.command("write_data %sx%d_GB.lammpstrj" % (GB_type, nx))


def FS_eng():
    global fse_0
    global fsemJm2_0
    global area_fs
    global pote_fs
    # read_lmp_line("FS_create.lammps")
    read_lmp_line("min.lammps")
    pote_fs = lmp.extract_variable("pot", "all", 0)
    area_fs = lmp.extract_variable("area", "all", 0)
    fse_0 = (pote_fs - (fe_bulk * num_fe)
             ) / area_fs - gbe_0
    fsemJm2_0 = fse_0 * 16021.7733
    # lmp.command("write_data %sx%d_FS.lammpstrj" % (GB_type, nx))


def interface_H_min(pos_h, gamma_int, data_pre):
    for pos in pos_h:
        lmp.command("create_atoms 2 single %5f %5f %5f units box" %
                    (pos[0], pos[1], pos[2]))
    num_h = lmp.get_natoms() - num_fe
    print "number of H is %d" % (num_h)
    read_lmp_line("min.lammps")
    pote_h_face = lmp.extract_variable("pot", "all", 0)
    ed = pote_h_face - gamma_int - h_coh * num_h
    area_face_H = lmp.extract_variable("area", "all", 0)
    data_name = f + '/Summary/Results-excess/' + "%s_%d_H.lammpstrj" % (
        data_pre, num_h)
    lmp.command("write_data %s" %
                (data_name))
    # if myid == 0:
    #     rename_p(data_name, f + '/Summary/' + "Results-excess/%s" %
    #              (data_name))
    return num_h, ed, area_face_H, pote_h_face


def rich_poor_eng(data_prefix, posH_eff):
    if data_prefix == 'GB':
        # +-+-+ +-+-+-+-+-+-+-+-+-+-+-+
        # |G|B| |c|a|l|c|u|l|a|t|i|o|n|
        # +-+-+ +-+-+-+-+-+-+-+-+-+-+-+
        muh_poor_gb = min(data['EdH']) + h_coh
        muh_rich_gb = -2.087830167184861
        muh_poor_r_gb = muh_poor_gb - h_coh
        muh_rich_r_gb = muh_rich_gb - h_coh
        print 'On GB\nMin[Ed] is %.5f, Max[Ed] is %.5f ' % (
            muh_poor_r_gb, muh_rich_r_gb)
        lammps_initial(gb_init)
        GB_eng()
        # add H

        num_h, ed_h_gb, h_area_gb, pote_h_gb = interface_H_min(
            posH_eff,
            pote_gb,
            type_name + "GB")
        # print num_h
        # Ed_H_GB.append(ed_h_gb)
        # H_area_GB.append(h_area_gb)
        # Pote_h_GB.append(pote_h_gb)
        gbe_h_poor = 1000 * ((pote_h_gb - (fe_bulk * num_fe)
                              - (muh_poor_gb * num_h)
                              ) / h_area_gb - gbe_0)
        gbe_h_rich = 1000 * ((pote_h_gb - (fe_bulk * num_fe)
                              - (muh_rich_gb * num_h)
                              ) / h_area_gb - gbe_0)
        gbemJm2_h_poor = gbe_h_poor * 16021.7733 / 1000
        gbemJm2_h_rich = gbe_h_rich * 16021.7733 / 1000
        ecess_h = num_h / h_area_gb * 100
        if myid == 0:
            print "##########################################"
            print('GB energy without H is %.5f mJm2\
            \nPotential energy of GB with (%.5f nm^-2) H is %.5f eV \
            \nGB energy (poor H) = %.5f mJm2 (or %.5f meV/Ang2)\
            \nGB energy (rich H) = %.5f mJm2 (or %.5f meV/Ang2)\
            \nDissolution E of H(GB) is %.5f eV') % (
                gbemJm2_0, ecess_h, pote_h_gb,
                gbemJm2_h_poor, gbe_h_poor,
                gbemJm2_h_rich, gbe_h_rich, ed_h_gb)
            print "##########################################"
            fgb = open(record_file_gb, 'a')
            fgb.write('%.5f,%.5f,%.5f,%.5f,%.5f,%.5f\n' %
                     (ecess_h, muh_poor_r_gb,
                      gbemJm2_h_poor, gbe_h_poor,
                      area_gb, h_area_gb))
            fgb.write('%.5f,%.5f,%.5f,%.5f,%.5f,%.5f\n' %
                     (ecess_h, muh_rich_r_gb,
                      gbemJm2_h_rich, gbe_h_rich,
                      area_gb, h_area_gb))
            fgb.close

    if data_prefix == 'FS':
        # +-+-+ +-+-+-+-+-+-+-+-+-+-+-+
        # |F|S| |c|a|l|c|u|l|a|t|i|o|n|
        # +-+-+ +-+-+-+-+-+-+-+-+-+-+-+
        muh_poor_fs = min(data['EdH']) + h_coh
        muh_rich_fs = -2.087830167184861
        muh_poor_r_fs = muh_poor_fs - h_coh
        muh_rich_r_fs = muh_rich_fs - h_coh
        print 'On GB\nMin[Ed] is %.5f, Max[Ed] is %.5f ' % (
            muh_poor_r_fs, muh_rich_r_fs)
        lammps_initial(fs_init)
        FS_eng()
        num_h, ed_h_fs, h_area_fs, pote_h_fs = interface_H_min(
            posH_eff,
            pote_fs,
            type_name + "FS")
        # Ed_H_FS.append(ed_h_fs)
        # H_area_FS.append(h_area_fs)
        # Pote_h_FS.append(pote_h_fs)
        fse_h_poor = 1000 * ((pote_h_fs - (fe_bulk * num_fe)
                              - (muh_poor_fs * num_h)
                              ) / h_area_fs - gbe_0)
        fse_h_rich = 1000 * ((pote_h_fs - (fe_bulk * num_fe)
                              - (muh_rich_fs * num_h)
                              ) / h_area_fs - gbe_0)
        fsemJm2_h_p = fse_h_poor * 16021.7733 / 1000
        fsemJm2_h_r = fse_h_rich * 16021.7733 / 1000
        # Assuming up and two have the same H trapped
        fsemJm2_h_poor = 2 * (fsemJm2_h_p - fsemJm2_0 / 2)
        fsemJm2_h_rich = 2 * (fsemJm2_h_r - fsemJm2_0 / 2)
        ecess_h = num_h / h_area_fs * 100
        if myid == 0:
            print "##########################################"
            print('Surface energy without H is %.5f mJm2\
                \n Potential energy of FS with (%.5f nm^-2) H is %.5f eV \
                \nFS energy (poor H) = %.5f mJm2 (or %.5f meV/Ang2)\
                \nFS energy (rich H) = %.5f mJm2 (or %.5f meV/Ang2)\
                \nDissolution E of H(FS) is %.5f eV') % (
                fsemJm2_0, ecess_h, pote_h_fs,
                fsemJm2_h_poor, fse_h_poor,
                fsemJm2_h_rich, fse_h_rich, ed_h_fs)
            print "##########################################"
            fs = open(record_file_fs, 'a')
            fs.write('%.5f,%.5f,%.5f,%.5f,%.5f,%.5f\n' %
                     (ecess_h, muh_poor_r_fs,
                      fsemJm2_h_poor, fse_h_poor,
                      area_fs, h_area_fs))
            fs.write('%.5f,%.5f,%.5f,%.5f,%.5f,%.5f\n' %
                     (ecess_h, muh_rich_r_fs,
                      fsemJm2_h_rich, fse_h_rich,
                      area_fs, h_area_fs))
            fs.close
# +-+-+-+ +-+-+-+-+ +-+-+-+-+-+
# |r|u|n| |w|i|t|h| |p|y|p|a|r|
# +-+-+-+ +-+-+-+-+ +-+-+-+-+-+
import pypar
myid = pypar.rank()
nprocs = pypar.size()
node = pypar.get_processor_name()
# lammps_scr_path = "/home/swang/share/lammps_scripts/"
# pot_path = "/home/swang/share/potentials/eam/Ram_song_potA.fs"
lammps_scr_path = "/Users/swang/Simulation/LAMMPS/Py/lammps_scripts/"
pot_path = "/Users/swang/Simulation/potentials/EAM/Ram_song_potA.fs"
maxiter = 1000000
maxeval = 1000000
ftol = "1.0e-5"
etol = "0.0"
surfelist = []
msd = []
counter = 0
lmp = lammps("", ["-screen", "none"])
pair = "eam/fs"
potfile = "%s Fe H" % (pot_path)
T = 300
T2 = T * 2
dump_id = "mydump"
dump_step = 1
dump_name = "mydump"
nx = 1
ny = 1
nz = 1
h_coh = -2.36800333126474
fe_bulk = -4.155410449999999755
# dir_names = [os.path.basename(x)
#              for x in glob.glob('*Results')]
dir_names = ['1_2_0_Results']
type_names = ['_'.join(d.split('_')[0:3])
              for d in dir_names]  # grain boundary type with a form like 1_2_0
angle_list = []

for t in type_names:
    m = int(t.split('_')[0])
    n = int(t.split('_')[1])
    angle_list.append(2 * np.degrees(np.arctan2(m, n)))
    # print angle_list

if __name__ == "__main__":
    # lammps_initial(initial_trj)
    for i, f in enumerate(dir_names):
        gb_init = f + '/' + type_names[i] + '_GB.lammpstrj'
        fs_init = f + '/' + type_names[i] + '_FS.lammpstrj'
        record_dir = f + '/Summary/' + "Results-excess"
        record_file_gb = record_dir + '/excess_gb_H.csv'
        record_file_fs = record_dir + '/excess_fs_H.csv'
        if myid == 0:
            mkdir_p(record_dir)
            with open(record_file_gb, 'w') as figb:
                figb.write(
                    "#ecess_h,muh_r_gb,gbemJm2_h_poor,gbe_h,area_gb,area_gb_h\n")

            with open(record_file_fs, 'w') as fs:
                fs.write(
                    "#ecess_h,muh_r_fs,fsemJm2_h_poor,fse_h,area_fs,area_fs_h\n")

        type_name = type_names[i]
        angle = angle_list[i]
        file_names = [os.path.basename(x)
                      for x in glob.glob(f + '/Summary/*_EdH_sum.txt')]
        file_names = [file_names[-1], file_names[0]]
        for d in file_names:
            posH_eff = []
            data_name = d.split('.')[0]
            data_prefix = data_name.split('_')[3]  # what type of interface ?
            # print data_prefix
            data = np.genfromtxt(
                f + '/Summary/' + d, names=True, delimiter=",")
            data = data[np.argsort(data['EdH'])]
            # print data
            min_trap_id = np.where(data['EdH'] == data['EdH'].min())
            ok_trap_id = np.where(data['EdH'] == data['EdH'].min())
            low_eng_id = np.where(
                data['EdH'] - data['EdH'].min() <= 0.01)
            low_eng = data['EdH'][low_eng_id[0][0]]
            up_eng_id = np.where(data['EdH'] == data['EdH'].max())
            up_eng = data['EdH'][up_eng_id[0][0]]
            # Ed less than t site is OK trapping
            OK_eng_id = np.where(data['EdH'] <= 0.20)
            # print data[OK_eng_id]
            # Ed over o site is NG trapping
            NG_eng_id = np.where(data['EdH'] >= 0.33)
            np.sort(data, order='EdH')
            # transpose position x y z to a single list
            posH = zip(
                *[data['x'][OK_eng_id],
                  data['y'][OK_eng_id],
                  data['z'][OK_eng_id]])
            # posH = zip(
            #     *[data['x'][NG_eng_id],
            #       data['y'][NG_eng_id],
            #       data['z'][NG_eng_id]])
            # start of lammps
            if d == file_names[0]:
                np.savetxt(
                    'pos_H_excess_record_GB.txt',
                    map(np.array, posH),
                    fmt="%.5f", delimiter=',')
            for p in posH:
                posH_eff.append(p)
                rich_poor_eng(data_prefix, posH_eff)


    # if myid == 0:
    #     print '''
    #              +-+-+-+-+-+-+-+-+
    #              |F|i|n|i|s|h|e|d|
    #              +-+-+-+-+-+-+-+-+
    #      +-+-+-+-+ +-+ +-+-+-+-+ +-+-+-+-+
    #      |H|A|V|E| |A| |N|I|C|E| |D|A|Y|!|
    #      +-+-+-+-+ +-+ +-+-+-+-+ +-+-+-+-+'''
    # 	pypar.finalize()
