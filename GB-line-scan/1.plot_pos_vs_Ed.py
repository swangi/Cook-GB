#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (c) 2014 - Swang <swangi@outlook.com>
# Filename: plot_pos_vs_Ed.py
# +-+-+-+-+-+-+ +-+-+-+-+-+-+-+-+-+
# |I|m|p|o|r|t| |s|o|m|e|t|h|i|n|g|
# +-+-+-+-+-+-+ +-+-+-+-+-+-+-+-+-+
#
from pylab import *
import os
import glob
from scipy import spatial
# from sets import Set
#
# +-+-+-+-+-+ +-+-+ +-+-+-+-+
# |S|t|a|r|t| |o|f| |m|a|i|n|
# +-+-+-+-+-+ +-+-+ +-+-+-+-+
#


def fliter_dup(data_copy, tol=0.1):
    data_tmp = map(list, zip(data_copy['x'], data_copy['y'], data_copy['z']))
    tree = spatial.KDTree(data_tmp)
    to_del = np.array(map(list, tree.query_pairs(tol)))[:, 0]
    dup_list = np.unique(to_del)
    return dup_list
    # print len(np.unique(del_list))
    # data_new = np.delete(data_tmp, del_list, 0)
    # print len(data_new)


def plot_fig(x, pp_prefix, y, data_prefix, obd):
    fig = plt.figure()
    plt.ylabel('$\mathrm{Dissolution\ Energy\ (eV)}$')
    plt.xlabel('$\mathrm{Distance\ ' + '(\AA)}$')
    plt.plot(x, y, 'ro')
    plt.grid()
    fig.savefig(obd + '/' + data_name + '_Ed_vs_pos_' + pp_prefix + '.png',
                transparent=True,
                bbox_inches='tight',
                pad_inches=0,
                dpi=300)
    plt.close(fig)


def mkdir_p(path):
    try:
        os.makedirs(path)
    except os.error as e:
        if e.errno != 17:
            raise
        else:
            pass


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
        data_A.append(map(float, i.split(',')[0:8]))
    for i in linesb[16:num_B + 16]:
        data_B.append(map(float, i.split()[0:8]))
    num_A = len(data_A)
    num_all = num_A + num_B
    output_file = file_B + '_H_linescan_sum.lammpstrj'
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
    # rename_p(output_file, '%s/%s' % (Result_dir, output_file))


def rename_p(old, new):
    try:
        os.rename(old, new) 
    except os.error as e:
        if e.errno != 2:
            raise
        else:
            pass
dir_names = [os.path.basename(x)
             for x in glob.glob('*Results')]
# dir_names = ['2_13_0_Results', '2_15_0_Results']
type_names = ['_'.join(d.split('_')[0:3]) for d in dir_names]
# print type_names
angle_list = []
for t in type_names:
    m = int(t.split('_')[0])
    n = int(t.split('_')[1])
    angle_list.append(2 * np.degrees(np.arctan2(m, n)))
    # print angle_list
with open('summary_FS.csv', 'w') as s:
    s.write('Type,Angle,min_Ed,max_Ed,Nums_min,OK_trap,NG_trap,Area,trap_distance,trap_density(A^2)\n')
with open('summary_GB.csv', 'w') as s:
    s.write('Type,Angle,min_Ed,max_Ed,Nums_min,OK_trap,NG_trap,Area,trap_distance,trap_density(A^2)\n')
if __name__ == "__main__":
    for i, f in enumerate(dir_names):
        type_name = type_names[i]
        angle = angle_list[i]
        print f
        mkdir_p(f + '/verbose')
        verbose_name = [os.path.basename(x)
                        for x in glob.glob(f + '/*z.lammpstrj')]
        # print verbose_name
        for v in verbose_name:
            # print v
            rename_p(f + '/' + v, f + '/verbose/' + v)
        file_names = [os.path.basename(x)
                      for x in glob.glob(f + '/*_Dissolution*.txt')]
        for d in file_names:
            data_name = d.split('.')[0]
            data_prefix = data_name.split('_')[-1]
            print d
            data = np.genfromtxt(f + '/' + d, names=True)
            np.sort(data, order='EdH')
            del_list = fliter_dup(data)
            data = np.delete(data, del_list, 0)
            low_eng_id = np.where(
                data['EdH'] - data['EdH'].min() <= 0.01)
            low_eng = data['EdH'][low_eng_id[0][0]]
            up_eng_id = np.where(data['EdH'] == data['EdH'].max())
            up_eng = data['EdH'][up_eng_id[0][0]]
            # Ed less than t site is OK trapping
            OK_eng_id = np.where(data['EdH'] <= 0.20)
            # Ed over o site is NG trapping
            NG_eng_id = np.where(data['EdH'] >= 0.33)
            # get the prefered trapping range
            prefered_trap_pos = [(data['y'][x]) for x in OK_eng_id[0]]
            # up limit of traping in y direction
            tr_up = max(prefered_trap_pos)
            # low limit of traping in y direction
            tr_low = min(prefered_trap_pos)
            tr_distance = tr_up - tr_low
            print 'effective traping distance is %.5f' % (tr_distance)
            print 'Summarying................'
            file_B = f + '/%s_GB.lammpstrj' % (type_name)
            with open(file_B, 'r') as fb:
                linesb = fb.readlines()
            num_B = int(linesb[2].split()[0])
            summary_FS = f + '/%s_FS_EdH_sum.txt' % (type_name)
            summary_GB = f + '/%s_GB_EdH_sum.txt' % (type_name)
            # uncomment below if modifying of grain boundary position is needed
            # pos_cor = data['y'][low_eng_id[0][0]]
            # data['y'][:] = data['y'][:] - pos_cor
            # plt.axis([0, 90, 0, 1250])
            # +-+-+-+-+ +-+-+-+-+-+-+-+
            # |p|l|o|t| |f|i|g|u|r|e|s|
            # +-+-+-+-+ +-+-+-+-+-+-+-+
            plot_fig(data['x'], 'x', data['EdH'], data_prefix, f)
            plot_fig(data['y'][:], 'y', data['EdH'], data_prefix, f)
            plot_fig(data['z'], 'z', data['EdH'], data_prefix, f)

            if data_prefix == 'FS':
                with open(summary_FS, 'w') as fs:
                    fs.write('id,type,x,y,z,EdH,area,layer_n\n')
                trap_density = len(
                    OK_eng_id[0]) / data['area'][0]  # / tr_distance
                with open('summary_FS.csv', 'a') as s:
                    s.write(
                        type_name + '_' + data_prefix + ',' +
                        '%.5f,%.5f,%.5f,%d,%d,%d,%.5f,%.5f,%.5f\n' %
                        (angle, low_eng, up_eng,
                            len(low_eng_id[0]),
                            len(OK_eng_id[0]),
                            len(NG_eng_id[0]),
                            data['area'][0],
                            tr_distance,
                            trap_density))
                n_H = 0
                with open(summary_FS, 'a') as fs:
                    for da in data:
                        n_H += 1
                        id_H = num_B + n_H
                        fs.write('%d,%d,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f\n' %
                                 (id_H, 2, da[2], da[3],
                                  da[4], da[5], da[6], da[7]))
                make_traj(f + '/%s_FS.lammpstrj' % (type_name), summary_FS)
            if data_prefix == 'GB':
                with open(summary_GB, 'w') as gb:
                    gb.write('id,type,x,y,z,EdH,area,layer_n\n')
                trap_density = len(
                    OK_eng_id[0]) / data['area'][0]  # / tr_distance
                with open('summary_GB.csv', 'a') as s:
                    s.write(
                        type_name + '_' + data_prefix + ',' +
                        '%.5f,%.5f,%.5f,%d,%d,%d,%.5f,%.5f,%.5f\n' %
                        (angle, low_eng, up_eng,
                         len(low_eng_id[0]),
                         len(OK_eng_id[0]),
                         len(NG_eng_id[0]),
                         data['area'][0],
                         tr_distance,
                         trap_density))
                n_H = 0
                with open(summary_GB, 'a') as fs:
                    for da in data:
                        n_H += 1
                        id_H = num_B + n_H
                        fs.write('%d,%d,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f\n' %
                                 (id_H, 2, da[2], da[3],
                                  da[4], da[5], da[6], da[7]))
                make_traj(f + '/%s_GB.lammpstrj' % (type_name), summary_GB)
            mkdir_p(f + '/Summary')
            file_sums = [os.path.basename(x)
                         for x in glob.glob(f + '/*sum*')]
            for sf in file_sums:
                rename_p(f + '/' + sf, f + '/Summary/' + sf)
