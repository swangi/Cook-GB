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



def x0_reader(data_file):
    data = np.reshape(np.genfromtxt(data_file), (-1, 4))
    x = np.linspace(data[0][0], data[0][2], num)
    return x


def line_reader(x_0, data_file):
    data = np.reshape(np.genfromtxt(data_file), (-1, 4))
    y_array = []
    k_array = []
    x_array = []
    for d in data:
        x1, y1, x2, y2 = d
        k = (y2 - y1) / (x2 - x1)
        y = [k * (x - x1) + y1 for x in x_0]
        y_array.append(y)
        x_array.append(x_0)
        k_array.append(k)
    return np.array(x_array), np.array(y_array),  np.array(k)


def main():
    data_files = ['excess_gb_H.txt', 'excess_fs_H.txt']
    # data GB
    x0 = x0_reader(data_files[0])
    x_array_GB, y_array_GB, k_array_GB = line_reader(x0, 'excess_gb_H.txt')
    x_array_GB = x_array_GB[0:45]
    y_array_GB = y_array_GB[0:len(x_array_GB)]
    y_min_GB = np.array([min(y) for y in y_array_GB.T])
    # data FS
    x_array_FS, y_array_FS, k_array_FS = line_reader(x0, 'excess_fs_H.txt')
    print len(x_array_GB)
    x_array_FS, y_array_FS = (x_array_FS[0:len(x_array_GB) / 2],
                              y_array_FS[0:len(y_array_GB) / 2])
    y_min_FS = np.array([min(y) for y in y_array_FS.T])
    # plot
    for i, x in enumerate(x_array_GB):
        ax[0].plot(x, y_array_GB[i], 'g--', alpha=0.2)
    for i, x in enumerate(x_array_FS):
        ax[0].plot(x, y_array_FS[i], 'g.-', alpha=0.2)
    ax[0].plot(x0, y_min_GB, '-', lw=5)
    ax[0].plot(x0, y_min_FS, '-', lw=5)
    ax[1].plot(x_array_GB[0], y_min_FS -
               y_min_GB,
               'o',
               ms=8)

    # set plot and save
    ax[0].set_xlabel(r'Hydrogen chemical potential (eV)')
    ax[0].set_ylabel(r'Interface energy (mJ/m$^2$)')
    ax[1].set_xlabel(r'Hydrogen chemical potential (eV)')
    ax[1].set_ylabel(r'Grain boundary cohesive energy (mJ/m$^2$)')
    for i, x in enumerate(ax):
        if i == 0:
            label_lc = 'a'
        else:
            label_lc = 'b'
        x.set_xlim(-0.22, 0.25)
        x.text(0.013, 0.981, label_lc,
               ha='left',
               va='top',
               fontsize=24,
               bbox=bbox_props,
               transform=x.transAxes
               )
    # x.set_ylim(-1000, 1000)
    fig.set_size_inches(11, 16)
    fig.savefig('summary.pdf',
                # transparent=True,
                bbox_inches='tight',
                pad_inches=0,
                dpi=300)
    # plt.show()
if __name__ == "__main__":
    num = 100
    bbox_props = dict(boxstyle="square", fc="w", ec="w", alpha=1.0)
    fig, ax = plt.subplots(2, 1)
    main()
