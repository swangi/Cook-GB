#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (c) 2014 - Swang <swangi@outlook.com>
# Filename:
# +-+-+-+-+-+-+ +-+-+-+-+-+-+-+-+-+
# |I|m|p|o|r|t| |s|o|m|e|t|h|i|n|g|
# +-+-+-+-+-+-+ +-+-+-+-+-+-+-+-+-+
#
from pylab import *
import re
#
# +-+-+-+-+-+ +-+-+ +-+-+-+-+
# |S|t|a|r|t| |o|f| |m|a|i|n|
# +-+-+-+-+-+ +-+-+ +-+-+-+-+
#
lines_needed = []
if __name__ == "__main__":
    with open('GB-100-Linescan.o96', 'r') as f:
        for line in f:
            if re.match("(.*)FS\s\\((.*)", line):
                # temp_line = line.split()[7:9]
                # temp_line[2].split('(')
                # ready_line = temp_line.split('(')
                # print temp_line[1].split('(')
                type_fs = line.split()[7]
                eng_fs = line.split()[8].split('(')[1]
                if [type_fs, eng_fs] not in lines_needed:
                    lines_needed.append([type_fs, eng_fs])
    with open('fs_sum.csv', 'w') as fs_sum:
        for ls in lines_needed:
            fs_sum.write('%s,%s\n' % (ls[0], ls[1]))
    # print lines_needed
