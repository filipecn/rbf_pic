#!/usr/bin/env python

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import sys
import utils


def render_frame(num, data, line):
    xyz_file = open(files[num])

    _x = []
    _y = []
    lines = xyz_file.readlines()
    for i in range(1, len(lines)):
        xyz = lines[i]
        _x.append(float(xyz.split()[0]))
        _y.append(float(xyz.split()[1]))

    line.set_data(_x, _y)
    xyz_file.close()

    return line,


# Simulation Path
path = sys.argv[1]
output_name = sys.argv[2] if len(sys.argv) > 2 else 'result.mp4'
files = utils.get_all_files(path, "*.pv")
files.sort()
size_files = len(files)

# Set up formatting for the movie files
Writer = animation.writers['ffmpeg']
writer = Writer(fps=30, bitrate=1800)

fig1 = plt.figure(figsize=(14, 14))

data = []
l, = plt.plot([], [], 'b.')

plt.xlim(0, 1)
plt.ylim(0, 1)
plt.xlabel('x')
plt.title('Simulation Result')

line_ani = animation.FuncAnimation(fig1, render_frame, size_files,
                                   fargs=(data, l),
                                   interval=1, blit=True)

line_ani.save(output_name, writer=writer, dpi=72)
