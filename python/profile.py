#!/usr/bin/env python

import matplotlib.pyplot as plt
import sys
import utils

# get fields
snapshots = []
snapshotNames = []
with open(sys.argv[1], "r") as f:
    lines = f.readlines()
    for i in range(len(lines)):
        l = lines[i]
        snapshot = l.split()
        snapshotData = {}
        for i in range(int(snapshot[1])):
            statsName = snapshot[i * 2 + 2]
            statsValue = float(snapshot[i * 2 + 3])
            snapshotData[statsName] = statsValue
        snapshotNames.append(snapshot[0])
        snapshots.append(snapshotData)
indices = {}
statsNum = 0
for i in snapshots[-1]:
    indices[i] = statsNum
    statsNum += 1
# setup data
data = []
for i in range(statsNum):
    data.append(len(snapshotNames) * [0])

statSum = statsNum * [0]
for stat in indices:
    for i in range(len(snapshotNames)):
        if stat in snapshots[i]:
            data[indices[stat]][i] = snapshots[i][stat]
        statSum[indices[stat]] += data[indices[stat]][i]

fig, ax = plt.subplots(figsize=(16, 9))
for stat in indices:
    if 'Eigen' in stat:
        ax.plot(range(len(data[indices[stat]])), data[indices[stat]], label=stat)
ax.set_xlabel('Frame')
ax.set_ylabel('Value')
ax.legend()
plt.show()

# eigen count bars
fig, ax = plt.subplots(figsize=(16, 9))
numbers = []
for stat in indices:
    if 'Eigen' in stat:
        numbers.append(int(stat.split('Solve')[1]))
numbers.sort()
d = []
for n in numbers:
        d.append(statSum[indices['EigenSolve' + str(n)]])
plt.bar(numbers, d)
plt.ylabel('# of systems')
plt.xlabel('System Size')
plt.title('Count of Eigen system solves')
plt.show()
