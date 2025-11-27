# python db4.py 4 6 2 40 10 25
import sys

if len(sys.argv)  != 7:
    print("Usage (all integers): min_level max_level frames domain_size box_size box_height")
    exit(1)

min_level = int(sys.argv[1])
max_level = int(sys.argv[2])
frames = int(sys.argv[3])
# in grid coordinates based on max level:
domain_size = int(sys.argv[4])
box_size = int(sys.argv[5])
box_height = int(sys.argv[6])

full_size = int(2**max_level)

print("minLevel " + str(min_level))
print("maxLevel " + str(max_level))
print("outputDir db4")
print("frames " + str(frames))
print("dt 0.005")
print("particleSpacing", 1.0 / (full_size * 2))
print("scene NONE")
t = 3 * [max_level]
s = "iSolidBoxes " + ' '.join([str(x) for x in t]) + ' '
s += ' '.join(['0 0 0', str(full_size), "1", str(full_size)]) + ' '
s += ' '.join(['0 0 0 1', str(full_size), str(domain_size)]) + ' '
s += ' '.join([str(domain_size - 1), '0 0', str(domain_size), str(full_size), str(domain_size)]) + ' '
s += ' '.join(['0 0 0', str(domain_size), str(full_size),'1']) + ' '
s += ' '.join(['0 0', str(domain_size - 1), str(full_size), str(full_size), str(domain_size)])
print(s)
s = "iFluidBoxes " + ' '.join([str(x) for x in t]) + ' '
s += ' '.join(["1 1 1", str(box_size), str(box_height), str(box_size), '1']) + ' '
s += ' '.join([str(domain_size - box_size - 1), '1 1', str(domain_size-1), str(box_height), str(box_size), '2']) + ' '
s += ' '.join([str(domain_size - box_size - 1), '1', str(domain_size - box_size - 1), str(domain_size-1), str(box_height), str(domain_size-1), '3']) + ' '
s += ' '.join(['1 1', str(domain_size - box_size - 1), str(box_size), str(box_height), str(domain_size-1), "4"])
print(s)
print("rbfKernel PHS")
print("base QUADRATIC")
print("graded 1")
print("algorithm advectParticles updateTree markSolids markBoundariesOnFaces particlesToFaces  addGravityToFaces divergenceOnCellsFromFaces  pressureOnCellsRing pressureGradientOnFacesFromCells correctFaceVelocitiesFromFacePressureGradient reseedParticles facesToParticlesAllComponents")
