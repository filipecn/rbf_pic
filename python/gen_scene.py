import sys
import math


scene = sys.argv[1]
if scene == "well":
    res, h, h_, H, H_, w_, n = int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5]), int(sys.argv[6]), int(sys.argv[7]), int(sys.argv[8])
    Y = int(math.pow(2, res))
    w = 2
    s = 'iSolidBoxes ' + 3 * (str(res) + ' ')
    walls = []
    walls.append([0,H_ + w,0,w,Y,H + 2 * w])
    walls.append([0,H_ + w,0,3 * h_ + 4 * w, Y, w])
    walls.append([2 * h_ + 2 * w, H_ + w, H + w - h, 2 * h_ + 3 * w, Y, H + 2 * w])
    walls.append([h_ + w, H_ + w, 0, h_ + 2 * w, Y, h + w])
    walls.append([3 * h_ + 3 * w, H_ + w, 0, 3 * h_ + 4 * w, Y, H + 2 * w])
    walls.append([0, H_ + w, H + w, 3 * h_ + 4 * w, Y, H + 2 * w])
    walls.append([0,H_ + w,0, 2 * h_ + 3 * w, H_ + w + w, H + 2 * w])
    offset = H - h // 2
    for i in range(1, n):
        walls.append([2 * h_ + 2*2, H_ + w, offset - n + i, 3 * h_ + 4*2, H_ + w + i + w, offset])
    walls.append([2 * h_ + 2 * w, 0, offset - w, 3 * h_ + 4 * w, H_ + 2 * w, offset])
    walls.append([2 * h_ + 2 * w, 0, offset + w_, 3 * h_ + 4 * w, H_ + 2 * w, offset + w_ + w])
    walls.append([2 * h_ + 2 * w, 0, offset - w, 3 * h_ + 4 * w, w, offset + w_ + w])
    walls.append([2 * h_ + 2 * w, 0, offset - w, 2 * h_ + 3 * w, H_ + 2 * w, offset + w_ + w  ])
    walls.append([3 * h_ + 3 * w, 0, offset - w, 3 * h_ + 4 * w, H_ + 2 * w, offset + w_ + w ])
    walls.append([2 * h_ + 2 * w, H_ + w,0, 3 * h_ + 4 * w, H_ + w + w, offset])
    walls.append([2 * h_ + 2 * w, H_ + w,offset + w_, 3 * h_ + 4 * w, H_ + w + w, H + 2 * w])
    
    for w in walls:
        s += ' '.join([str(x) for x in w]) + ' '
    print(s)

elif scene == "kimaze":
    res, h, h_, H, H_, n = int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5]), int(sys.argv[6]), int(sys.argv[7])
    Y = int(math.pow(2, res))
    s = 'iSolidBoxes ' + 3 * (str(res) + ' ')
    w = 2
    walls = []
    walls.append([0,w,w,w,Y,H + w])
    walls.append([0,w,0,3 * h_ + 4*w, Y, w])
    walls.append([2 * h_ + 2*w, w, H + w - h, 2 * h_ + 3*w, Y, H + w])
    walls.append([h_ + w, w, w, h_ + 2*w, Y, h + w])
    walls.append([3 * h_ + 3*w, w, w, 3 * h_ + 4*w, Y, H + w])
    walls.append([0, w, H + w, 3 * h_ + 4*w, Y, H + 2*w])
    walls.append([0,0,0, 3 * h_ + 4*w, w, H + 2*w])
    offset = H - h // 2
    for i in range(1, n):
        walls.append([2 * h_ + 3*w, w, offset - n + i + w, 3 * h_ + 3*w, i + w, offset - n + i - 1 + w])
    for w in walls:
        s += ' '.join([str(x) for x in w]) + ' '
    print(s)
elif scene == "maze":
    res, h, h_, H, H_, n = int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5]), int(sys.argv[6]), int(sys.argv[7])
    Y = int(math.pow(2, res))
    w = 2
    s = 'iSolidBoxes ' + 3 * (str(res) + ' ')
    walls = []
    walls.append([0,0,0,w,Y,H + 2 * w])
    walls.append([0,0,0,3 * h_ + 4 * w, Y, w])
    walls.append([2 * h_ + 2 * w, 0, H + w - h, 2 * h_ + 3 * w, Y, H + 2 * w])
    walls.append([h_ + w, 0, 0, h_ + 2 * w, Y, h + w])
    #            walls.append([2 * h_ + 2, 0, H_ + H + 2, 3 * h_ + 4, Y, H_ + H + 3])
    walls.append([3 * h_ + 3 * w, 0, 0, 3 * h_ + 4 * w, Y, H + 2 * w])
    walls.append([0, 0, H + w, 3 * h_ + 4 * w, Y, H + 2 * w])
    walls.append([0,0,0, 3 * h_ + 4 * w, w, H + 2 * w])
    offset = H - h // 2
    for i in range(1, n):
        walls.append([2 * h_ + 2*2, 0, offset - n + i, 3 * h_ + 4*2, i + w, offset])
    for w in walls:
        s += ' '.join([str(x) for x in w]) + ' '
    print(s)
elif scene == "drop":
    # h - tank size m = tank height 
    res, h, m = int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4])
    maxCoord = int(math.pow(2, res))
    s = 'iSolidBoxes ' + 3 * (str(res) + ' ')
    floor = [0, 0, 0, h, 1, h]
    s += ' '.join([str(x) for x in floor]) + ' '
    left = [0, 0, 0, 1, maxCoord, h]
    s += ' '.join([str(x) for x in left]) + ' '
    right = [h-1, 0, 0, h, maxCoord, h]
    s += ' '.join([str(x) for x in right]) + ' '
    back = [0, 0, 0, h, maxCoord, 1]
    s += ' '.join([str(x) for x in back]) + ' '
    front = [0, 0, h-1, h, maxCoord, h]
    s += ' '.join([str(x) for x in front]) + ' '
    print('drop3')
    print(s)
    s = 'iFluidBoxes ' + 3 * (str(res) + ' ')
    tank = [1, 1, 1, h - 1, m + 1, h - 1]
    s += ' '.join([str(x) for x in tank]) + ' '
    print(s)
    s = 'iFluidSpheres ' + 3 * (str(res) + ' ')
    drop = [h // 2, m + 2* m + h // 8, h // 2, h // 8]
    s += ' '.join([str(x) for x in drop]) + ' '
    print(s)
    s = 'iSolidBoxes ' + 2 * (str(res) + ' ')
    floor = [0, 0, h, 1]
    s += ' '.join([str(x) for x in floor]) + ' '
    left = [0, 0, 1, maxCoord]
    s += ' '.join([str(x) for x in left]) + ' '
    right = [h-1, 0, h, maxCoord]
    s += ' '.join([str(x) for x in right]) + ' '
    print('drop2')
    print(s)
    s = 'iFluidBoxes ' + 2 * (str(res) + ' ')
    tank = [1, 1, h - 1, m + 1]
    s += ' '.join([str(x) for x in tank]) + ' '
    print(s)
    s = 'iFluidSpheres ' + 2 * (str(res) + ' ')
    drop = [h // 2, m + 2 * m + h // 8, h // 8]
    s += ' '.join([str(x) for x in drop]) + ' '
    print(s)

elif scene == "teaser":
    res, h, m = int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4])
    maxCoord = int(math.pow(2, res))
    s = 'iSolidBoxes ' + 3 * (str(res) + ' ')
    n = 2
    # floor 1
    floor1 = [h-1, 0, 0, 3 * h, 1, 2 * h]
    s += ' '.join([str(x) for x in floor1]) + ' '
    floor2 = [0, 2 * h, 0, h, 2*h + 1, 2 * h]
    s += ' '.join([str(x) for x in floor2]) + ' '
    left1 = [h - 1, 0, 0, h, 2 * h, 2 * h]
    s += ' '.join([str(x) for x in left1]) + ' '
    left2 = [0, 2 * h, 0, 1, 5 * h, 2 * h]
    s += ' '.join([str(x) for x in left2]) + ' '
    back1 = [h-1, 0, 0, 3 * h, 5 * h, 1]
    s += ' '.join([str(x) for x in back1]) + ' '
    back2 = [0, 2 * h, 0, h, 5 * h, 1]
    s += ' '.join([str(x) for x in back2]) + ' '
    front1 = [h, 0, 2 * h-1, 3 * h, 5 * h, 2 * h]
    s += ' '.join([str(x) for x in front1]) + ' '
    front2 = [0, 2 * h, 2 * h-1, h, 5 * h, 2 * h]
    s += ' '.join([str(x) for x in front2]) + ' '
    right1 = [3 * h - 1, 0, 0, 3 * h, 5 * h, 2 * h]
    s += ' '.join([str(x) for x in right1]) + ' '
    right2 = [h - 1, 2 * h + n * 2, 0, h, 5 * h, 2 * h]
    s += ' '.join([str(x) for x in right2]) + ' '
    # holes
    box = [0,0,0,1,n*2+1,9]
    for i in range(h):
        offset = [h-1, 2 * h, 5 * i - 4]
        obox = []
        if i % 2:
            continue
        for j in range(6):
            obox.append(offset[j % 3] + box[j])
        s += ' '.join([str(x) for x in obox]) + ' '
    box = [0,0,0,1,4,2*h]
    for i in range(1, n):
        offset = [h-1, 2 * i + 2 * h, 0]
        obox = []
        for j in range(6):
            obox.append(offset[j % 3] + box[j])
        s += ' '.join([str(x) for x in obox]) + ' '
    print(s)
    s = 'iFluidBoxes ' + 3 * (str(res) + ' ')
    box = [1,2*h, 1, h-1, m + 2 * h, 2 * h - 1]
    s += ' '.join([str(x) for x in box]) + ' '
    box = [h, 1, 1, 3 * h-1, 4, 2 * h -1]
    s += ' '.join([str(x) for x in box]) + ' '
    print(s)

