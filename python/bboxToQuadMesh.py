boxes = []
boxes.append([[0, 0, 0 ],[0.03125, 1, 0.6875 ]])
boxes.append([[0, 0, 0 ],[0.6875, 1, 0.03125 ]])
boxes.append([[0.4375, 0, 0.1875 ],[0.46875, 1, 0.6875 ]])
boxes.append([[0.21875, 0, 0 ],[0.25, 1, 0.5 ]])
boxes.append([[0.65625, 0, 0 ],[0.6875, 1, 0.6875 ]])
boxes.append([[0, 0, 0.65625 ],[0.6875, 1, 0.6875 ]])
boxes.append([[0, 0, 0 ],[0.6875, 0.03125, 0.6875 ]])
boxes.append([[0, 0.96875, 0 ],[0.6875, 1, 0.6875 ]])
print(len(boxes))
indices = []
for i in range(len(boxes)):
    box = boxes[i]
    #       g6       h7
    #   e4        f5
    #
    #       c2       d3
    #   a0        b1
    a = [box[0][0], box[0][1], box[0][2]]
    b = [box[1][0], box[0][1], box[0][2]]
    c = [box[0][0], box[0][1], box[1][2]]
    d = [box[1][0], box[0][1], box[1][2]]
    e = [box[0][0], box[1][1], box[0][2]]
    f = [box[1][0], box[1][1], box[0][2]]
    g = [box[0][0], box[1][1], box[1][2]]
    h = [box[1][0], box[1][1], box[1][2]]
    vertices = [a,b,c,d,e,f,g,h]
    for k in range(len(vertices)):
        print(str(i * 8 + k),' '.join([str(x) for x in vertices[k]]))
    indices.append([0 + i * 8, 1 + i * 8, 3 + i * 8, 2 + i * 8])
    indices.append([4 + i * 8, 5 + i * 8, 7 + i * 8, 6 + i * 8])
    indices.append([0 + i * 8, 1 + i * 8, 5 + i * 8, 4 + i * 8])
    indices.append([2 + i * 8, 3 + i * 8, 7 + i * 8, 6 + i * 8])
    indices.append([0 + i * 8, 2 + i * 8, 6 + i * 8, 4 + i * 8])
    indices.append([1 + i * 8, 3 + i * 8, 7 + i * 8, 5 + i * 8])
print(len(indices))
for index in indices:
    print(' '.join(str(i) for i in index))