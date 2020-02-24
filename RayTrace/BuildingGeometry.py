import numpy as np

FaceNormalNo = 5
FaceNormals = [(-1, 0, 0), (0, 1, 0), (1, 0, 0), (0, -1, 0), (0, 0, 1)]
# ^Will's Code

BoxNumber = 1
BoxArrayNear = np.zeros([BoxNumber, 3])
BoxArrayFar = np.zeros([BoxNumber, 3])
BoxArrayNear[0] = [10, 10, 0]
BoxArrayFar[0] = [64.4322, 46.9316, 8.2423]
TriangleNumber = 0
SquareNumber = 0
PolyBuilding = 0
