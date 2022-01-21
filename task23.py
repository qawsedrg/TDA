from copy import deepcopy
from itertools import combinations
from typing import List

from task1 import SphereMin, CircleMin
from utils import Point, contains


def Cech(P: List[Point], k: int, l: float, dim=2):
    result = dict()
    for i in range(min(k + 1, len(P))):
        for tup in list(combinations(range(len(P)), i + 1)):
            if len(tup) == 1:
                print("({:})->[{:}]".format(tup[0], 0))
                result[tup] = 0
            else:
                if dim == 2:
                    c = CircleMin(deepcopy([P[i] for i in tup]), [])
                if dim == 3:
                    c = SphereMin(deepcopy([P[i] for i in tup]), [])
                print(tup, end="")
                if c.radius > l:
                    print("->[out]")
                else:
                    print("->[{:.5f}]".format(c.radius))
                    result[tup] = c.radius ** 2
    return result


def CechOptimized(P: List[Point], k: int, l: float, dim=2):
    out = []
    result = dict()
    for i in range(min(k + 1, len(P))):
        for tup in list(combinations(range(len(P)), i + 1)):
            if len(tup) == 1:
                print("({:})->[{:}]".format(tup[0], 0))
                result[tup] = 0
            else:
                if dim == 2:
                    c = CircleMin(deepcopy([P[i] for i in tup]), [])
                if dim == 3:
                    c = SphereMin(deepcopy([P[i] for i in tup]), [])
                print(tup, end="")
                if contains(tup, out):
                    print("->[out]")
                elif c.radius > l:
                    out.append(tup)
                    print("->[out]")
                else:
                    print("->[{:.5f}]".format(c.radius))
                    result[tup] = c.radius ** 2
    return result


if __name__ == "__main__":
    points = [Point(5, 0, 1), Point(-1, -3, 4), Point(-1, -4, -3), Point(-1, 4, -3)]
    # points = [Point(-10, 0,0),Point(10, 0, 0),Point(0, 1, 0)]
    # points = [Point(-5, 0,0),Point(3, -4, 0),Point(3, 4, 0)]
    CechOptimized(points, k=3, l=4, dim=3)
