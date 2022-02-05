from copy import deepcopy
from itertools import combinations
from typing import List

import numpy as np
from scipy.spatial import Delaunay

from utils import Point, Circle, get_filtration_2d, getSphere_3points, get_filtration_3d, is_Gabriel


def Alpha_2d(P: List[Point]):
    """
        return simplexes of the corresponding alpha complexe
    """
    fil = dict()
    tri = Delaunay(np.array([[p.x, p.y] for p in P]))
    S = [deepcopy(set()) for i in range(len(tri.simplices[0]))]
    S[-1] = set([tuple(sorted(simplice)) for simplice in tri.simplices])
    for i in range(len(tri.simplices[0]), 0, -1):
        while len(S[i - 1]) > 0:
            C = S[i - 1].pop()
            if not C in fil.keys():
                fil[C] = get_filtration_2d(P, C)
            for tup in list(combinations(C, len(C) - 1)):
                if len(tup) == 1:
                    fil[tup] = 0
                    continue
                S[i - 2].add(tup)
                if tup in fil.keys():
                    fil[tup] = min(fil[tup], fil[C])
                else:
                    if len(tup) == 2:
                        p0 = P[tup[0]]
                        p1 = P[tup[1]]
                        if not is_Gabriel([P[j] for j in C], Circle(Point((p0.x + p1.x) / 2, (p0.y + p1.y) / 2),
                                                                    p0.calculate_distance_to(p1) / 2)):
                            fil[tup] = fil[C]
    return fil


def Alpha_3d(P: List[Point]):
    """
    return simplexes of the corresponding alpha complexe
    """
    fil = dict()
    tri = Delaunay(np.array([[p.x, p.y, p.z] for p in points]))
    S = [deepcopy(set()) for i in range(len(tri.simplices[0]))]
    S[-1] = set([tuple(sorted(simplice)) for simplice in tri.simplices])
    for i in range(len(tri.simplices[0]), 0, -1):
        while len(S[i - 1]) > 0:
            C = S[i - 1].pop()
            if not C in fil.keys():
                fil[C] = get_filtration_3d(P, C)
            for tup in list(combinations(C, len(C) - 1)):
                if len(tup) == 1:
                    fil[tup] = 0
                    continue
                S[i - 2].add(tup)
                if tup in fil.keys():
                    fil[tup] = min(fil[tup], fil[C])
                else:
                    if len(tup) == 2:
                        p0 = P[tup[0]]
                        p1 = P[tup[1]]
                        if not is_Gabriel([P[j] for j in C],
                                          Circle(Point((p0.x + p1.x) / 2, (p0.y + p1.y) / 2, (p0.z + p1.z) / 2),
                                                 p0.calculate_distance_to(p1) / 2)):
                            fil[tup] = fil[C]
                    if len(tup) == 3:
                        p0 = P[tup[0]]
                        p1 = P[tup[1]]
                        p2 = P[tup[2]]
                        if not is_Gabriel([P[j] for j in C],
                                          Circle(*getSphere_3points(p0, p1, p2))):
                            fil[tup] = fil[C]
    return fil


def Alpha(P: List[Point], k: int, l: float, dim: int):
    """
    :param P: points
    :param k: max dim of simplexe
    :param l: max filtration value
    :param dim: dimension of the space of points (2 or 3)
    :return: dict of simplexes and corresponding filtration value
    """
    fil = Alpha_2d(P) if dim == 2 else Alpha_3d(P)
    d = dict()
    for (key, value) in fil.items():
        if len(key) <= k + 1 and value <= l:
            d[key] = value
    return d


if __name__ == "__main__":
    points = [Point(1, 1), Point(7, 0), Point(4, 6), Point(9, 6), Point(0, 14), Point(2, 19), Point(9, 17)]
    d = Alpha(points, 2, 10, 2)
    for (k, v) in d.items():
        if v == 0:
            print("{:}->[out]".format(k))
        else:
            print("{:}->[{:.5f}]".format(k, v))
