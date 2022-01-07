from copy import deepcopy
from itertools import combinations
from typing import List

import gudhi
import numpy as np
from scipy.spatial import Delaunay

from utils import Point, Circle, is_Gabriel, get_filtration_2d


def task3_2d(P: List[Point]):
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


if __name__ == "__main__":
    points = [Point(1, 1), Point(7, 0), Point(4, 6), Point(9, 6), Point(0, 14), Point(2, 19), Point(9, 17)]
    d = task3_2d(points)
    alpha_complex = gudhi.AlphaComplex(points=[[p.x, p.y] for p in points])
    simplex_tree = alpha_complex.create_simplex_tree()
    for filtered_value in simplex_tree.get_filtration():
        print('%s -> %.2f' % tuple(filtered_value), end=" ")
        print("%.2f" % (
            d.get(tuple(sorted([points.index(Point(*alpha_complex.get_point(i))) for i in filtered_value[0]])), 0)))
