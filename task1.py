import random
from copy import deepcopy
from typing import List

from utils import Point, Circle, draw_circle_points, getcircle_3points

'''
def task1(P: List[Point], R: List[Point]):
    print(len(P),len(R))
    if len(P) == 0 or len(R) >= 4:
        if len(R) == 1:
            p = R[0]
            return Circle(p, 0)
        elif len(R) == 2:
            p0 = R[0]
            p1 = R[1]
            center = Point((p0.x + p1.x) / 2, (p0.y + p1.y) / 2, (p0.z + p1.z) / 2)
            diameter = p0.calculate_distance_to(p1)
            return Circle(center, diameter / 2)
        elif len(R) == 3:
            p0 = R[0]
            p1 = R[1]
            p2 = R[2]
            return Circle(*getSphere_3points(p0, p1, p2))
        elif len(R) >= 4:
            c = Circle(*getSphere_4points(*R[:4]))
            for r in R:
                if abs(r.calculate_distance_to(c.centre) - c.radius) > 1e-3:
                    return None
            return c
        else:
            return None
    p = random.choice(P)
    P.remove(p)
    D1 = task1(deepcopy(P), deepcopy(R))
    R.append(p)
    D2 = task1(deepcopy(P), deepcopy(R))
    if D1 is not None and p.calculate_distance_to(D1.centre) > D1.radius:
        D1 = None
    if D1 is not None and D2 is not None:
        return D1 if D1.radius < D2.radius else D2
    if D1 is None and D2 is None:
        return None
    if D1 is None and D2 is not None:
        return D2
    else:
        return D1
'''
'''
def task1(L: List[Point], R: List[Point],n:int):
    mb=Circle(Point(0,0,0),0)
    if len(R) == 1:
        p = R[0]
        mb=Circle(p, 0)
    elif len(R) == 2:
        p0 = R[0]
        p1 = R[1]
        center = Point((p0.x + p1.x) / 2, (p0.y + p1.y) / 2, (p0.z + p1.z) / 2)
        diameter = p0.calculate_distance_to(p1)
        mb=Circle(center, diameter / 2)
    elif len(R) == 3:
        p0 = R[0]
        p1 = R[1]
        p2 = R[2]
        mb= Circle(*getSphere_3points(p0, p1, p2))
    elif len(R) >= 4:
        mb = Circle(*getSphere_4points(*R[:4]))
        for r in R:
            if abs(r.calculate_distance_to(mb.centre) - mb.radius) > 1e-3:
                mb=Circle(Point(0,0,0),0)
    if len(R)==4:
        return mb
    for i in range(1,n+1):
        p=L[i-1]
        if p.calculate_distance_to(mb.centre) > mb.radius:
            R_tmp=deepcopy(R)
            R_tmp.append(p)
            mb=task1(deepcopy(L),R_tmp,i-1)
            L.remove(p)
            L.insert(0,p)
    return mb
'''
'''
def task1(P: List[Point], R: List[Point]):
    if len(P) == 0 or len(R) >= 4:
        if len(R) == 1:
            p = R[0]
            mb=Circle(p, 0)
            return mb
        elif len(R) == 2:
            p0 = R[0]
            p1 = R[1]
            center = Point((p0.x + p1.x) / 2, (p0.y + p1.y) / 2, (p0.z + p1.z) / 2)
            diameter = p0.calculate_distance_to(p1)
            mb=Circle(center, diameter / 2)
            return mb
        elif len(R) == 3:
            p0 = R[0]
            p1 = R[1]
            p2 = R[2]
            mb= Circle(*getSphere_3points(p0, p1, p2))
            return mb
        elif len(R) >= 4:
            mb = Circle(*getSphere_4points(*R[:4]))
            for r in R:
                if abs(r.calculate_distance_to(mb.centre) - mb.radius) > 1e-3:
                    mb=Circle(Point(0,0,0),0)
            return mb
        else:
            return Circle(Point(0,0,0),0)
    p = random.choice(P)
    P.remove(p)
    D1 = task1(deepcopy(P), deepcopy(R))
    if p.calculate_distance_to(D1.centre) <= D1.radius:
        return D1
    R.append(p)
    return task1(deepcopy(P), deepcopy(R))
'''


def task1(P: List[Point], R: List[Point]):
    if len(P) == 0 or len(R) >= 4:
        if len(R) == 1:
            p = R[0]
            mb = Circle(p, 0)
            return mb
        elif len(R) == 2:
            p0 = R[0]
            p1 = R[1]
            center = Point((p0.x + p1.x) / 2, (p0.y + p1.y) / 2)
            diameter = p0.calculate_distance_to(p1)
            mb = Circle(center, diameter / 2)
            return mb
        elif len(R) == 3:
            p0 = R[0]
            p1 = R[1]
            p2 = R[2]
            mb = Circle(*getcircle_3points(p0, p1, p2))
            return mb
        else:
            return Circle(Point(0, 0), 0)
    p = random.choice(P)
    P.remove(p)
    D1 = task1(deepcopy(P), deepcopy(R))
    if p.calculate_distance_to(D1.centre) <= D1.radius:
        return D1
    R.append(p)
    return task1(deepcopy(P), deepcopy(R))


if __name__ == "__main__":
    points = [Point(5, 0, 1), Point(-1, -3, 4), Point(-1, -4, -3), Point(-1, 4, -3)]
    # points = [Point(-10, 0,0),Point(10, 0, 0),Point(0, 1, 0)]
    # points = [Point(-5, 0,0),Point(3, -4, 0),Point(3, 4, 0)]
    c = task1(deepcopy(points), [])
    print(c)
    draw_circle_points(c, points)
