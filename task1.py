import random
from copy import deepcopy
from typing import List

from utils import Point, Circle, draw_circle_points, getcircle_3points, getSphere_3points, getSphere_4points


def SphereMin(P: List[Point], R: List[Point]):
    '''
    SphereMin(points, []) returns the minimum bounding sphere of points
    :param P: Points in the Sphere
    :param R: Points on the surface of the Sphere
    :return: Minimum Sphere
    '''
    # base case
    if len(P) == 0 or len(R) >= 4:
        if len(R) == 1:
            p = R[0]
            mb = Circle(p, 0)
            return mb
        elif len(R) == 2:
            p0 = R[0]
            p1 = R[1]
            center = Point((p0.x + p1.x) / 2, (p0.y + p1.y) / 2, (p0.z + p1.z) / 2)
            diameter = p0.calculate_distance_to(p1)
            mb = Circle(center, diameter / 2)
            return mb
        elif len(R) == 3:
            p0 = R[0]
            p1 = R[1]
            p2 = R[2]
            mb = Circle(*getSphere_3points(p0, p1, p2))
            return mb
        elif len(R) >= 4:
            mb = Circle(*getSphere_4points(*R[:4]))
            for r in R:
                if abs(r.calculate_distance_to(mb.centre) - mb.radius) > 1e-3:
                    mb = Circle(Point(0, 0, 0), 0)
            return mb
        else:
            return Circle(Point(0, 0, 0), 0)
    # Linear Programming
    p = random.choice(P)
    P.remove(p)
    D1 = SphereMin(deepcopy(P), deepcopy(R))
    if p.calculate_distance_to(D1.centre) <= D1.radius:
        return D1
    R.append(p)
    return SphereMin(deepcopy(P), deepcopy(R))


def CircleMin(P: List[Point], R: List[Point]):
    '''
    CircleMin(points, []) returns the minimum bounding circle of points
    :param P: Points in the Cirle
    :param R: Points on the surface of the Circle
    :return: Minimum Circle
    '''
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
    D1 = CircleMin(deepcopy(P), deepcopy(R))
    if p.calculate_distance_to(D1.centre) <= D1.radius:
        return D1
    R.append(p)
    return CircleMin(deepcopy(P), deepcopy(R))


if __name__ == "__main__":
    # points = [Point(5, 0, 1), Point(-1, -3, 4), Point(-1, -4, -3), Point(-1, 4, -3)]
    # points = [Point(-10, 0,0), Point(10, 0, 0), Point(0, 1, 0)]
    # points = [Point(-5, 0,0), Point(3, -4, 0), Point(3, 4, 0)]
    points = [Point(random.randint(-100, 100), random.randint(-100, 100), random.randint(-100, 100)) for _ in
              range(100)]
    c = SphereMin(deepcopy(points), [])
    print(c)
    draw_circle_points(c, points)
