from collections import defaultdict
from typing import List, Dict, Tuple

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.widgets import Slider
from mpl_toolkits.mplot3d import Axes3D
from scipy.linalg import solve


class Point:
    """
    A representation of a point in 2D or 3d
    """

    def __init__(self, *coordinates):
        self.x = coordinates[0]
        self.y = coordinates[1]
        self.z = None
        if len(coordinates) == 3:
            self.z = coordinates[2]

    def __eq__(self, other):
        if abs(self.x - other.x) <= 1e-5 and abs(self.y - other.y) <= 1e-5 and (
                (self.z is not None and abs(self.z - other.z) <= 1e-5) or self.z is None):
            return True
        return False

    def __repr__(self):
        if self.z is not None:
            return "({:.2f}, {:.2f},{:.2f})".format(self.x, self.y, self.z)
        else:
            return "({:.2f}, {:.2f})".format(self.x, self.y)

    def __str__(self):
        if self.z is not None:
            return "({:.2f}, {:.2f},{:.2f})".format(self.x, self.y, self.z)
        else:
            return "({:.2f}, {:.2f})".format(self.x, self.y)

    def calculate_distance_to(self, p):
        x12 = (self.x - p.x) ** 2
        y12 = (self.y - p.y) ** 2
        euclidean_distance = np.sqrt(x12 + y12)
        if self.z is not None and p.z is not None:
            z12 = (self.z - p.z) ** 2
            euclidean_distance = np.sqrt(x12 + y12 + z12)

        return euclidean_distance

    def d(self):
        if self.z is not None:
            return self.x ** 2 + self.y ** 2 + self.z ** 2
        else:
            return self.x ** 2 + self.y ** 2


class Circle:
    """
    A representation of a sphere in all dimensions
    """

    def __init__(self, centre: Point, radius: float):
        self.centre = centre
        self.radius = radius

    def __str__(self):
        return self.centre.__str__() + " " + str(self.radius)


def getSphere_4points(p1: Point, p2: Point, p3: Point, p4: Point):
    """
    return the sphere defined by these points.
    Donot use on points which are on the same plane
    :param p1: point
    :param p2: point
    :param p3: point
    :param p4: point
    :return: Sphere defined by these points
    """
    a1 = 2 * (p1.x - p2.x)
    b1 = 2 * (p1.y - p2.y)
    c1 = 2 * (p1.z - p2.z)
    d1 = p2.d() - p1.d()

    a2 = 2 * (p2.x - p3.x)
    b2 = 2 * (p2.y - p3.y)
    c2 = 2 * (p2.z - p3.z)
    d2 = p3.d() - p2.d()

    a3 = 2 * (p3.x - p4.x)
    b3 = 2 * (p3.y - p4.y)
    c3 = 2 * (p3.z - p4.z)
    d3 = p4.d() - p3.d()

    coeff1 = [[a1, b1, c1], [a2, b2, c2], [a3, b3, c3]]
    coeff2 = [-d1, -d2, -d3]

    p = Point(*solve(coeff1, coeff2))
    r = p1.calculate_distance_to(p)

    return p, r


def getSphere_3points(p1: Point, p2: Point, p3: Point):
    """
    return minimum sphere defined by these points.
    Donot use on points which are on the same line
    :param p1: point
    :param p2: point
    :param p3: point
    :return: minimum sphere defined by these points
    """
    a3 = (p2.y - p1.y) * (p3.z - p1.z) - (p2.z - p1.z) * (p3.y - p1.y)
    b3 = (p2.z - p1.z) * (p3.x - p1.x) - (p2.x - p1.x) * (p3.z - p1.z)
    c3 = (p2.x - p1.x) * (p3.y - p1.y) - (p2.y - p1.y) * (p3.x - p1.x)
    d3 = -a3 * p1.x - b3 * p1.y - c3 * p1.z

    a1 = 2 * (p1.x - p2.x)
    b1 = 2 * (p1.y - p2.y)
    c1 = 2 * (p1.z - p2.z)
    d1 = p2.d() - p1.d()

    a2 = 2 * (p2.x - p3.x)
    b2 = 2 * (p2.y - p3.y)
    c2 = 2 * (p2.z - p3.z)
    d2 = p3.d() - p2.d()

    coeff1 = [[a1, b1, c1], [a2, b2, c2], [a3, b3, c3]]
    coeff2 = [-d1, -d2, -d3]

    p = Point(*solve(coeff1, coeff2))
    r = p1.calculate_distance_to(p)

    return p, r


def getcircle_3points(p1: Point, p2: Point, p3: Point):
    """
    return the circle defined by these points.
    Donot use on points which are on the same line
    :param p1: point
    :param p2: point
    :param p3: point
    :return: circle defined by these points
    """
    a1 = 2 * (p1.x - p2.x)
    b1 = 2 * (p1.y - p2.y)
    d1 = p2.d() - p1.d()

    a2 = 2 * (p2.x - p3.x)
    b2 = 2 * (p2.y - p3.y)
    d2 = p3.d() - p2.d()

    coeff1 = [[a1, b1], [a2, b2]]
    coeff2 = [-d1, -d2]

    p = Point(*solve(coeff1, coeff2))
    r = p1.calculate_distance_to(p)

    return p, r


def draw_circle_points(c: Circle, points: List[Point]):
    """
    draw in a 3D space a circle and some points
    """

    def drawSphere(c: Circle):
        xCenter = c.centre.x
        yCenter = c.centre.y
        zCenter = c.centre.z
        r = c.radius
        u, v = np.mgrid[0:2 * np.pi:20j, 0:np.pi:10j]
        x = np.cos(u) * np.sin(v)
        y = np.sin(u) * np.sin(v)
        z = np.cos(v)
        x = r * x + xCenter
        y = r * y + yCenter
        z = r * z + zCenter
        return (x, y, z)

    fig = plt.figure()
    ax = Axes3D(fig)

    ax.set_xlim(left=c.centre.x - 2 * c.radius, right=c.centre.x + 2 * c.radius)
    ax.set_ylim(bottom=c.centre.y - 2 * c.radius, top=c.centre.y + 2 * c.radius)
    ax.set_zlim(bottom=c.centre.z - 2 * c.radius, top=c.centre.z + 2 * c.radius)

    ax.scatter(c.centre.x, c.centre.y, c.centre.z, c='r')
    for point in points:
        ax.scatter(point.x, point.y, point.z, c='b', s=40)

    ax.plot_wireframe(*drawSphere(c), color="r")
    plt.show()


def is_Gabriel(P: List[Point], C: Circle):
    """
    return True id C is Gabriel
    """
    for p in P:
        if p.calculate_distance_to(C.centre) < C.radius:
            return False
    return True


def get_filtration_2d(P: List[Point], C: List[int]):
    """
    :param P: points
    :param C: indexes of points
    :return: squared radius of the minimum bounding circle
    """
    if len(C) == 1:
        return 0
    if len(C) == 2:
        return (P[C[0]].calculate_distance_to(P[C[1]]) / 2) ** 2
    if len(C) == 3:
        return Circle(*getcircle_3points(P[C[0]], P[C[1]], P[C[2]])).radius ** 2


def get_filtration_3d(P: List[Point], C: List[int]):
    """
    :param P: points
    :param C: indexes of points
    :return: squared radius of the minimum bounding sphere
    """
    if len(C) == 1:
        return 0
    if len(C) == 2:
        return (P[C[0]].calculate_distance_to(P[C[1]]) / 2) ** 2
    if len(C) == 3:
        return Circle(*getSphere_3points(P[C[0]], P[C[1]], P[C[2]])).radius ** 2
    if len(C) == 4:
        return Circle(*getSphere_4points(P[C[0]], P[C[1]], P[C[2]], P[C[3]])).radius ** 2


def contains(tup: Tuple, outs: List[Tuple]):
    """
    :param tup: a simplexe
    :param outs: list of simplexes
    :return: if there is a simplexe in outs which is sous-simplexe of tup
    """
    for out in outs:
        if set(out).issubset(set(tup)):
            return True
    return False


def show(d: Dict, points: List[Point]):
    """
    show the complexe filtered
    :param d: simplexes and corresponding filtration value
    :param points: points
    """
    d_reverse = defaultdict(list)
    for (k, v) in d.items():
        if v == 0:
            continue
        d_reverse[np.sqrt(v)].append(k)
    dis = sorted(d_reverse.keys())

    fig, ax = plt.subplots()

    plt.subplots_adjust(bottom=0.3)
    ax.set_aspect(1)
    axcolor = 'lightgoldenrodyellow'
    for point in points:
        ax.scatter(point.x, point.y, c='b', s=40)
    axfreq = plt.axes([0.25, 0.1, 0.65, 0.03], facecolor=axcolor)
    sfreq = Slider(axfreq, 'Filtration', 0.0, max(dis) + 1, valfmt='% .2f', valinit=0, valstep=0.001)

    def update(val):
        # recalculate the complexe filtered after changing the filtration value (sfreq)
        freq = sfreq.val
        ax.clear()
        for point in points:
            ax.scatter(point.x, point.y, c='b', s=40)
            ax.add_artist(plt.Circle((point.x, point.y), radius=freq, alpha=.1))
        for di in dis:
            if di <= freq:
                simplexes = d_reverse[di]
                for simplex in simplexes:
                    if len(simplex) == 2:
                        p1 = points[simplex[0]]
                        p2 = points[simplex[1]]
                        ax.plot([p1.x, p2.x], [p1.y, p2.y], color='r')
                    if len(simplex) == 3:
                        p1 = points[simplex[0]]
                        p2 = points[simplex[1]]
                        p3 = points[simplex[2]]
                        ax.tripcolor([p1.x, p2.x, p3.x], [p1.y, p2.y, p3.y], [[0, 1, 2]],
                                     facecolors=np.array([100000]))
        fig.canvas.draw_idle()

    sfreq.on_changed(update)
    sfreq.reset()
    sfreq.set_val(0)

    plt.show()
