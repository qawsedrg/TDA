import argparse
import sys
from collections import defaultdict
from typing import List, Dict

import gudhi
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.widgets import Slider
from mpl_toolkits.mplot3d import Axes3D
from scipy.linalg import solve


class Point():
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
            return ("({:.2f}, {:.2f},{:.2f})".format(self.x, self.y, self.z))
        else:
            return ("({:.2f}, {:.2f})".format(self.x, self.y))

    def __str__(self):
        if self.z is not None:
            return ("({:.2f}, {:.2f},{:.2f})".format(self.x, self.y, self.z))
        else:
            return ("({:.2f}, {:.2f})".format(self.x, self.y))

    def calculate_distance_to(self, p):
        x12 = (self.x - p.x) ** 2
        y12 = (self.y - p.y) ** 2
        euclidean_distance = np.sqrt(x12 + y12)
        if self.z is not None:
            z12 = (self.z - p.z) ** 2
            euclidean_distance = np.sqrt(x12 + y12 + z12)

        return euclidean_distance

    def d(self):
        if self.z is not None:
            return self.x ** 2 + self.y ** 2 + self.z ** 2
        else:
            return self.x ** 2 + self.y ** 2


class Circle():
    def __init__(self, centre: Point, radius: float):
        self.centre = centre
        self.radius = radius

    def __str__(self):
        return self.centre.__str__() + " " + str(self.radius)

    def V(self):
        return 4 / 3 * np.pi * self.radius ** 3


def getSphere_4points(p1: Point, p2: Point, p3: Point, p4: Point):
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
    for p in P:
        if p.calculate_distance_to(C.centre) < C.radius:
            return False
    return True


def get_filtration_2d(P: List[Point], C: List[int]):
    if len(C) == 1:
        return 0
    if len(C) == 2:
        return (P[C[0]].calculate_distance_to(P[C[1]]) / 2) ** 2
    if len(C) == 3:
        return Circle(*getcircle_3points(P[C[0]], P[C[1]], P[C[2]])).radius ** 2


def get_filtration_3d(P: List[Point], C: List[int]):
    if len(C) == 1:
        return 0
    if len(C) == 2:
        return (P[C[0]].calculate_distance_to(P[C[1]]) / 2) ** 2
    if len(C) == 3:
        return (Circle(*getSphere_3points(P[C[0]], P[C[1]], P[C[2]])).radius) ** 2
    if len(C) == 4:
        return (Circle(*getSphere_4points(P[C[0]], P[C[1]], P[C[2]], P[C[3]])).radius) ** 2


def show(d: Dict, points: List[Point]):
    d_reverse = defaultdict(list)
    for (k, v) in d.items():
        if v == 0:
            continue
        d_reverse[np.sqrt(v)].append(k)
    dis = sorted(d_reverse.keys())
    plt.subplots_adjust(bottom=0.3)

    fig, ax = plt.subplots()

    plt.subplots_adjust(bottom=0.3)
    ax.set_aspect(1)
    axcolor = 'lightgoldenrodyellow'
    for point in points:
        ax.scatter(point.x, point.y, c='b', s=40)
    axfreq = plt.axes([0.25, 0.1, 0.65, 0.03], facecolor=axcolor)
    sfreq = Slider(axfreq, 'Freq', 0.0, max(dis) + 1, valfmt='% .2f', valinit=0, valstep=0.001)

    def update(val):
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


def compute_diagram(action, file_complex):
    st = gudhi.SimplexTree()
    for line in file_complex:
        fields = line.split()
        st.insert([int(i) for i in fields])
    st.compute_persistence()
    for i, n in enumerate(st.betti_numbers()):
        if n != 0:
            print("dimension", i, ":", n)


def plot3d(action, file_complex, file_coordinates):
    points = np.loadtxt(file_coordinates)
    simplexes = [[], [], []]
    for line in file_complex:
        fields = line.split()
        k = len(fields) - 1
        if k <= 2:
            simplexes[k].append([int(i) for i in fields])

    from mpl_toolkits.mplot3d.art3d import Line3DCollection
    import matplotlib.pyplot as plt
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    # Plot triangles
    ax.plot_trisurf(points[:, 0], points[:, 1], points[:, 2], triangles=simplexes[2])
    # Plot points
    points2 = points[np.array(simplexes[0]).reshape(-1)]
    ax.scatter3D(points2[:, 0], points2[:, 1], points2[:, 2])
    # Plot edges
    ax.add_collection3d(Line3DCollection(segments=[points[e] for e in simplexes[1]]))
    # plt.savefig('myplot.pdf')
    plt.show()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Process a simplicial complex.')
    parser.add_argument('--complex', default="./cplx.txt", type=argparse.FileType('r'))
    parser.add_argument('--coordinates', default="./coord.txt", type=argparse.FileType('r'))
    parser.add_argument('--action', default="plot3d", choices=['betti', 'plot3d'])
    args = parser.parse_args()
    if args.action in ['betti']:
        if not args.complex:
            print("missing complex")
            sys.exit()
        compute_diagram(args.action, args.complex)
    elif args.action in ['plot3d']:
        if not args.complex:
            print("missing complex")
            sys.exit()
        if not args.coordinates:
            print("missing coordinates")
            sys.exit()
        plot3d(args.action, args.complex, args.coordinates)
