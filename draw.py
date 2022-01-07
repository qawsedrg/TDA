import numpy as np

from task2 import task2
from utils import Point, show
from task3_2d import task3_2d

points = [Point(1, 1), Point(7, 0), Point(4, 6), Point(9, 6), Point(0, 14), Point(2, 19), Point(9, 17)]
show(task3_2d(points), points)
points = [Point(1, 0), Point(0, 1), Point(2, 1), Point(3, 2), Point(0, 3), Point(3 + np.sqrt(3), 3), Point(1, 4),
          Point(3, 4), Point(2, 4 + np.sqrt(3)), Point(0, 4), Point(-0.5, 2)]
show(task2(points), points)
