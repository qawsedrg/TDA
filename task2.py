from copy import deepcopy
from itertools import combinations
from typing import List

from task1 import task1
from utils import Point


def task2(P: List[Point]):
    result = dict()
    for i in range(len(P)):
        for tup in list(combinations(range(len(P)), i + 1)):
            if len(tup) == 1:
                print("({:})->[{:}]".format(tup[0], 0))
                result[tup] = 0
            else:
                c = task1(deepcopy([P[i] for i in tup]), [])
                print(tup, end="")
                print("->[{:.5f}]".format(c.radius))
                result[tup] = c.radius ** 2
    return result


if __name__ == "__main__":
    points = [Point(5, 0, 1), Point(-1, -3, 4), Point(-1, -4, -3), Point(-1, 4, -3)]
    # points = [Point(-10, 0,0),Point(10, 0, 0),Point(0, 1, 0)]
    # points = [Point(-5, 0,0),Point(3, -4, 0),Point(3, 4, 0)]
    task2(points)
