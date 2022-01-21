from Cech import CechOptimized
from Alpha import Alpha
from utils import Point, show

# points = [Point(1, 0), Point(0, 1), Point(2, 1), Point(3, 2), Point(0, 3), Point(3 + np.sqrt(3), 3), Point(1, 4),
#          Point(3, 4), Point(2, 4 + np.sqrt(3)), Point(0, 4), Point(-0.5, 2)]
points = [Point(1, 1), Point(7, 0), Point(4, 6), Point(9, 6), Point(0, 14), Point(2, 19), Point(9, 17)]
show(Alpha(points, k=10, l=10, dim=2), points)
show(CechOptimized(points, k=10, l=10, dim=2), points)
