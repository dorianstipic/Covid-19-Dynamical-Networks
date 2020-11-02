import numpy as np
import matplotlib.pyplot as plt
import math
import random
import bisect

def get_point(C):
    points_on_curve = []
    for x in np.arange(C, 1, (1 - C) / 10000):
        points_on_curve.append((x, C / x))
        points_on_curve.append((C / x, x))
    points_on_curve = sorted(points_on_curve)
    arc_length_sum = [0.0]
    for i in range(1, len(points_on_curve)):
        x1, y1 = points_on_curve[i - 1]
        x2, y2 = points_on_curve[i]
        d = math.sqrt((x1 - x2) ** 2 + (y1 - y2) ** 2)
        arc_length_sum.append(d + arc_length_sum[i - 1])
    l = random.uniform(0, arc_length_sum[-1])
    idx = bisect.bisect_left(arc_length_sum, l)
    prob_transmission, prob_i_to_ic = points_on_curve[idx]
    return prob_transmission, prob_i_to_ic

xs = []
ys = []
for i in range(100):
    x, y = get_point(0.001)
    xs.append(x)
    ys.append(y)

d = []
for i in range(1, len(xs)):
    x1, y1 = xs[i - 1], ys[i - 1]
    x2, y2 = xs[i], ys[i]
    d.append(math.sqrt((x1 - x2) ** 2 + (y1 - y2) ** 2))
plt.scatter(xs[1:], d)
plt.show()
