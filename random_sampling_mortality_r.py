from __future__ import print_function
import numpy as np
import json
from multiprocessing import Pool
import datetime
import sys
from argparse import ArgumentParser

from matplotlib import rc

import math
import random
import subprocess
import sys
import os
import bisect

parser = ArgumentParser()
parser.add_argument('num_samples', metavar='num_samples', type=int)
parser.add_argument('-np', dest='num_processes', default=1, type=int,
                    help='Num of processors to use')
parsed = parser.parse_args()


def model_cluster_trip(config_file_name, seed, devnull):
    if sys.version_info > (3, 0):
        return subprocess.run(["./model_cluster_trip_v2", config_file_name, str(seed)],
                stdout=subprocess.PIPE, stderr=devnull).stdout
    else:
        p = subprocess.Popen(["./model_cluster_trip_v2", config_file_name, str(seed)],
                stdout=subprocess.PIPE, stderr=devnull)
        return "".join(p.stdout.readlines())


def f(sample_id):
    devnull = open(os.devnull, 'w')
    with open("config_for_fat_tail.json") as f:
        config = json.load(f)

    C = 0.002209375
    while True:
        prob_transmission = random.uniform(0, 1)
        prob_i_to_ic = random.uniform(0, 1)
        if prob_transmission * prob_i_to_ic < C: break

    #C = 0.1 * (0.72 * 0.00471 * 0.286 + 0.28 * 0.11285 * 0.286)

    ##x = math.exp(random.uniform(C * math.log(C), C * math.log(1)) / C)
    ##prob_transmission = x
    ##prob_i_to_ic = C / x

    #while True:
    #    prob_transmission = random.uniform(0, 1)
    #    prob_i_to_ic = random.uniform(0, 1)
    #    if prob_transmission * prob_i_to_ic < C: break

    #while True:
    #    C = np.random.normal(0.0001, 7.03125e-4)
    #    if C > 0: break
    #points_on_curve = []
    #for x in np.arange(C, 1, (1 - C) / 100000):
    #    points_on_curve.append((x, C / x))
    #    points_on_curve.append((C / x, x))
    #points_on_curve = sorted(points_on_curve)
    #arc_length_sum = [0.0]
    #for i in range(1, len(points_on_curve)):
    #    x1, y1 = points_on_curve[i - 1]
    #    x2, y2 = points_on_curve[i]
    #    d = math.sqrt((x1 - x2) ** 2 + (y1 - y2) ** 2)
    #    arc_length_sum.append(d + arc_length_sum[i - 1])
    #l = random.uniform(0, arc_length_sum[-1])
    #idx = bisect.bisect_left(arc_length_sum, l)
    #prob_transmission, prob_i_to_ic = points_on_curve[idx]

    config["simulation"]["prob_transmission"] = prob_transmission
    config["simulation"]["initial_params"][0]["prob_i_to_ic"] = prob_i_to_ic

    config_file_name = "tmp/random_sampling_{}.json".format(sample_id)
    with open(config_file_name, "w") as f:
        json.dump(config, f, indent=4)

    print("Running {}/{} with params: C = {:.5f}, prob_transmission = {:.5f}, prob_i_to_ic = {:.5f}"
            .format(sample_id, parsed.num_samples, C, prob_transmission,
                prob_i_to_ic), file=sys.stderr)

    stdout = model_cluster_trip(config_file_name, 0, devnull)
    output=json.loads(stdout)
    os.remove(config_file_name)

    with open("outputs/random_sampling/{}.json".format(sample_id), "w") as f:
        json.dump(output, f, indent=4)

    devnull.close()


if parsed.num_processes == 1:
    print("Running ONE process", file=sys.stderr)
    for sample_id in range(parsed.num_samples):
        f(sample_id)
else:
    print("Running {} processes".format(parsed.num_processes), file=sys.stderr)
    # This won't work on isabella, because python2 has different API for Pool
    # object. However, that is fine since we run parallelize our code with SGE
    # scripts.
    with Pool(parsed.num_processes) as p:
        p.map(f, range(parsed.num_samples))
