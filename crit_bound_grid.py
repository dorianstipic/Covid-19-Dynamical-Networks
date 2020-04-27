# Autor: Dorian
from __future__ import print_function
import numpy as np
import json
from multiprocessing import Pool
import datetime
import sys
from argparse import ArgumentParser

from matplotlib import rc

import subprocess
import sys
import os

parser = ArgumentParser()
parser.add_argument('cluster_sizes', metavar='cluster_size', type=int, nargs='+')
parser.add_argument('-np', dest='num_processes', default=1, type=int,
                    help='Num of processors to use')
parser.add_argument('-mu', dest='mu', default=100, type=int,
                    help='mu parametar in model config')
parsed = parser.parse_args()

def model_cluster_trip(config_file_name, seed, devnull):
    if sys.version_info > (3, 0):
        return subprocess.run(["./model_cluster_trip_v2", config_file_name, str(seed)],
                stdout=subprocess.PIPE, stderr=devnull).stdout
    else:
        p = subprocess.Popen(["./model_cluster_trip_v2", config_file_name, str(seed)],
                stdout=subprocess.PIPE, stderr=devnull)
        return "".join(p.stdout.readlines())


def f_crit(cluster_size):
    with open("config_for_grid_search.json") as f:
        config = json.load(f)

    mu = parsed.mu

    scale = 1 #1 if real simul
    num_icus=int(200/scale)#200 baseline
    num_nodes=int(1000000/scale)#1,000,000 baseline
    num_clusters=int(num_nodes/cluster_size)

    scaledays = 1 # 1 if real simul
    num_days=int(1200/scaledays) #1200 baseline, we want most, even controlled pandemic to end.

    config["graph_generation"]["num_people_per_cluster"]=cluster_size
    config["graph_generation"]["num_clusters"]=num_clusters
    config["simulation"]["stopping_conditions"]["num_days"]=num_days
    config["simulation"]["num_icus"]=num_icus
    config["simulation"]["mu"] = mu
    config["simulation"]["events"][0]["prob_s_to_i"]=[el * scale for el in config["simulation"]["events"][0]["prob_s_to_i"]]

    config["simulation"]["stopping_conditions"]["on_icu_overflow"]=True # Important!

    h=0.5 #0.5 for precise grid.
    ptrip=np.arange(0,1.00001,0.01*h)
    #pdisobedient=np.arange(0,1.00001,0.01*h)
    #seeds=np.arange(0,1,1) #step 1 - 100 possible
    step=0.01*h
    seed=0

    bounds_dict = {}
    p1_list = []
    p2_list = []
    p2 = 1 # we start from the top.

    start = datetime.datetime.now()
    devnull = open(os.devnull, 'w')
    for p1 in ptrip:
        while p2>-0.1*step: #suppose float approx is 10% or less of step
            config["simulation"]["prob_goes_on_trip"] = p1
            config["simulation"]["initial_params"][0]["prob_c_neighbour_trip_candidate"] = p2
            config["simulation"]["initial_params"][1]["prob_c_neighbour_trip_candidate"] = p2

            config_file_name = "tmp/tmp_config_cbg_{}_{}_{}.json".format(
                    mu, cluster_size, seed)
            with open(config_file_name, "w") as f:
                json.dump(config, f, indent=4)

            print("Running model with params: cluster_size = {:.3f}".format(cluster_size),
                  ", prob_goes_on_trip = {:.3f}".format(p1),
                  ", prob_c_neighbour_trip_candidate = {:.3f}".format(p2),
                  "seed = {}".format(seed), file=sys.stderr)
            stdout = model_cluster_trip(config_file_name, seed, devnull)

            output=json.loads(stdout)
            os.remove(config_file_name)

            if output["stopping_condition"]=="icu_overflow":
                p2=p2-step # we go one step down.
                if p2>-0.1*step and p2<0: # if -0.1*step<p2<0 means we need to put p2=0, next time after p2=p2-step it will exit while
                    p2=0
                continue
            else:
                p1_list.append(p1)
                p2_list.append(p2)
                break
    bounds_dict["p1_vrijednosti"]=p1_list
    bounds_dict["p2_vrijednosti"]=p2_list

    with open("outputs/crit_bound_search_mu{}_cluster_size{}.json".format(
            parsed.mu, cluster_size), "w") as f:
        json.dump(bounds_dict, f, indent=4)

    devnull.close()

    end = datetime.datetime.now()
    print("Time elapsed during the calculation:", end - start)


if parsed.num_processes == 1:
    print("Running ONE process", file=sys.stderr)
    for cluster_size in parsed.cluster_sizes:
        f_crit(cluster_size)
else:
    print("Running {} processes".format(parsed.num_processes), file=sys.stderr)
    # This won't work on isabella, because python2 has different API for Pool
    # object. However, that is fine since we run parallelize our code with SGE
    # scripts.
    with Pool(parsed.num_processes) as p:
        p.map(f_crit, parsed.cluster_sizes)
