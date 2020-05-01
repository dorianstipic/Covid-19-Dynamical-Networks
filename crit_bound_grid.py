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
parser.add_argument('-k', dest='k', default=10, type=float, help='k_trip parametar in model config')
parser.add_argument('-ext', dest='ext', default=0, type=int, help='extension type in simulation, 0:base model, 1:superspreaders, 2:domovi')
parser.add_argument('-extpop', dest='extpop', default=0, type=int, help='population in extra subgraph')
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

    ext=parsed.ext
    mu=parsed.mu
    k=parsed.k
    extpop=parsed.extpop
    
    scale = 1 #1 if real simul
    scaledays = 1 # 1 if real simul
    
    if ext==0: # base model
        with open("config_for_grid_search.json") as f:
            config = json.load(f)

        num_icus=int(200/scale)#200 baseline
        num_nodes=int(1000000/scale)#1,000,000 baseline
        num_clusters=int(num_nodes/cluster_size)

        num_days=int(1200/scaledays) #1200 baseline, we want most, even controlled pandemic to end.

        config["graph_generation"][0]["num_people_per_cluster"]=cluster_size
        config["graph_generation"][0]["num_clusters"]=num_clusters
        config["simulation"]["stopping_conditions"]["num_days"]=num_days
        config["simulation"]["num_icus"]=num_icus
        config["simulation"]["mu"] = mu
        config["simulation"]["k_trip"] = k
        config["simulation"]["events"][0]["update_params"]["prob_s_to_i"]=[el * scale for el in config["simulation"]["events"][0]["update_params"]["prob_s_to_i"]]

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
                config["simulation"]["initial_params"][0]["prob_goes_on_trip"] = p1
                config["simulation"]["initial_params"][1]["prob_goes_on_trip"] = p1
                config["simulation"]["initial_params"][0]["prob_c_neighbour_trip_candidate"] = p2
                config["simulation"]["initial_params"][1]["prob_c_neighbour_trip_candidate"] = p2

                config_file_name = "tmp/BASE_tmp_config_cbg_{}_{}_{}_{}.json".format(
                        k, mu, cluster_size, seed)
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

        with open("outputs/base_model/crit_bound_search_k_trip{}_mu{}_cluster_size{}.json".format(
                parsed.k, parsed.mu, cluster_size), "w") as f:
            json.dump(bounds_dict, f, indent=4)

        devnull.close()

        end = datetime.datetime.now()
        print("Time elapsed during the calculation:", end - start)
        
    if ext==1: # superspreaders model
        with open("config_for_grid_search_superspreaders.json") as f:
            config = json.load(f)

        num_icus=int(200/scale)#200 baseline

        num_nodes0=int((1000000-extpop)/scale)#1,000,000 baseline
        num_clusters0=int(num_nodes0/cluster_size)

        num_nodes1=int(extpop/scale)#1,000,000 baseline
        num_clusters1=int(num_nodes1/cluster_size)

        num_days=int(1200/scaledays) #1200 baseline, we want most, even controlled pandemic to end.

        config["graph_generation"][0]["num_people_per_cluster"]=cluster_size
        config["graph_generation"][0]["num_clusters"]=num_clusters0

        config["graph_generation"][1]["num_people_per_cluster"]=cluster_size
        config["graph_generation"][1]["num_clusters"]=num_clusters1

        config["graph_generation"][0]["category_ratios"][0]=int(1000000*0.72-num_nodes1) # diminish <60 because some are superspreaders
        config["graph_generation"][0]["category_ratios"][1]=int(1000000*0.28)
        
        config["simulation"]["stopping_conditions"]["num_days"]=num_days
        config["simulation"]["num_icus"]=num_icus
        config["simulation"]["mu"] = mu
        config["simulation"]["k_trip"] = k
        config["simulation"]["events"][0]["update_params"]["prob_s_to_i"]=[el * scale for el in config["simulation"]["events"][0]["update_params"]["prob_s_to_i"]]

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
                config["simulation"]["initial_params"][0]["prob_goes_on_trip"] = p1
                config["simulation"]["initial_params"][1]["prob_goes_on_trip"] = p1
                config["simulation"]["initial_params"][0]["prob_c_neighbour_trip_candidate"] = p2
                config["simulation"]["initial_params"][1]["prob_c_neighbour_trip_candidate"] = p2

                config_file_name = "tmp/Superspreader_tmp_config_cbg_{}_{}_{}_{}_extpop{}.json".format(
                        k, mu, cluster_size, seed, extpop)
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

        with open("outputs/superspreaders_model/crit_bound_search_extpop{}_k_trip{}_mu{}_cluster_size{}.json".format(
                parsed.extpop, parsed.k, parsed.mu, cluster_size), "w") as f:
            json.dump(bounds_dict, f, indent=4)

        devnull.close()

        end = datetime.datetime.now()
        print("Time elapsed during the calculation:", end - start)
        
    if ext==2: # domovi model
        with open("config_for_grid_search_domovi.json") as f:
            config = json.load(f)

        num_icus=int(200/scale)#200 baseline

        num_nodes0=int((1000000-extpop)/scale)#1,000,000 baseline
        num_clusters0=int(num_nodes0/cluster_size)

        num_nodes1=int(extpop/scale)#1,000,000 baseline
        num_clusters1=int(num_nodes1/cluster_size)

        num_days=int(1200/scaledays) #1200 baseline, we want most, even controlled pandemic to end.

        config["graph_generation"][0]["num_people_per_cluster"]=cluster_size
        config["graph_generation"][0]["num_clusters"]=num_clusters0

        config["graph_generation"][1]["num_people_per_cluster"]=cluster_size
        config["graph_generation"][1]["num_clusters"]=num_clusters1

        config["graph_generation"][0]["category_ratios"][0]=int(1000000*0.72) 
        config["graph_generation"][0]["category_ratios"][1]=int(1000000*0.28-num_nodes1) # diminish >=60 because some are in domovi
        
        config["simulation"]["stopping_conditions"]["num_days"]=num_days
        config["simulation"]["num_icus"]=num_icus
        config["simulation"]["mu"] = mu
        config["simulation"]["k_trip"] = k
        config["simulation"]["events"][0]["update_params"]["prob_s_to_i"]=[el * scale for el in config["simulation"]["events"][0]["update_params"]["prob_s_to_i"]]

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
                config["simulation"]["initial_params"][0]["prob_goes_on_trip"] = p1
                config["simulation"]["initial_params"][1]["prob_goes_on_trip"] = p1
                config["simulation"]["initial_params"][0]["prob_c_neighbour_trip_candidate"] = p2
                config["simulation"]["initial_params"][1]["prob_c_neighbour_trip_candidate"] = p2

                # >=60 in domovi 10 times more likely to be in quarantine
                config["simulation"]["initial_params"][2]["prob_goes_on_trip"] = p1/10 
                config["simulation"]["initial_params"][2]["prob_c_neighbour_trip_candidate"] = p2/10

                config_file_name = "tmp/Domovi_tmp_config_cbg_{}_{}_{}_{}_extpop{}.json".format(
                        k, mu, cluster_size, seed, extpop)
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

        with open("outputs/domovi_model/crit_bound_search_extpop{}_k_trip{}_mu{}_cluster_size{}.json".format(
                parsed.extpop, parsed.k, parsed.mu, cluster_size), "w") as f:
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
