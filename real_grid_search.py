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

# LOGICAL STRUCTURE OF THE GRID SEARCH: (keys:elements) is a dictionary , [] is a list.

        # cluster_dict : p1_dict : p2_dict : [seed_dict : [list of 8 elements], ratio_of_success_seeds]

        # cluster_dict is just logical, implementationally non necessary (we save the whole below structure using the cluster value)
            #(we parallelize therefore always only one value cluster.)
        # p1_dict is a dictionary which contains all values p1
        # p2_dict is a dictionary which contains all values p2
        # for each value of p2 we have a list of 2 elements a seed_dict and ratio_of_success_seeds
        # seed_dict is a dictionary which contains all values seed
        # for each value of seed we have a list of 5 elements (relevant info we want to keep track of) :
            #1 first day of infectious case,2 length of pandemic (first infectious, last mild case),
            #3 peak of pandemic (max (infectious+mild+ic) ) 4 peak of pandemic known cases,
            #5 total corona deaths, 6 total no_corona deaths
            #7 1st day system overload,8 n days in system overload.
        # ratio_of_success_seeds is simply the ratio of seeds for which the system is never in fail mode.

def model_cluster_trip(config_file_name, seed, devnull):
    if sys.version_info > (3, 0):
        return subprocess.run(["./model_cluster_trip_v2", config_file_name, str(seed)],
                stdout=subprocess.PIPE, stderr=devnull).stdout
    else:
        p = subprocess.Popen(["./model_cluster_trip_v2", config_file_name, str(seed)],
                stdout=subprocess.PIPE, stderr=devnull)
        return "".join(p.stdout.readlines())

def graph_generation(baseline_icu,baseline_nodes,baseline_days,scale,scaledays,ext,cluster_size,mu,k,extpop):
    num_icus=baseline_icu // scale
    num_days=baseline_days // scaledays 

    config_filenames = [
            "config_for_grid_search.json",
            "config_for_grid_search_superspreaders.json",
            "config_for_grid_search_domovi.json"]
    with open(config_filenames[ext]) as f:
        config = json.load(f)

    if ext==0:
        num_nodes=baseline_nodes // scale
        num_clusters=num_nodes // cluster_size

        config["graph_generation"][0]["num_people_per_cluster"]=cluster_size
        config["graph_generation"][0]["num_clusters"]=num_clusters
        
    if ext==1 or ext==2: # superspreaders model and/or domovi model
        num_nodes0=(baseline_nodes-extpop) // scale
        num_clusters0=num_nodes0 // cluster_size

        num_nodes1=extpop // scale
        num_clusters1=num_nodes1 // cluster_size

        config["graph_generation"][0]["num_people_per_cluster"]=cluster_size
        config["graph_generation"][0]["num_clusters"]=num_clusters0

        config["graph_generation"][1]["num_people_per_cluster"]=cluster_size
        config["graph_generation"][1]["num_clusters"]=num_clusters1
        
        a=config["graph_generation"][0]["category_ratios"][0]
        b=config["graph_generation"][0]["category_ratios"][1]
        
        if ext==1:
            config["graph_generation"][0]["category_ratios"][0]=int(max(baseline_nodes*a/(a+b)/scale-num_nodes1,0)) # diminish <60 because some are superspreaders
            config["graph_generation"][0]["category_ratios"][1]=int(max(baseline_nodes*b/(a+b)/scale,0))
        if ext==2:
            config["graph_generation"][0]["category_ratios"][0]=int(max(baseline_nodes*a/(a+b)/scale,0)) 
            config["graph_generation"][0]["category_ratios"][1]=int(max(baseline_nodes*b/(a+b)/scale-num_nodes1,0)) # diminish >=60 because some are in domovi
            
    config["simulation"]["stopping_conditions"]["num_days"]=num_days
    config["simulation"]["num_icus"]=num_icus
    config["simulation"]["mu"] = mu
    config["simulation"]["k_trip"] = k
    config["simulation"]["events"][0]["update_params"]["prob_s_to_i"]=[el * scale for el in config["simulation"]["events"][0]["update_params"]["prob_s_to_i"]]
    
    return config
    
def grid_search_parameters(config,p1,p2,ext,k,mu,cluster_size,seed,extpop):
    config["simulation"]["initial_params"][0]["prob_goes_on_trip"] = p1
    config["simulation"]["initial_params"][1]["prob_goes_on_trip"] = p1
    config["simulation"]["initial_params"][0]["prob_c_neighbour_trip_candidate"] = p2
    config["simulation"]["initial_params"][1]["prob_c_neighbour_trip_candidate"] = p2
    
    baseline_file_name = "_tmp_config_rgs{}_{}_{}_{}"
    extra_param = "_extpop{}"
    if ext==0:
        config_file_name = ("tmp/Base" + baseline_file_name + ".json").format(
            k, mu, cluster_size, seed)
    if ext==1:
        config_file_name = ("tmp/Superspreaders" + baseline_file_name + extra_param + ".json").format(
                k, mu, cluster_size, seed, extpop)
    if ext==2:
        # >=60 in domovi 10 times more likely to be in quarantine and 10 times more likely not to go on trip.
        config["simulation"]["initial_params"][2]["prob_goes_on_trip"] = p1/10 
        config["simulation"]["initial_params"][2]["prob_c_neighbour_trip_candidate"] = p2/10

        config_file_name = ("tmp/Domovi" + baseline_file_name + extra_param + ".json").format(
                k, mu, cluster_size, seed, extpop)
    return [config, config_file_name]
    
def save_json(ext,p1_dict, k, mu, cluster_size, extpop):
    baseline_file_name1 = "_real_grid_search"
    baseline_file_name2 = "_k_trip{}_mu{}_cluster_size{}.json"
    extra_param = "_extpop{}"
    if ext==0:
        with open(("outputs/base_model/Base" + baseline_file_name1 + baseline_file_name2).format(
                k, mu, cluster_size), "w") as f:
            json.dump(p1_dict, f, indent=4)
    if ext==1:
        with open(("outputs/superspreaders_model/Superspreaders" + baseline_file_name1 + extra_param + baseline_file_name2).format(
                extpop, k, mu, cluster_size), "w") as f:
            json.dump(p1_dict, f, indent=4)
    if ext==2:
        with open(("outputs/domovi_model/Domovi" + baseline_file_name1 + extra_param + baseline_file_name2).format(
                extpop, k, mu, cluster_size), "w") as f:
            json.dump(p1_dict, f, indent=4)


def f(cluster_size):

    ext=parsed.ext
    mu=parsed.mu
    k=parsed.k
    extpop=parsed.extpop
    scale=1
    scaledays=1

    config=graph_generation(200,1000000,1200,scale,scaledays,ext,cluster_size,mu,k,extpop)

    h=5 #put 1 for very precise grid
    ptrip=np.arange(0,1.00001,0.01*h)
    pdisobedient=np.arange(0,1.00001,0.01*h)
    seeds=np.arange(0,1,1) 

    p1_dict={}
    start = datetime.datetime.now()
    devnull = open(os.devnull, 'w')
    for p1 in ptrip:
        p2_dict={}
        for p2 in pdisobedient:
            list_2len=[]
            seed_dict={}
            ratio_succ=-1
            succ=0
            for seed in seeds:
                config_list=grid_search_parameters(config,p1,p2,ext,k,mu,cluster_size,seed,extpop)
                config=config_list[0]
                config_file_name=config_list[1]
                
                with open(config_file_name, "w") as f:
                    json.dump(config, f, indent=4)

                print("Running model with params: cluster_size = {:.3f}".format(cluster_size),
                      ", prob_goes_on_trip = {:.3f}".format(p1),
                      ", prob_c_neighbour_trip_candidate = {:.3f}".format(p2),
                      "seed = {}".format(seed), file=sys.stderr)
                stdout = model_cluster_trip(config_file_name, seed, devnull)

                output=json.loads(stdout)
                os.remove(config_file_name)
                try:
                    beginning_pandemic = next(x for x, val in enumerate(output["stats"]["infectious"]) if val > 0)
                except StopIteration:
                    beginning_pandemic = -1
                    end_pandemic = -1
                    len_pandemic = -1
                else:
                    beginning_pandemic = next(x for x, val in enumerate(output["stats"]["infectious"]) if val > 0)
                    end_pandemic = len((output["stats"]["confirmed"])) - 1 - next(x for x, val
                                                    in enumerate(reversed(output["stats"]["confirmed"])) if val > 0 )
                    len_pandemic = end_pandemic - beginning_pandemic + 1

                peak_corona_total = max([sum(x) for x in zip( (output["stats"]["infectious"]),(output["stats"]["confirmed"]),
                                                             (output["stats"]["icu"]))])
                peak_corona_system_load = max([sum(x) for x in zip( (output["stats"]["confirmed"]), (output["stats"]["icu"]))])
                corona_deaths = output["stats"]["dead"][-1]
                no_corona_deaths = output["stats"]["nocorona_dead"][-1]
                
                total_immune = output["stats"]["immune"][-1]

                n_days_icu_overflow = output["num_days_icu_overflow"]
                first_day_icu_overflow = output["first_day_icu_overflow"]

                if n_days_icu_overflow==0:
                    succ=succ+1
                
                seed_dict[str(seed)]=[beginning_pandemic, len_pandemic, peak_corona_total, peak_corona_system_load, 
                                        corona_deaths, no_corona_deaths, n_days_icu_overflow, first_day_icu_overflow, total_immune]

            ratio_succ=float(succ/len(seeds))

            list_2len.append(seed_dict)
            list_2len.append(ratio_succ)

            p2_dict[str(p2)]=list_2len

        p1_dict[str(p1)]=p2_dict

    save_json(ext,p1_dict, k, mu, cluster_size, extpop)

    devnull.close()

    end = datetime.datetime.now()
    print("Time elapsed during the calculation:", end - start)
    

if parsed.num_processes == 1:
    print("Running ONE process", file=sys.stderr)
    for cluster_size in parsed.cluster_sizes:
        f(cluster_size)
else:
    print("Running {} processes".format(parsed.num_processes), file=sys.stderr)
    # This won't work on isabella, because python2 has different API for Pool
    # object. However, that is fine since we run parallelize our code with SGE
    # scripts.
    with Pool(parsed.num_processes) as p:
        p.map(f, parsed.cluster_sizes)
