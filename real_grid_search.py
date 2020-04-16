# Autor: Dorian
import numpy as np
import matplotlib.pyplot as plt
import json
import sys
from multiprocessing import Pool
import datetime

from matplotlib import rc

import subprocess
import sys
import os

from mpi4py import MPI

#mu_range=np.concatenate((np.arange(20,100,20), np.arange(1,3,0.4),np.arange(3,5,0.5),np.arange(5,10.5,1) )
#                        , axis=None)
#mu_range=np.concatenate((np.arange(20,100,20), np.arange(100,1000,100),np.arange(1000,5000,500),np.arange(5000,10005,1000) )
#                        , axis=None)
#mu_range=np.concatenate((np.arange(1250.,5000.,500),np.arange(12500.,20005.,2500) ), axis=None)

def f(cluster_size):
    with open("config_for_grid_search.json") as f:
        config = json.load(f)

    mu=5 #FIX mu to 5, was already more extensively analyzed than other mu. #mu=5000 is also somewhat interesting

    scale=1 #1 if real simul
    num_icus=int(200/scale)#200 baseline
    num_nodes=int(1000000/scale)#1,000,000 baseline
    num_clusters=int(num_nodes/cluster_size)

    scaledays=1
    num_days=int(1200/scaledays) #1200 baseline, we want most, even controlled pandemic to end.

    config["graph_generation"]["num_people_per_cluster"]=cluster_size
    config["graph_generation"]["num_clusters"]=num_clusters
    config["simulation"]["stopping_conditions"]["num_days"]=num_days
    config["simulation"]["num_icus"]=num_icus
    config["simulation"]["mu"] = mu
    config["simulation"]["events"][0]["prob_s_to_i"]=[el * scale for el in config["simulation"]["events"][0]["prob_s_to_i"]]

    h=5 #1 for precise grid.
    ptrip=np.arange(0,1.00001,0.01*h)
    pdisobedient=np.arange(0,1.00001,0.01*h)
    seeds=np.arange(0,1,1) #step 1 - 100 possible


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
                list_8len=[]
                config["simulation"]["prob_goes_on_trip"] = p1
                config["simulation"]["initial_params"][0]["prob_c_neighbour_trip_candidate"] = p2
                config["simulation"]["initial_params"][1]["prob_c_neighbour_trip_candidate"] = p2

                #NOT necessary?
                config_file_name = "tmp/tmp_config{}.json".format(cluster_size)
                with open(config_file_name, "w") as f:
                    json.dump(config, f, indent=4)

                print("Running model with params: cluster_size = {:.3f}".format(cluster_size),
                      ", prob_goes_on_trip = {:.3f}".format(p1),
                      ", prob_c_neighbour_trip_candidate = {:.3f}".format(p2),
                      "seed = {}".format(seed), file=sys.stderr)
                rez = subprocess.run(["./model_cluster_trip", config_file_name, str(seed)],
                        stdout=subprocess.PIPE, stderr=devnull)
                os.remove(config_file_name)

                # chance of StopIteration is ((0.999999)^{10^6})^30 = 10^{-13}.
                output=json.loads(rez.stdout)
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

                n_days_icu_overflow=output["num_days_icu_overflow"]
                first_day_icu_overflow=output["first_day_icu_overflow"]

                if n_days_icu_overflow==0:
                    succ=succ+1

                list_8len.append(beginning_pandemic)
                list_8len.append(len_pandemic)
                list_8len.append(peak_corona_total)
                list_8len.append(peak_corona_system_load)
                list_8len.append(corona_deaths)
                list_8len.append(no_corona_deaths)
                list_8len.append(n_days_icu_overflow)
                list_8len.append(first_day_icu_overflow)

                seed_dict[str(seed)]=list_8len

            ratio_succ=float(succ/len(seeds))

            list_2len.append(seed_dict)
            list_2len.append(ratio_succ)

            p2_dict[str(p2)]=list_2len

        p1_dict[str(p1)]=p2_dict



    with open("outputs/{0}.json".format(cluster_size), "w") as f:
        json.dump(p1_dict, f, indent=4)

    devnull.close()

    end = datetime.datetime.now()
    print("Time elapsed during the calculation:", end - start)


#if __name__ ==  '__main__':
# num_processors = 4
# p=Pool(processes = num_processors)
# p.map(f,mu_range)

comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()

if rank == 0:
    cluster_dim_range = np.arange(1.0, size + 0.1, 1)
else:
    cluster_dim_range = None

cluster_size = comm.scatter(cluster_dim_range, root=0)
f(cluster_size)
