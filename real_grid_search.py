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

#scale=200
#config["simulation"]["num_icus"]=config["simulation"]["num_icus"]/scale
#config["graph_generation"]["num_clusters"]=config["graph_generation"]["num_clusters"]/scale
#config["simulation"]["events"][0]["prob_s_to_i"]=[el * scale for el in config["simulation"]["events"][0]["prob_s_to_i"]]

mu_range=np.concatenate((np.arange(0,1,0.2), np.arange(1,3,0.4),np.arange(3,5,0.5),np.arange(5,10.5,1) )
                        , axis=None)

def f(mu):
    with open("config_for_grid_search.json") as f:
        config = json.load(f)

    h=5
    seeds=np.arange(0,100,10)

    #mu_range=np.arange(0,0.6,0.2)
    ptrip=np.arange(0,1,0.01*h)
    pdisobedient=np.arange(0,1,0.01*h)

    p1=0 # ptrip
    p2=0 # pdisobedient

    step=0.01*h # ovo nam je step za p2 u while-u

    p1_dict={} # ovo ce biti dictionary, key su vrijednosti p1
    p2=1
    start = datetime.datetime.now()
    devnull = open(os.devnull, 'w')
    for p1 in ptrip:
        p2_list=[] # ovo ce biti lista. Tu spremamo sve vrijednosti od p2 (za svaki seed dakle) koje ce pripadati jednom p1.
        for seed in seeds:
            while p2>=0:
                config["simulation"]["mu"] = mu
                config["simulation"]["prob_goes_on_trip"] = p1
                config["simulation"]["initial_params"][0]["prob_c_neighbour_trip_candidate"] = p2
                config["simulation"]["initial_params"][1]["prob_c_neighbour_trip_candidate"] = p2
                config_file_name = "tmp/tmp_config{}.json".format(mu)
                with open(config_file_name, "w") as f:
                    json.dump(config, f, indent=4)

                print("Running model with params: mu = {:.3f}".format(mu), ", prob_goes_on_trip = {:.3f}".format(p1),
                      ", prob_c_neighbour_trip_candidate = {:.3f}".format(p2), "seed = {}".format(seed), file=sys.stderr)
                rez = subprocess.run(["./model_cluster_trip", config_file_name, str(seed)],
                        stdout=subprocess.PIPE, stderr=devnull)
                os.remove(config_file_name)

                # ako nisi jos uvik nasao rub smanji za step i nastavi.
                if json.loads(rez.stdout)["stopping_condition"]=="icu_overflow":
                    p2=p2-step
                    continue
                # ako si nasao, prirodaj taj p2 listi, za iduci seed povecaj za odredjeni broj koraka.
                else:
                    p2_list.append(p2)
                    ptemp=min(1,p2+10*step) #dignemo ga natrag za neke steppove za iduci SEED.
                    break
            p2=ptemp
        if (len(p2_list)==0 or np.mean(p2_list)<=0):
            break
        else:
            p2=max(p2_list) # za iducu iteraciju, korak udesno od p1 krecemo od max, to je sigurno dovoljno dobro.

        p1_dict[p1]=p2_list # listu od p2, pripajamo dictionaryju koji ima p1 kao key.

    with open("outputs/{0}.json".format(mu), "w") as f:
        json.dump(p1_dict, f, indent=4)

    devnull.close()
    end = datetime.datetime.now()
    print("Time elapsed during the calculation:", end - start)

with Pool(4) as p:
    p.map(f, mu_range)
