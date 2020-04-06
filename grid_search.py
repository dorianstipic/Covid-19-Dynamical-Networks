import numpy as np
import json
import subprocess
import matplotlib.pyplot as plt
import sys
import os

with open("config_cluster_trip.json") as f:
    config = json.load(f)

mu_range = np.arange(0.5, 1.5, 0.1)

devnull = open(os.devnull, 'w')

results = []
for mu in mu_range:
    config["simulation"]["mu"] = mu
    with open("tmp_config.json", "w") as f:
        json.dump(config, f, indent=4)
    print("Running model with params: mu={}".format(mu), file=sys.stderr)
    completed = subprocess.run(["./model_cluster_trip", "tmp_config.json", "0"],
            stdout=subprocess.PIPE, stderr=devnull)
    results.append(json.loads(completed.stdout))

devnull.close()

os.remove("tmp_config.json")

deads = []
for result in results:
    deads.append(result["dead"][-1])

plt.plot(mu_range, deads, label="dead")
plt.show()
