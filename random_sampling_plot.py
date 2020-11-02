import numpy as np
import matplotlib.pyplot as plt
import json
from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument('output_files', metavar='output_file', nargs='+')
parsed = parser.parse_args()

prob_transmissions = []
prob_i_to_ics = []
deads = []

for output_file in parsed.output_files:
    with open(output_file) as f:
        output = json.load(f)
        simulation_conf = output["config"]["simulation"]
        prob_transmissions.append(simulation_conf["prob_transmission"])
        prob_i_to_ics.append(simulation_conf["initial_params"][0]["prob_i_to_ic"])
        deads.append(output["stats"]["dead"][-1])

print(json.dumps({
    "prob_i_to_ic": prob_i_to_ics,
    "dead": deads,
    "prob_transmission": prob_transmissions,
}))
#plt.scatter(prob_transmissions, prob_i_to_ics)
#plt.plot(prob_transmissions, deads)
#pandemics = list(filter(lambda x: x > 1000, deads))
#plt.hist(deads, bins=10)
#plt.hist(pandemics, bins=10)
#plt.show()
