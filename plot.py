import numpy as np
import matplotlib.pyplot as plt
import json
import sys

# Expecting args: config.json output.json

with open(sys.argv[1]) as f:
    config = json.load(f)
with open(sys.argv[2]) as f:
    data = json.load(f)

for state in data:
    plt.plot(np.arange(len(data[state])), data[state], label=state)
for event in config["simulation"]["events"]:
    plt.axvline(x=event["day"], label=event["label"])

plt.legend()
plt.show()
