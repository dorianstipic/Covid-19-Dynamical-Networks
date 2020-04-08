import numpy as np
import matplotlib.pyplot as plt
import json
import sys

# Expecting args: config.json output.json

with open(sys.argv[1]) as f:
    config = json.load(f)
with open(sys.argv[2]) as f:
    data = json.load(f)
data = data["stats"]

# "mild" is better name than "confirmed" since people in "icu" state are also
# confirmed cases in the real world
data["mild"] = data["confirmed"]
del data["confirmed"]

for state in data:
    plt.plot(np.arange(len(data[state])), data[state], label=state)
for event in config["simulation"]["events"]:
    plt.axvline(x=event["day"], label=event["label"])

plt.legend()
plt.show()
