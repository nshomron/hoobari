import numpy as np
from sys import argv

MAX_VALUE=1000

def calculatePhred(joint_probabilities):
    index = tuple(joint_probabilities)
    joint_probabilities -= np.min(joint_probabilities)
    joint_probabilities = np.exp(joint_probabilities)
    sum_probabilities = joint_probabilities.sum()
    posteriors = joint_probabilities / sum_probabilities
    max_probability = np.max(posteriors)

    return ((np.log10(1-max_probability))*(-10))

def randCalc(arr):

    # Reset seed
    np.random.seed(None)
    arr = [arr]
    arr.extend(np.random.random_sample(2) * MAX_VALUE)
    return (tuple(arr),calculatePhred(arr))

# Testing
import pandas as pd
from multiprocessing import Pool
from collections import OrderedDict

pool = Pool(int(argv[2]))

# Construct testing array
probs = [0.0]*(int(argv[1]))

# Process array
results = pool.map(randCalc, probs)

res=OrderedDict({})
for key, value in results:
    res[key] = value

# Plotting results
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt

# Parameter to check
diffs=[]
maxs=[]
for key in res.keys():
    diffs.append(np.fabs(key[1] - key[2]))
    maxs.append(np.max(key))
print(list(res.values()))

# Plot
colors=list(res.values)
plt.scatter(max, diff, c=)
ax0.plot(np.array(maxs),np.array(list(res.values())))
# ax0.ylabel("PHRED")
ax0.title("Max values")
ax1.plot(np.array(diffs),np.array(list(res.values())))
ax1.title("Diffrence betweeen values")
fig.tight_layout()
fig.savefig('ha.png')
