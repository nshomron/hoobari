import numpy as np
from sys import argv
from decimal import *

getcontext().prec = 45
MAX_VALUE=100

def calculatePhred(joint_probabilities):
    index = tuple(joint_probabilities)
    joint_probabilities -= np.min(joint_probabilities)
    joint_probabilities = np.exp(joint_probabilities)
    sum_probabilities = joint_probabilities.sum()
    posteriors = joint_probabilities / sum_probabilities
    max_probability = np.max(posteriors)

    if max_probability == 1:
        print('ERROR:',index,'product is zero')

    return ((np.log10(1-max_probability))*(-10))

def randCalc(arr):

    # Reset seed
    np.random.seed(None)
    arr = [arr]
    arr.extend(np.random.random_sample(2) * MAX_VALUE)

    # Turn to extreme precision in cases of possible overflow
    if np.fabs(arr[2] - arr[1]) > 35 and 'nodecimal' not in argv:
        return (tuple(arr),calculatePhred(np.array([Decimal(arr[0]),Decimal(arr[1]),Decimal(arr[2])])))
    else:
        return (tuple(arr),calculatePhred(np.array(arr)))

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
import matplotlib.patches as mpc
matplotlib.use('agg')
import matplotlib.pyplot as plt

# Parameter to check
diffs=[]
maxs=[]
for key in res.keys():
    diffs.append(np.fabs(key[1] - key[2]))
    maxs.append(np.max(key))

# Plot
colors=[]
for phred in res.values():
    if np.isnan(float(phred)):
        colors.append('#000000')
    elif np.isinf(float(phred)):
        colors.append('#ff0000')
    else:
        colors.append('#00ff04')
plt.scatter(maxs, diffs, c=colors, marker='.')
plt.xlabel("Max values")
plt.ylabel("Diffrence betweeen values")
plt.xlabel("Max values")
plt.minorticks_on()
plt.legend(handles = [mpc.Patch(color='#000000', label='NaN'),
                      mpc.Patch(color='#ff0000', label='Infinity'),
                      mpc.Patch(color='#00ff04', label='Valid number')])
plt.savefig('ha.png', dpi=200.0)
