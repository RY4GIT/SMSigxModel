import multiprocessing as mp
print("Number of processors: ", mp.cpu_count())
import numpy as np
from time import time

# Prepare data
np.random.RandomState(100)
arr = np.random.randint(0, 10, size=[2000, 5])
data = arr.tolist()
print(data[:5])

# https://www.machinelearningplus.com/python/parallel-processing-python/
results = []
pool = mp.Pool(mp.cpu_count())

def test_run(a, b):
    c = a+b
    return c

def collect_result(result):
    global results
    results.append(result)
