import multiprocessing as mp
print("Number of processors: ", mp.cpu_count())
import numpy as np
from time import time

# Prepare data
np.random.RandomState(100)
arr = np.random.randint(0, 10, size=[10, 5])
data = arr.tolist()
# print(data[:5])

# https://www.machinelearningplus.com/python/parallel-processing-python/


def collect_result(result):
    global results
    results.append(result)

def howmany_within_range(row, minimum=4, maximum=8):
    """Returns how many numbers lie within `maximum` and `minimum` in a given `row`"""
    count = 0
    for n in row:
        if minimum <= n <= maximum:
            count = count + 1
    return count

def print_arg():
    """Returns how many numbers lie within `maximum` and `minimum` in a given `row`"""
    return print('test')

if __name__ == '__main__':
    # freeze_support()
    """
    results = []
    pool = mp.Pool(mp.cpu_count())
    results = pool.map(howmany_within_range)
    pool.close()
    print(results[:10])
    """

    process1 = mp.Process(target=print_arg)
    process2 = mp.Process(target=print_arg)
    process1.start()
    process2.start()
