# https://github.com/thouska/spotpy/blob/master/src/spotpy/examples/spot_setup_rosenbrock.py
# https://github.com/thouska/spotpy/blob/master/tutorials/tutorial_rosenbrock.py
# https://github.com/thouska/spotpy/blob/master/src/spotpy/algorithms/lhs.py
# https://github.com/thouska/spotpy/blob/master/src/spotpy/algorithms/_algorithm.py

# -*- coding: utf-8 -*-
"""
Copyright (c) 2018 by Tobias Houska
This file is part of Statistical Parameter Optimization Tool for Python (SPOTPY).
:author: Tobias Houska
"""
# %%
import random

import numpy as np

from spotpy.algorithms import lhs
import spotpy
from spotpy.describe import describe


# %%
class mylhs(lhs):
    """
    The Latin Hypercube algorithm generates random parameters from their respective
    distribution functions.
    """

    def __init__(self, *args, **kwargs):
        super(lhs, self).__init__(*args, **kwargs)

    def sample(self, repetitions):
        """
        Parameters
        ----------
        repetitions: int
            maximum number of function evaluations allowed during optimization
        """
        self.set_repetiton(repetitions)
        print(
            "Starting the LHS algotrithm with " + str(repetitions) + " repetitions..."
        )
        print("Creating LatinHyperCube Matrix")
        # Get the names of the parameters to analyse
        names = self.parameter()["name"]
        print(names)
        # Define the jump size between the parameter
        segment = 1 / float(repetitions)
        # Get the minimum and maximum value for each parameter from the
        # distribution
        parmin, parmax = (
            self.parameter()["minbound"],
            self.parameter()["maxbound"],
        )

        # Create an matrx to store the parameter sets
        matrix = np.empty((repetitions, len(parmin)))
        # Create the LatinHypercube matrix as given in McKay et al. (1979)
        for i in range(int(repetitions)):
            segmentMin = i * segment
            pointInSegment = segmentMin + (random.random() * segment)
            parset = pointInSegment * (parmax - parmin) + parmin
            matrix[i] = parset
        for i in range(len(names)):
            random.shuffle(matrix[:, i])

        # A generator that produces the parameters
        param_generator = ((rep, matrix[rep]) for rep in range(int(repetitions)))
        for rep, randompar, simulations in self.repeat(param_generator):
            # A function that calculates the fitness of the run and the manages the database
            print(randompar)
            self.postprocessing(rep, randompar, simulations)
        self.final_call()


# %%
"""
Copyright 2015 by Tobias Houska
This file is part of Statistical Parameter Optimization Tool for Python (SPOTPY).
:author: Tobias Houska

This example implements the Rosenbrock function into a SPOTPY class.
"""

import numpy as np

from spotpy.objectivefunctions import rmse
from spotpy.parameter import Uniform


class spot_setup(object):
    def __init__(self, obj_func=None):
        self.dim = 3
        # self.params = [Uniform("a", -32.768, 32.768, 2.5, -20.0)]
        self.params = [
            Uniform("bb", low=2, high=15),
            Uniform("slop", low=0, high=1),
            Uniform("satdk", low=0.001, high=0.002),
        ]

    def parameters(self):
        return spotpy.parameter.generate(self.params)

    def simulation(self, vector):
        return [np.zeros(1)]

    def evaluation(self):
        return [np.zeros(1)]

    def objectivefunction(self, simulation, evaluation, params=None):
        return np.zeros(1)


# %%

# Create samplers for every algorithm:
results = []
rep = 5
timeout = 10  # Given in Seconds
parallel = "seq"
dbformat = "csv"
# %%
sampler = mylhs(
    spot_setup(),
    parallel=parallel,
    dbname="RosenLHS",
    dbformat=dbformat,
    sim_timeout=timeout,
)
# %%
sampler.sample(rep)


# %%
sampler

# %%
spot_setup().parameter()["name"]
# %%
sampler.get_parameters()
# %%
