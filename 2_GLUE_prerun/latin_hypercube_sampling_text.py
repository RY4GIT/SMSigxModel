# https://github.com/thouska/spotpy/blob/master/src/spotpy/examples/spot_setup_rosenbrock.py
# https://github.com/thouska/spotpy/blob/master/tutorials/tutorial_rosenbrock.py
# https://github.com/thouska/spotpy/blob/master/src/spotpy/algorithms/lhs.py
# https://github.com/thouska/spotpy/blob/master/src/spotpy/algorithms/_algorithm.py

# -*- coding: utf-8 -*-
"""

"""
# %%

import random
import numpy as np
from spotpy.algorithms import lhs
import numpy as np
from spotpy.parameter import Uniform
import spotpy

np.random.seed(0)


# %%
class mylhs(lhs):
    """
    The Latin Hypercube algorithm generates random parameters from their respective
    distribution functions.

    Copyright (c) 2018 by Tobias Houska
    This file is part of Statistical Parameter Optimization Tool for Python (SPOTPY).
    :author: Tobias Houska

    Ryoko modified this class just to get sampled parameters
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

        ################ CHANGES MADE HERE ################
        # A generator that produces the parameters
        sampled_params = []
        param_generator = ((rep, matrix[rep]) for rep in range(int(repetitions)))
        for _, randompar, _ in self.repeat(param_generator):
            sampled_params.append(randompar)
            print(randompar)
        ###################################################

        return sampled_params


# %%
class spot_setup(object):
    def __init__(self):
        self.dim = 3
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
# Example run
rep = 5
sampler = mylhs(spot_setup())
sampler.sample(rep)
