from matplotlib import *
import itertools as it
from numpy import *
from math import *
import random as r

class FitnessFunctions(object):
    
    def __init__(self):
        pass
        
    def normal_distribution(self, x, center, sigma):
        #sigma squared is the variance
        if sigma is None:
                sigma = 2
        return exp(-((x-center)**2)/(2*sigma**2))
        
    def gumbel_distribution(self):
        return exp(-exp(-x) - x)
        
    def nk_model(self, n, k, distribution, *args):
        """Returns a k"""
        # generate nk model fitness table
        nk_table = dict()
        interactions = ["".join(r) for r in it.product('01', repeat=k)]
        for s in interactions:
            m = s.count('1')
            f = distribution(m,k,.5)
            nk_table[s] = f

        sequences = [list("".join(r)) for r in it.product('01', repeat=n)]

        l = len(sequences[0])
        fitness = dict()
        neighbor = int(k/2)
        for i in range(len(sequences)):
            f_total = 0
            for j in range(l):
                if j-neighbor < 0:
                    pre = sequences[i][-neighbor:]
                    post = sequences[i][j:neighbor+j+1]
                    f = "".join(pre) + "".join(post)
                elif j+neighbor > l-1:
                    pre = sequences[i][j-neighbor:j+1]
                    post = sequences[i][0:neighbor]
                    f = "".join(pre) + "".join(post)
                else:
                    f = "".join(sequences[i][j-neighbor:j+neighbor+1])
                f_total += nk_table[f]
            fitness["".join(sequences[i])] = f_total/l
        
        return sequences, fitness
        
        
    def random_distribution(self, n=None):
        """
        Create a library of sequences with n sites and assign
        random weights to each sequence.

        Parameters:
        ----------
        n: int
            Number os sites in each sequences.

        Returns:
        -------
        sequences: dict
            Sequence (key) and weight(value) pair
        frequencies: dict
            Sequence (key) and initial # of counts == 0 (value)
        """
        if n is None:
            n = self.n
        # Create a reference string of all zeros, with len(n)
        
        sequences = dict()
        frequencies = dict()
        for seq in itp("01", repeat=n):
            s = "".join(seq)
            sequences["".join(seq)] = r.random()
            frequencies[s] = 0
        return sequences, frequencies