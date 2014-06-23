from matplotlib import *
import itertools as it
from numpy import *
from math import *
import random as r

class FitnessFunctions(object):
    
    def __init__(self):
        pass
        
    def normal_distribution(self, x, center, sigma):
        """ Returns a normal distribution """
        #sigma squared is the variance
        height = .2
        if sigma is None:
                sigma = 2
        if isinstance(x,(list,ndarray)):
            y = list()
            for i in x:
                y.append(height*exp(-((i-center)**2)/(2*float(sigma)**2)))
            return y
        else:
            return height*exp(-((x-center)**2)/(2*float(sigma)**2))
        
    def gumbel_distribution(self, x, center):
        """ Returns a gumbel_distribution """
        if isinstance(x,(list,ndarray)):
            y = list()
            for i in x:
                y.append(exp(-exp(-i) - i))
            return y
        else:
            return exp(-exp(-x) - x)

    def nk_fitness(self, n, k, distribution):
        """ Returns an nk fitness distribution """
        # generate nk model fitness table
        nk_table = dict()
        interactions = ["".join(r) for r in it.product('01', repeat=k)]
        for s in interactions:
            m = s.count('1')
            f = distribution(m,k,.4)
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
        
    def bias_interactions(self, interactions, fitness):
        """ Add bias to particular interactions """

        for interaction in interactions:
            bias = list()
            for sequence in fitness:
                for i in interaction:
                    if sequence[int(i)-1] is '1':
                        keep = True
                    else:
                        keep = False
                        break
                if keep is True:
                    bias.append(sequence)
            for b in bias:
                fitness[b] = (1-fitness[b])*interactions[interaction] + fitness[b]
                
        return fitness
                
        