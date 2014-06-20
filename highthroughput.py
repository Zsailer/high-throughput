from matplotlib import *
from itertools import product as itp
from numpy import *
from math import *
import random as r

class HighThroughput(object):

    def __init__(self, sequences, system_size, sequencing_size):
        """
        Modules necessary to run a simulation of high through-put experiment

        Parameters:
        ----------
        n: int
            Number of sites in each sequences.
        sequencing_size: int
            Number of protein that can be experimentally sequenced.
        """
        self.sequences = sequences
        self.system_size = system_size
        self.sequencing_size = sequencing_size

    def build_system(self, sequences, system_size=None):
        """
        Creates a system of sequences with randomly generated weights. 
        
        Parameters:
        ----------
        sequences: dict
            Sequence-fitness (key, value) pair for the system.
        system_size: int
            The number of sequences to include in the initial system (before screening).
            
        Return:
        ------
        sys: list
            A lists of all sequences in the initial system.
        frequencies: dict
            Sequence-counts (key,value) pair that describes the number of times
            the sequence appears in the list of sequences.
        """
        if system_size is None:
            system_size = self.system_size
        
        seq = sequences.keys()
        system = list()
        frequencies = dict()
        for s in seq:
            frequencies[s] = 0
            
        for i in range(system_size):
            s = seq[r.randint(0,len(sequences)-1)]
            system.append(s)
            frequencies[s] = frequencies[s] + 1

        return system, frequencies
        
    def sequencing(self, sequences, system, sequencing_size=None):
        """
        """
        if sequencing_size is None:
            sequencing_size = self.sequencing_size
            
        sampling_freq = dict()
        for s in sequences.keys():
            sampling_freq[s] = 0
        
        print "Beginning Sequencing..."
        
        for i in range(int(sequencing_size)):
            s = system[r.randint(0,len(system)-1)]
            sampling_freq[s] = sampling_freq[s] + 1
            
        print "Completed Sequencing."
            
        return sampling_freq
        
    def screening(self, sequences, system):
        """
        Run a screening simulation on the initial set of sequences.
        """
        # Initialize a dict to hold sequence-count (key-value) pairs
        frequencies = dict()
        for s in sequences.keys():
            frequencies[s] = 0
            
        print "Beginning Screening..."
        # Track sequences that survive screening
        screened_system = list()
        # Screening is done using a random number generator
        for i in range(len(system)):
            if sequences[system[i]] >= r.random():
                screened_system.append(system[i])
                frequencies[system[i]] = frequencies[system[i]] + 1
        
        print "Completed Screening."
        
        return screened_system, frequencies
        
    def run_experiment(self, sequences, system_size=None, sequencing_size=None):
        """
        """
        if system_size is None:
            system_size = self.system_size
        if sequencing_size is None:
            sequencing_size = self.sequencing_size
        
        # Build the initial system of sequences    
        sys, freq = self.build_system(sequences, system_size)
        
        # Sequence this system with our limited sequencing size
        samp0 = self.sequencing(sequences, sys, sequencing_size)
        
        # Screen the system and store new system
        sc_sys, sc_freq = self.screening(sequences, sys)
        
        # Sequence the screened system
        samp1 = self.sequencing(sequences, sc_sys, sequencing_size)
        
        return sys, freq, samp0, sc_sys, sc_freq, samp1