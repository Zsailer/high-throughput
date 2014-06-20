from matplotlib import *
import itertools as it
from numpy import *
from math import *
import random as r

class HighThroughputAnalysis(object):
    
    def __init__(self):
        pass
        
    def kmer_counting(self, char, k, frequencies, strict=False):
        """
        Count the frequency of particular kmer sequences
        
        Parameter:
        ---------
        char: 'string'
            A single character, either '0' and '1', that will be counted for frequency
        k: int
            The order to which the frequency is counted
        frequencies: dict
            Includes sequence-frequency (key,pair) for kmer counting
            
        Returns:
        -------
        kmer: dict
            Dictionary with kmer frequencies. A kmer is a string with all combinations
            of char in the sequence.
        """
        if char is not '0' and char is not '1':
            raise NameError('char must be "0" or "1"')
            
        sequences = frequencies.keys()            
        #split_seq = [list(s) for s in sequences]
        
        n = len(sequences[0])
        kmer = dict()
        
        initial = "".join([str(i) for i in range(1,n+1)])        
        kmers = [list(r) for r in it.combinations(initial, k)]
        
        for i in range(len(kmers)):
            contributions = list()
            for s in sequences:
                for index in kmers[i]:
                    if s[int(index)-1] is not char:
                        keep = True
                    else:
                        keep = False
                        break
                if keep is True:
                    contributions.append(frequencies[s])
            kmer["".join(kmers[i])] = sum(contributions)
            
        return kmer


    def k_order_contribution(self, k, freq):
        """
        Compares the freq0/freq1 to k order.
        """
        kmer0 = self.kmer_counting('0', k, freq)        
        kmer1 = self.kmer_counting('1', k, freq)

        kmer_contribution = dict()
        # Calculate the ratio of frequencies f0/f1
        for k in kmer0:
            kmer_contribution[k] = float(kmer1[k])/float(kmer0[k])

        return kmer_contribution


    def k_system_calculation(self, k, freq_before, freq_after):
        """
        Determine all interactions between sites to k order
        
        Parameters:
        ----------
        k: int
            The highest order interaction to calculate.
        freq_before: dict
            A dictionary which holds sequence-frequency (key,value) pair for initial
            comparison (i.e. before a high-throughput screening).
        freq_after: dict
            A dictionary which holds sequence-frequency (key,value) pair for final
            comparison (i.e. after a high-throughput screening).
        
        Return:
        ------
        system: dict
            A dictionary which holds combination-frequency (key,value) pair, where 
            combination is the index of all sites with '1'
        """        
        system = dict()
        for i in range(1,k+1):
            system[str(i)] = self.k_order_contribution(i, freq_before, freq_after)
        return system


    def sequence_frequency_function(self, sequence, kmer_frequencies):
        """
        Enter a sequence and kmer frequency dict, 
        return the frequency contributions for function
        """
            
        k = len(kmer_frequencies)
        sites = list()
        for i in range(len(sequence)):
            if sequence[i] is '1':
                sites.append(str(i+1))
        sites = "".join(sites)
        
        freq_func = dict()
        list_freq = list()
        for i in range(len(sites)):
            kmer_contr = ["".join(list(s)) for s in it.combinations(sites, i+1)]
            freq_func[str(i+1)] = dict()
            for j in range(len(kmer_contr)):
                freq = kmer_frequencies[str(i+1)][kmer_contr[j]]
                freq_func[str(i+1)][kmer_contr[j]] = freq
                list_freq.append(freq)
        
        total_freq = sum(list_freq)
        return freq_func, list_freq, total_freq
