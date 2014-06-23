from matplotlib import *
import itertools as it
from numpy import *
from math import *
import random as r
import heapq as h

class HighThroughputAnalysis(object):
    
    def __init__(self):
        pass
        
    def kmer_counting(self, char, kmer, frequencies):
        """ Count the frequency of particular kmer """
        
        if char is not '0' and char is not '1':
            raise NameError('char must be "0" or "1"')
        
        sequences = frequencies.keys()            
        contributions = list()
        for s in sequences:
            for index in kmer:
                if s[int(index)-1] is char:
                    keep = True
                else:
                    keep = False
                    break
            if keep is True:
                contributions.append(frequencies[s])
        freq_total = sum(contributions)
        
        return freq_total
            
    
    def k_order_counting(self, char, k, frequencies):
        """
        Count the frequency of all k_order sequences
        
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
            kmer["".join(kmers[i])] = self.kmer_counting(char, kmers[i], frequencies)
            
        return kmer


    def k_order_contribution(self, k, freq):
        """
        Compares the freq0/freq1 to k order.
        """
        kmer0 = self.k_order_counting('0', k, freq)        
        kmer1 = self.k_order_counting('1', k, freq)

        kmer_contribution = dict()
        # Calculate the ratio of frequencies f0/f1
        for k in kmer0:
            kmer_contribution[k] = float(kmer1[k])#/float(kmer0[k])*normalization

        return kmer_contribution


    def k_system_calculation(self, k, freq):
        """
        Determine all interactions between sites to k order
        
        Parameters:
        ----------
        k: int
            The highest order interaction to calculate.
        freq: dict
            A dictionary which holds sequence-frequency (key,value) pair for initial
            comparison (i.e. before a high-throughput screening).
        
        Return:
        ------
        system: dict
            A dictionary which holds combination-frequency (key,value) pair, where 
            combination is the index of all sites with '1'
        """        
        system = dict()
        for i in range(1,k+1):
            system[str(i)] = self.k_order_contribution(i, freq)
            
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
     
    def site_interactions(self, k, site, freq):
        """ Specify a site and return all kmer interactions """
        kmer_freqs = self.k_system_calculation(k, freq)
        monomers = kmer_freqs['1']
        
        contributions = dict()
        for i in range(2, k+1):
            kmers = kmer_freqs[str(i)]
            contributions[str(i)] = dict()
            for key in kmers.keys():
                if str(site) in key:
                    contributions[str(i)][key] = float(kmers[key])/monomers[str(i)]
        
        return contributions
        
    def strongest_site_interactions(self, top_tier, k, site, freq):
        """ Find the strongest (top-tier fraction) kmer interations """
        site_stuff = self.site_interactions(k, site, freq)
        
        strongest = dict()
        for i in range(2, len(site_stuff)+1):
            strongest[str(i)] = dict()
            top = int(len(site_stuff[str(i)])*top_tier)
            vals = h.nlargest(top, site_stuff[str(i)].values())
            keys = [site_stuff[str(i)].keys()[site_stuff[str(i)].values().index(v)] for v in vals]
            for j in range(len(keys)):
                strongest[str(i)][keys[j]] = vals[j]
                
        return strongest

