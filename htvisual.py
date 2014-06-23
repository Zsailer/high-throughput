from matplotlib import pyplot as plt
import numpy as np
import scipy as sp
import mpld3

class HighThroughputVisualization(object):
    
    def __init__(self):
        self.title_size = 20
        self.x_label_size = 16
        self.y_label_size = 16
        
    def plot_sequence_freq(self, full0, full1, freq0, freq1, ref_string=None, normalization=1, labels_on=True):
        """ 
        Plot a bar graph of all sequences and their ratio of frequencies before/after screening 
        
        Parameters:
        ----------
        full0: dict
            Full set of sequencing data before screening, sequence-count (key, value) pair
        full1: dict 
            Full set of sequencing data after screening, sequence-count (key, value) pair
        freq0: dict
            Sampled set of sequencing data before screening, sequence-count (key, value) pair
        freq1: dict
            Sampled set of sequencing data after screening, sequence-count (key, value) pair
        ref_string: str
            Reference string to normalize the data, must have length(sequence)
        normalization: anything
            Determines whether to do a normalization
        labels_on: bool (default True)
            Labels bars with the corresponding sequence
            
        Returns:
        -------
        seq_freq_plot: matplotlib bar plot
        """

        seq_labels = freq0.keys()
        x_marks = np.arange(len(freq0))
        
        if normalization is not 1:
            if ref_string is None:
                # Dirty hack to get string of ones with len(sequences)
                ref_string = "".join([str(int(one)) for one in np.ones(len(seq_labels[0]))])
            elif len(ref_string) is not len(seq_labels[0]):
                raise NameError('ref_string must have length: ' + str(len(seq_labels[0])))
            normalization = float(freq1[ref_string]*full0[ref_string])/(freq0[ref_string]*full1[ref_string])
            
        ratio = list()
        for s in seq_labels:
            ratio.append(float(freq1[s])/(freq0[s])*(1/normalization))
        
        seq_freq_plot = plt.bar(x_marks, ratio, align='center', alpha=.7, color='red')
        plt.tick_params(axis='x', which='both', bottom='off', labelbottom='off')
        
        if labels_on:
            for x in x_marks:
                plt.text(x_marks[x], ratio[x], seq_labels[x])
                
        plt.title('Sequence Frequencings after Screening', fontsize=self.title_size)
        plt.ylabel('Frequency ratio ($f_1/f_0$)', fontsize=self.y_label_size)
        
        return seq_freq_plot
        
    def plot_kmer_interactions(self, kmer_interactions):
        """ Plot a bar chart of all k-order interactions for the system """
        
        seq_labels = kmer_interactions.keys()
        seq_values = kmer_interactions.values()
        x_marks = np.arange(len(seq_labels))
        
        kmer_plot = plt.bar(x_marks, seq_values, align='center', alpha=0.7, color='blue')
        plt.xticks(x_marks, seq_labels, rotation=50)
        
        plt.title('Frequency of sequences with '+ str(len(seq_labels[0])) + '-mer present', fontsize=self.title_size)
        plt.ylabel('Frequency ratio ($f_1/f_0$)', fontsize=self.y_label_size)
        
        return kmer_plot