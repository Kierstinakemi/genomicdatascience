# course: CIS 5900
# Professor: Dr. Giri
# date: 07/27/2018
# name: Kierstin Matsuda
# description: helper functions for genome.py

import random

##############################################################
# sequencerSimulator(genome, numReads, readLength)
# This function will act like a sequencer by creating reads, 
# except not from slides it will be from an input genome
# inputs: genome, how many reads, length of read desired)
# returns: the reads
def sequencerSimulator(genome, numReads, readLength):
    reads = []
    for _ in range(numReads):
        start = random.randint(0, len(genome)-readLength) -1
        reads.append(genome[start : start + readLength])
    return reads
