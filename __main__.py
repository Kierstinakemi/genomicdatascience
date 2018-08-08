# course: CIS 5900
# Professor: Dr. Giri
# date: 07/27/2018
# name: Kierstin Matsuda
# description: testing various sequencing functions on 3 genomes

#Imports
#longestCommonPrefix, isMatch, reverseNucleotideComplement, readGenome, countNumBases, qToPhred, phredToQ, readFastqFile, createHistogram, findGCPosition
from genome import * 
from genomeHelper import sequencerSimulator
import random
import collections
import matplotlib.pyplot as plt



"""The main function."""
def main(args=None):
   
    #if args is None:
    #    args = sys.argv[1:]

    print("\nHere are tests of some genome sequencing functions.")
    print("The purpose of this program is to test the functions on 3 different genomes.\n")


    print "\nLAMBDA VIRUS GENOME"
    lambdaVirus = readGenome('lambda_virus.fa')
    print lambdaVirus[0:5]

    #defining dna sequences
    dnaSeq1 =''
    dnaSeq2 =''
    #defining nucleotide base pairs
    basepairs = {'A':'T', 'T':'A', 'C':'G', 'G':'C'}

    #filling sequences with random nucleotides
    for _ in range(50):
        dnaSeq1 += random.choice('ACTG')
        dnaSeq2 += random.choice('ACTG')

    print '\nThe DNA Sequences to test longest Common Prefix:'
    print("DNA Seq 1:" + dnaSeq1)
    print("DNA Seq 2:" + dnaSeq2)

    print '\nThe longest common prefix of DNA Sequence 1 and DNA Sequence 2 is:'
    print longestCommonPrefix(dnaSeq1, dnaSeq2) +'\n'

    print '\nCreating the second strand of the double helix using one strand'
    print 'The strand to complete: AACTTGTAGCCTA'
    print 'The opposite side:      '+reverseNucleotideComplement('AACTTGTAGCCTA')

    print 'The number of A''s T''s G''S and C''s in the lambda virus genome'
    print countNumBases(lambdaVirus)


    reads, readQualities = readFastqFile('ERR037900_1.first1000.fastq')
    print '\nSome sequencing machine reads of a human genome'
    print reads[0:5]
    print '\nThe read qualities of each of the reads'
    print readQualities[0:5]

    print '\nA histogram of the read qualities of the previous human genome should pop up'
    histogram = createHistogram(readQualities)
    print histogram

    """
    print phredToQ('H')
    print phredToQ('J')
    print phredToQ('#')
    """

    print '\nA bar graph of the previous histogram data'
    plt.bar(range(len(histogram)), histogram)
    plt.show()

    print '\nFinding the ---'
    GC = findGCPosition(reads)
    plt.plot(range(len(GC)), GC)
    plt.show()

    #print countNumBases(reads)

    
    phixGenome = readGenome('phix.fa')

    templateSequence = 'ACGTGCCCAGTCAGTCGCAA'
    read = 'GC'

    print naive(read, templateSequence)

    #Creating a mock fresh DNA read and match
    print '\nCreating mock reads from an input genome and printing them:'
    mockReads = sequencerSimulator(lambdaVirus, 10, 100)
    print mockReads

    #reading mock reads to see if it matches it with its original genome
    matches = 0
    for read in mockReads:
        match = naive(read, lambdaVirus)
        if len(match) > 0:
            matches += 1
    print('\n%d / %d of the reads matched the genome exactly.' % (matches, len(mockReads)))

    """
    #Trying he same thing with real reads
    phixReads, _ = readFastqFile('phix.fa')

    matches = 0
    for read in phixReads:
        match = naive(read, )
        match.extend(naive(reverseNucleotideComplement(), phixReads))
        if len(match) > 0:
            matches += 1
    print('\n%d / %d of the reads matched the genome exactly.' % (matches, len(phixReads)))
    """
    
    



if __name__ == "__main__":
    main()