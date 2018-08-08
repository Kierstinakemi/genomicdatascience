# course: CIS 5900
# Professor: Dr. Giri
# date: 07/27/2018
# name: Kierstin Matsuda
# description: sequencing functions and helper functions

import random
import collections
import matplotlib.pyplot as plt

"""defining sequences
dnaSeq1 =''
dnaSeq2 =''

basepairs = {'A':'T', 'T':'A', 'C':'G', 'G':'C'}

#filling sequences with random nucleotides
for _ in range(50):
   dnaSeq1 += random.choice('ACTG')
   dnaSeq2 += random.choice('ACTG')

print("DNA Seq 1:" + dnaSeq1)
print("DNA Seq 2:" + dnaSeq2)"""


##############################################################
# longestCommonPrefix(dnaSeq1, dnaSeq2)
# This function calculates the Longest Common Prefix from two reads
# inputs: Two sequences of A's T's C's and G's
# returns: the prefix which both sequences have, from the first sequence
def longestCommonPrefix(dnaSeq1, dnaSeq2):
    nucleotide=0
    while nucleotide < len(dnaSeq1) and nucleotide < len(dnaSeq2) and dnaSeq1[nucleotide] == dnaSeq2[nucleotide]:
        nucleotide += 1
# return the section up to the last nucleotide not including it
    return dnaSeq1[0:nucleotide] 

"""print longestCommonPrefix(dnaSeq1, dnaSeq2)"""

##############################################################
# isMatch(dnaSeq1, dnaSeq2)
# This function determine if two reads match for every nucleotide
# inputs: Two sequences of A's T's C's and G's
# returns: True if match, false if no match
def isMatch(dnaSeq1, dnaSeq2):
    if not len(dnaSeq1) == len(dnaSeq2):
        return False

    for nucleotide in range(len(dnaSeq1)):
        if not dnaSeq1[nucleotide] == dnaSeq2[nucleotide]:
            return False
    return True


##############################################################
# reverseNucleotideComplement(dnaSeq)
# This function completes the second half of the string like a DNA polymerase
# inputs: A sequences of A's T's C's and G's
# returns: the opposite side of the dna string as sequence of A's T's C's and G's
def reverseNucleotideComplement(dnaSeq):
    nucleotideComplement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    complement = ''
    for nucleotide in dnaSeq:
        complement = complement + nucleotideComplement[nucleotide] # + complement
    return complement

"""print 'AACTTGTAGCCTA'
print reverseNucleotideComplement('AACTTGTAGCCTA')"""



##############################################################
# readGenome(file)
# This function will open file and read genome
# inputs: the file name in the same directory
# returns: the reads from within the file
def readGenome(file):
    genome = ''
    #instead of try finally
    with open(file, 'r') as genomeFile:
        for line in genomeFile:
            if not line[0] =='>':
                genome += line.rstrip()
    return genome

#print readGenome('lambda_virus.fa')


##############################################################
# countNumBases(genome)
# This function counts the number of bases
# (nucleotides) of input string of dna or genome
# inputs: sequence of A's T's C's and G's
# returns: a dictionary of the counts
def countNumBases(genome):
    numBases = { 'A': 0, 'T': 0, 'C': 0, 'G': 0, 'N': 0} # N for no confidence/null
    for base in genome:
        numBases[base] += 1
    return numBases

"""print countNumBases(readGenome('lambda_virus.fa'))"""



##############################################################
# calcHypotenuse(a, b)
# This function solves Pythagorean theorem a^2 + b^2 = c^2
# for the value of c.
# inputs: a and b are the lengths of sides of a right triangle.
# returns: the length of the hypotenuse.
# changing ascii characters into their corresponding base quality index
def qToPhred(qValue):
    return chr(qValue + 33) #adds 33 then changes it into corresponding ASCII integer


##############################################################
# calcHypotenuse(a, b)
# This function solves Pythagorean theorem a^2 + b^2 = c^2
# for the value of c.
# inputs: a and b are the lengths of sides of a right triangle.
# returns: the length of the hypotenuse.
# changing ascii characters into their corresponding base quality index
def phredToQ(phredValue):
    return ord(phredValue) - 33 #changes int to ASCII character then subtracts 33


##############################################################
# readFastqFile(file)
# This function reads a file of fastq file - a file format for
# a dna sequence and its quality reads
# inputs: the name of the .fastq or .fa file in the same directoy
# returns: the separated sequences of reads and their quality scores
def readFastqFile(file):
    sequences = []
    qualities = []
    
    with open(file) as fastq:
        while True:
            fastq.readline()
            read = fastq.readline().rstrip()
            fastq.readline()
            quality = fastq.readline().rstrip()
            if len(read) == 0:
                break
            sequences.append(read)
            qualities.append(quality)
        return sequences, qualities

"""reads, readQualities = readFastqFile('ERR037900_1.first1000.fastq')
print reads[0:15]
print readQualities[0:15]"""



##############################################################
# createHistogram(qualities)
# This function creates a histogram of the read qualities
# inputs: the sequences of qualities of reads from a .fastq file
# returns: the histogram data
def createHistogram(qualities):
    histogram = [0] * 50
    for quality in qualities:
        for phred in quality:
            q = phredToQ(phred)
            histogram[q] +=1
    return histogram

"""histogram = createHistogram(readQualities)
print histogram

print phredToQ('H')
print phredToQ('J')
print phredToQ('#')

plt.bar(range(len(histogram)), histogram)
plt.show()"""


##############################################################
# findGCPosition(reads)
# This function finds where the G's and C's are by position
# inputs: stripped dna reads without any metadata
# returns: the positions of the g's and c's in an array
def findGCPosition(reads):
    gcPos = [0] * 100
    totals = [0] * 100

    for read in reads:
        for base in range(len(read)):
            if read[base] == 'C' or read[base] == 'G':
                gcPos[base] += 1
            totals[base] += 1

        for base in range(len(gcPos)):
            if totals[base] > 0:
                gcPos[base] /= float(totals[base])
    return gcPos


"""GC = findGCPosition(reads)
plt.plot(range(len(GC)), GC)
plt.show()

print countNumBases(reads)"""

##############################################################
# naive(patternToMatch, exampleGenome)
# This function will try and match reads to a genome
# inputs: the sequence or read to match, the genome to try and match
# returns: the indexes where there are matches to the example genome
def naive(patternToMatch, exampleGenome):
    occurrences = [] #to hold indexes of where there are matches
    for i in range(len(exampleGenome) - len(patternToMatch) + 1):  # loop over alignments
        match = True
        for j in range(len(patternToMatch)):  # loop over characters
            if exampleGenome[i+j] != patternToMatch[j]:  # compare characters
                match = False
                break
        if match:
            occurrences.append(i)  # all chars matched; record
    return occurrences