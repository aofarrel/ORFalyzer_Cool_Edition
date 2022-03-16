# finalproject.py
# Name: Ash O'Farrell (aofarrel)
# Group members: Blaise Nasri (banasri)

'''
Program Overview: This program takes a prokaryotic genome 
in the form of a FASTA file as and gives a breakdown of Open 
Reading frames and an analysis of the proteins that are coded 
for within them. This analysis includes: amino acid composition
and break down, pI, pH, molar extinction coefficient, mass 
extinction coefficient, and molecular weight. There are also 
user friendly options incorporated to accommodate 
preferences of the individual running the program. Relevant
documentation for these preferences can be seen below.

Input: FASTA file
Output: text file with ORFs, size of ORFs, amino acid 
composition and break down, pI, pH, molar extinction 
coefficient, mass extinction coefficient, & molecular weight.

Abbreviations: 
PP   = ProteinParams
NP   = NucParams
ORF  = ORF Finder

'''


import sequenceanalysis
import argparse


class CommandLine():
    '''
    Handle the command line, usage and help requests.
    CommandLine uses argparse, now standard in 2.7 and beyond.
    it implements a standard command line argument parser with various argument options,
    a standard usage and help.
    attributes:
    all arguments received from the commandline using .add_argument will be
    avalable within the .args attribute of object instantiated from CommandLine.
    For example, if myCommandLine is an object of the class, and requiredbool was
    set as an option using add_argument, then myCommandLine.args.requiredbool will
    name that option.
'''
    def __init__(self, inOpts=None) :
        '''Implement a parser to read command line argv string using argparse'''
        self.parser = argparse.ArgumentParser(
            description='orfalyzer.py -- Find ORFs in prokaryotic DNA/RNA sequence and output protein properties',
            epilog='',
            add_help=True, prefix_chars='-',
            usage='%(prog)s [options] input.fna output.txt')
        self.parser.add_argument('inFile',
                                 action='store',
                                 help='input FASTA file name')
        self.parser.add_argument('outFile',
                                 action='store',
                                 help='output file name')

        self.parser.add_argument('-mG', '--minGene',
                                 type=int,
                                 choices=(-1, 100, 200, 300, 500, 1000),
                                 default=100,
                                 action='store',
                                 help='minimum gene length, -1 = no minimum')

        self.parser.add_argument('-s', '--start',
                                 action='append',
                                 default=['ATG'],
                                 nargs='?',
                                 help='start codon(s)')  # allows multiple list options

        self.parser.add_argument('-t', '--stop',
                                 action='append',
                                 default=['TAG', 'TGA', 'TAA'],
                                 nargs='?',
                                 help='stop codon(s)')  # allows multiple list options

        self.parser.add_argument('-c', '--cytosine',
                                 action='store',
                                 default=True,
                                 help='If True (default), cytosine under oxidizing conditions. If false, reducing conditions.')

        self.parser.add_argument('-v', '--version',
                                 action='version',
                                 version='Orf-A-Lyzer 1.0')

        if inOpts is None :
            self.args = self.parser.parse_args()
        else :
            self.args = self.parser.parse_args(inOpts)

def main(inCL=None):
    '''Create OrfFinder object to find some genes'''
    if inCL is None:
        myCommandLine = CommandLine()
    else :
        myCommandLine = CommandLine(inCL)
    thisOrfFinder = sequenceanalysis.OrfFinder()
    thisOrfFinder._processArgs_(myCommandLine.args.minGene,
                                True,
                                myCommandLine.args.outFile,
                                myCommandLine.args.start,
                                myCommandLine.args.stop,
                                myCommandLine.args.cytosine)
    readInFile = sequenceanalysis.FastaReader(myCommandLine.args.inFile)
    for head, seq in readInFile.readFasta():
        thisOrfFinder._processHeader_(head, readInFile.headerNumber)
        thisOrfFinder._addsequence_(seq)
    # the remainder is performed in sequenceanalysis.py

if __name__ == "__main__":
    main()