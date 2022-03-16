#!/usr/bin/env python3
# Name: Ash O'Farrell (aofarrel)
# Group Members: Blaise Nasri (banasri)
'''sequenceanalysis.py -- Orf-A-Lyzer Edition

This program consists of four classes:
* FastaReader, used for parasing fasta files
* NucParams, used to analyze DNA/RNA data
* ProteinParams, used to analyze amino acid chains
* OrfFinder, used to analyze open reading frames

OrfFinder finds open reading frames (ORFs) in a DNA sequence that is
being parsed by FastaReader. It scans all three possible reading frames
per strand, and will search both the plus and minus strands.

In this version of sequenceanalysis, the contents of the entire
ORF is saved. In order to prevent ORFs from being generated and
saved to the disk unnecessarily, slowing computational time, 
the reverse strand method was re-written compared to what
Ash turned in for Lab 5. It now essentially reverse-compliments
the input string and runs on that instead of trying to get clever
with just reverse-complimenting the start and stop codons.

The output file will NOT be cleared before the results are written.
This was an intentional design decision to allow for the possibility
that multiple header-sequence pairs may exist in the inputfile.

Input: sequence data (DNA or RNA) as string
Output: ORF location, ORF content, predicted protein parameters'''

import sys

class FastaReader:
    ''' 
    Define objects to read FastA files.

    instantiation: 
    thisReader = FastaReader ('testTiny.fa')
    usage:
    for head, seq in thisReader.readFasta():
        print (head,seq)
    '''

    def __init__(self, fname=''):
        '''contructor: saves attribute fname '''
        self.fname = fname
        self.headerNumber = 0

    def doOpen(self):
        ''' Handle file opens, allowing STDIN.'''
        if self.fname == '':
            return sys.stdin
        else:
            return open(self.fname)

    def readFasta(self):
        ''' Read an entire FastA record and return the sequence header/sequence'''
        header = ''
        sequence = ''

        with self.doOpen() as fileH:

            header = ''
            sequence = ''

            # skip to first fasta header
            line = fileH.readline()
            while not line.startswith('>'):
                line = fileH.readline()
            header = line[1:].rstrip()

            for line in fileH:
                if line.startswith('>'):
                    self.headerNumber += 1
                    yield header, sequence
                    header = line[1:].rstrip()
                    sequence = ''
                else:
                    sequence += ''.join(line.rstrip().split()).upper()

        yield header, sequence

class OrfFinder:
    def __init__(self):
        '''Intialize the class.'''
        self.bigResultsList = []  # required to sort by len(ORF)

    def _addsequence_(self, inputSequence):
        '''Process "seq" passed from FastaReader'''
        self.scan(inputSequence)

    def _processHeader_ (self, header, headerNumber):
        '''Processes "header" (>) passed from FastaReader'''
        self.header = header
        self.headerNumberFromFR = headerNumber  # passed in from modified FASTA reader

    def _processArgs_(self, minOrf, biggestOrfOnly, outputfile, startcodonlist, stopcodonlist, cytosine):
        '''Ideally, arguments should be passed via argparse.'''
        self.minOrf = minOrf
        self.biggestOrfOnly = biggestOrfOnly
        self.outputfile = outputfile
        self.startcodons = startcodonlist
        self.stopcodons = stopcodonlist
        self.cytosine = cytosine

    def scan(self, inputSequence):
        '''Scan the input sequence for orfs'''

        with open(self.outputfile, 'a') as file:
            # a instead of w to allow for a file with multiple header/seq combos
            file.write("> %s" % self.header)  # add the > so fasta reader can skip the headers in outfile
            file.write("\n")

        lenInputSeq = len(inputSequence)
        self.bigResultsList = []
        self.debugThisResult = False

        # Forward direction
        for frame in [0, 1, 2]:

            self.noMoreHangingORFs = False
            '''
            Turns to True the first time we get a start codon, in order
            to prevent something like this from being reported as two ORFs:
               ATG NNN NNN TAA NNN TAA
               srt         stp     stp
            This should be reset every frame.
            '''

            self.tempResults = []  # clear every iteration and stop codon
            self.everyORFTempResults = []
            ORFasString = ""
            thisStartCodon = "None"  # in case we run into a stop codon before a start
            startPosition = 1
            revcom = False
            for ix in range(frame, lenInputSeq, 3):
                codon = inputSequence[ix:ix + 3]
                if codon in self.startcodons:
                    if thisStartCodon == "None":  # first start codon
                        if self.noMoreHangingORFs:
                            ORFasString = "" # clear ORFasString because now we start from here
                            ORFasString += codon
                            startPosition = ix + 1
                            thisStartCodon = codon

                        if self.biggestOrfOnly is not True:
                            # store every ORF we find until we find a stop codon
                            # end (1), ix (3) and endPosition (5) added once stop found
                            startPosition = ix + 1
                            self.everyORFTempResults.append([thisStartCodon, int(frame), startPosition, revcom, lenInputSeq, ORFasString])
                        ix += 3
                    else:
                        '''We already ran into a start codon earlier;
                        this is a little ORF within an ORF'''
                        ORFasString += codon
                        pass
                elif codon in self.stopcodons:
                    '''Stop codon found. This may be either a complete or a hanging ORF, but
                    in any case, this is as big as the ORF can get.'''

                    if thisStartCodon == "None":  # possible hanging ORF
                        if self.noMoreHangingORFs:
                            pass  # just a stop codon after a complete ORF
                        else: # print hanging ORF
                            ORFasString += codon
                            ix += 3  # do before storing ix or else it'll be off by 3
                            end = codon
                            endPosition = ix
                            results = [thisStartCodon, end, int(frame), ix,
                                                startPosition, endPosition, revcom, lenInputSeq, ORFasString]
                            self.writeResults(results)
                        self.noMoreHangingORFs = True

                    else:  # complete ORF
                        ORFasString += codon
                        ix += 3  # do before storing ix or else it'll be off by 3
                        end = codon
                        endPosition = ix
                        if self.biggestOrfOnly is not True:
                            for listy in self.tempResults:
                                listy.insert(1, end)
                                listy.insert(3, ix)
                                listy.insert(5, endPosition)
                                self.writeResults(listy)
                            self.tempResults = []
                        else:  # needs to be else to prevent printing twice
                            results = [thisStartCodon, end, int(frame), ix,
                                       startPosition, endPosition, revcom, lenInputSeq, ORFasString]
                            self.writeResults(results)
                    thisStartCodon = "None"  # needed for cases 2 and 3

                else: # neither stop nor start codon
                    if self.noMoreHangingORFs: # all ORFs must have a start codon now
                        if thisStartCodon != "None": # we have a start codon, keep saving until stop codon
                            ORFasString += codon
                    else: # might be a hanging ORF, save for now
                        ORFasString += codon
        revcomInput = self._revcom_(inputSequence)
        #reverse direction, very similiar to forward direction
        for frame in [0, 1, 2]:

            self.noMoreHangingORFs = False
            '''
            Turns to True the first time we get a start codon, in order
            to prevent something like this from being reported as two ORFs:
               ATG NNN NNN TAA NNN TAA
               srt         stp     stp
            This should be reset every frame.
            '''

            self.tempResults = []  # clear every iteration and stop codon
            self.everyORFTempResults = []
            ORFasString = ""
            thisStartCodon = "None"  # in case we run into a stop codon before a start
            startPosition = 1
            revcom = True
            for ix in range(frame, lenInputSeq, 3):
                codon = revcomInput[ix:ix + 3]
                if codon in self.startcodons:
                    if thisStartCodon == "None":  # first start codon
                        if self.noMoreHangingORFs:
                            ORFasString = "" # clear ORFasString because now we start from here
                            ORFasString += codon
                            startPosition = ix + 1
                            thisStartCodon = codon

                        if self.biggestOrfOnly is not True:
                            # store every ORF we find until we find a stop codon
                            # end (1), ix (3) and endPosition (5) added once stop found
                            startPosition = ix + 1
                            self.everyORFTempResults.append([thisStartCodon, int(frame), startPosition, revcom, lenInputSeq, ORFasString])
                        ix += 3
                    else:
                        '''We already ran into a start codon earlier;
                        this is a little ORF within an ORF'''
                        ORFasString += codon
                        pass
                elif codon in self.stopcodons:
                    '''Stop codon found. This may be either a complete or a hanging ORF, but
                    in any case, this is as big as the ORF can get.'''

                    if thisStartCodon == "None":  # possible hanging ORF
                        if self.noMoreHangingORFs:
                            pass  # just a stop codon after a complete ORF
                        else: # print hanging ORF
                            ORFasString += codon
                            ix += 3  # do before storing ix or else it'll be off by 3
                            end = codon
                            endPosition = ix
                            results = [thisStartCodon, end, int(frame), ix,
                                                startPosition, endPosition, revcom, lenInputSeq, ORFasString]
                            self.writeResults(results)
                        self.noMoreHangingORFs = True

                    else:  # complete ORF
                        ORFasString += codon
                        ix += 3  # do before storing ix or else it'll be off by 3
                        end = codon
                        endPosition = ix
                        if self.biggestOrfOnly is not True:
                            for listy in self.tempResults:
                                listy.insert(1, end)
                                listy.insert(3, ix)
                                listy.insert(5, endPosition)
                                self.writeResults(listy)
                            self.tempResults = []
                        else:  # needs to be else to prevent printing twice
                            results = [thisStartCodon, end, int(frame), ix,
                                       startPosition, endPosition, revcom, lenInputSeq, ORFasString]
                            self.writeResults(results)
                    thisStartCodon = "None"  # needed for cases 2 and 3

                else: # neither stop nor start codon
                    if self.noMoreHangingORFs: # all ORFs must have a start codon now
                        if thisStartCodon != "None": # we have a start codon, keep saving until stop codon
                            ORFasString += codon
                    else: # might be a hanging ORF, save for now
                        ORFasString += codon

        self.bigResultsList.sort(key=lambda x: (x[3], x[1]), reverse=True)
        with open(self.outputfile, 'a'):
            for listy in self.bigResultsList:
                if listy[4]:  # minus strand
                    with open(self.outputfile, 'a') as out:
                        out.write("-{:d} {:>5d}...{:>5d} {:>5d}\n".format(listy[0], listy[1], listy[2], listy[3]))
                        NP = NucParams(True)
                        NP.addSequence(listy[5])
                        PP = ProteinParam(NP.aaAsString(), self.cytosine)
                        out.write("┏━━━━━━━━━━━━━━━━┓\n")
                        out.write("┃  ORF in frame  ┃\n")
                        out.write("┡━━━━━━━━━━━━━━━━┩\n")
                        out.write(listy[5])
                        out.write("\n")
                        out.write("┏━━━━━━━━━━━━━━━━┓\n")
                        out.write("┃  Resulting AA  ┃\n")
                        out.write("┡━━━━━━━━━━━━━━━━┩\n")
                        if (NP.aaAsString() != '-'):
                            out.write("{}\n".format(NP.aaAsString()))
                        else:
                            out.write("None -- just a stop codon\n")
                        out.write("┏━━━━━━━━━━━━━━━━━━━━━━━┓\n")
                        out.write("┃  Protein properties   ┃\n")
                        out.write("┡━━━━━━━━━━━━━━━━━━━━━━━┩\n")
                        if (NP.aaAsString() != '-'):
                            out.write("pI value: %s\n" % PP.pI())
                            out.write("Molar Extinction: %s\n" % PP.molarExtinction())
                            out.write("Mass Extinction: %s\n" % PP.massExtinction())
                            out.write("Molecular Weight: %s\n" % PP.molecularWeight())
                        else:
                            out.write("None -- just a stop codon\n")
                        myAAcomposition = PP.aaComposition()
                        keys = list(myAAcomposition.keys())
                        keys.sort()
                        # handles the case where no AA are present
                        if PP.aaCount() == 0:
                            pass
                        elif PP.aaCount() == 1:
                            pass
                        else:
                            out.write("Amino acid composition:\n")
                            for key in keys:
                                out.write("\t{} = {:.2%}\n".format(key, myAAcomposition[key] / PP.aaCount()),
                                    )
                        out.write("\n^v^v^v^v^v^v^v^v^v^v^v^v^v^v^v^v^v^v^v^v^v^v^v^v^v^v^v^v^v^v^\n\n")
                else:  # plus strand
                    NP = NucParams(True)
                    NP.addSequence(listy[5])
                    PP = ProteinParam(NP.aaAsString(), self.cytosine)
                    with open(self.outputfile, 'a') as out:
                        out.write("{:+d} {:>5d}...{:>5d} {:>5d}\n".format(listy[0], listy[1], listy[2], listy[3]))
                        out.write("┏━━━━━━━━━━━━━━━━┓\n")
                        out.write("┃  ORF in frame  ┃\n")
                        out.write("┡━━━━━━━━━━━━━━━━┩\n")
                        out.write(listy[5])
                        out.write("\n")
                        out.write("┏━━━━━━━━━━━━━━━━┓\n")
                        out.write("┃  Resulting AA  ┃\n")
                        out.write("┡━━━━━━━━━━━━━━━━┩\n")
                        if (NP.aaAsString() != '-'):
                            out.write("{}\n".format(NP.aaAsString()))
                        else:
                            out.write("None -- just a stop codon\n")
                        out.write("┏━━━━━━━━━━━━━━━━━━━━━━━┓\n")
                        out.write("┃  Protein properties   ┃\n")
                        out.write("┡━━━━━━━━━━━━━━━━━━━━━━━┩\n")
                        if (NP.aaAsString() != '-'):
                            out.write("pI value: %s\n" % PP.pI())
                            out.write("Molar Extinction: %s\n" % PP.molarExtinction())
                            out.write("Mass Extinction: %s\n" % PP.massExtinction())
                            out.write("Molecular Weight: %s\n" % PP.molecularWeight())
                        else:
                            out.write("None -- just a stop codon\n")
                        myAAcomposition = PP.aaComposition()
                        keys = list(myAAcomposition.keys())
                        keys.sort()
                        # handles the case where no AA are present
                        if PP.aaCount() == 0:
                            pass
                        elif PP.aaCount() == 1:
                            pass
                        else:
                            out.write("Amino acid composition:\n")
                            for key in keys:
                                out.write("\t{} = {:.2%}\n".format(key, myAAcomposition[key]/PP.aaCount()))
                        out.write("\n^v^v^v^v^v^v^v^v^v^v^v^v^v^v^v^v^v^v^v^v^v^v^v^v^v^v^v^v^v^v^\n\n")
                    

    def writeResults(self, results):
        '''Process the results list passed in from scan(). Technically, the
        actually writing is done in scan().'''
        if results == []:
            pass
        else:

            resultsFrame = results[2]+1
            resultsStartPos = results[4]
            resultsStopPos = results[5]
            resultsStrand = results[6]
            resultslenInputSeq = results[7]

            if resultsStrand == False: # plus strand
                resultsORFasString = results[8]  # plus strand only for now
                resultsLen = results[5] - results[4] + 1
                if resultsLen >= self.minOrf:
                    self.bigResultsList.append([resultsFrame, resultsStartPos, resultsStopPos, resultsLen, False, resultsORFasString])
                else:
                    pass
            else:  # minus strand
                resultsORFasString = results[8]
                resultsLen = results[5] - results[4] + 1
                if resultsLen >= self.minOrf:
                    self.bigResultsList.append([resultsFrame, resultsStartPos, resultsStopPos, resultsLen, True, resultsORFasString])
                else:
                    pass

    def _revcom_(self, inputSequence):
        '''Generate reverse compliment of inputSequence'''
        revcomInput = inputSequence[::-1].translate(str.maketrans('TACG', 'ATGC'))
        return revcomInput


class NucParams:
    '''Analyzes a DNA or RNA string for aa compositon, nt composition,
    codon composition, and total number of nts. Requires sequence 
    data to either be passed in when constructing the object or by 
    using the addSequence() method.'''
    rnaCodonTable = {
        # RNA codon table
        # U
        'UUU': 'F', 'UCU': 'S', 'UAU': 'Y', 'UGU': 'C',  # UxU
        'UUC': 'F', 'UCC': 'S', 'UAC': 'Y', 'UGC': 'C',  # UxC
        'UUA': 'L', 'UCA': 'S', 'UAA': '-', 'UGA': '-',  # UxA
        'UUG': 'L', 'UCG': 'S', 'UAG': '-', 'UGG': 'W',  # UxG
        # C
        'CUU': 'L', 'CCU': 'P', 'CAU': 'H', 'CGU': 'R',  # CxU
        'CUC': 'L', 'CCC': 'P', 'CAC': 'H', 'CGC': 'R',  # CxC
        'CUA': 'L', 'CCA': 'P', 'CAA': 'Q', 'CGA': 'R',  # CxA
        'CUG': 'L', 'CCG': 'P', 'CAG': 'Q', 'CGG': 'R',  # CxG
        # A
        'AUU': 'I', 'ACU': 'T', 'AAU': 'N', 'AGU': 'S',  # AxU
        'AUC': 'I', 'ACC': 'T', 'AAC': 'N', 'AGC': 'S',  # AxC
        'AUA': 'I', 'ACA': 'T', 'AAA': 'K', 'AGA': 'R',  # AxA
        'AUG': 'M', 'ACG': 'T', 'AAG': 'K', 'AGG': 'R',  # AxG
        # G
        'GUU': 'V', 'GCU': 'A', 'GAU': 'D', 'GGU': 'G',  # GxU
        'GUC': 'V', 'GCC': 'A', 'GAC': 'D', 'GGC': 'G',  # GxC
        'GUA': 'V', 'GCA': 'A', 'GAA': 'E', 'GGA': 'G',  # GxA
        'GUG': 'V', 'GCG': 'A', 'GAG': 'E', 'GGG': 'G'  # GxG
    }
    dnaCodonTable = {key.replace('U', 'T'): value for key, value in rnaCodonTable.items()}

    def __init__(self, boolSaveString:bool, inString=''):
        self.saveAAString = boolSaveString  # less efficient
        self.cleanedNTstring = "" # modified if self.saveAAString is True
        self.nucComp = {'A':0, 'C':0, 'G':0, 'T':0, 'U':0, 'N':0} # modified by addSequence()
        self.nucNum = 0 # modified by nucCount()
        self.codonComp = dict.fromkeys(self.rnaCodonTable, 0) # modified by addSequence()
        aminoacids = {value for value in self.rnaCodonTable.values()}
        self.aaComp = dict.fromkeys(aminoacids, 0) # modified by aaComposition()


    def addSequence(self, inSeq):
        '''Adds sequence, generating nucComp and codonComp in the process'''
        # calculate nucleotide composition, ignoring invalid nts
        for nucleotide in inSeq:
            if nucleotide.isalpha():
                nucleotide = nucleotide.upper()
                if nucleotide in self.nucComp:
                    self.nucComp[nucleotide] += 1

        # calculate codon composition
        if self.nucCount() % 3 != 0:
            print("Warning: Input is not a multiple of three. ", end="")
            print("The resulting incomplete codon at the end will be ignored.")
        else:
            pass
        n = 0
        while n < len(inSeq) - 2:
            nucleotideIndex = -1
            currentCodon = [inSeq[n], inSeq[n+1], inSeq[n+2]]
            nFlag = False
            for nucleotide in currentCodon:
                nucleotideIndex += 1
                if nucleotide == 'N':
                    nFlag = True  # break statment doesn't work
                else:
                    if nucleotide == 'T':
                        currentCodon[nucleotideIndex] = 'U'
                    else:
                        pass
            if nFlag:
                pass  # codon doesn't get added
            else:
                currentCodonString = "".join(currentCodon)
                if self.saveAAString == True:
                    # We don't just pass in inseq because we still want to check for N
                    self.cleanedNTstring += ''.join(currentCodon)
                self.codonComp[currentCodonString] += 1
            n += 3

    def aaComposition(self):
        '''Outputs dictionary of AA composition based on codon composition'''
        for codon in self.codonComp:
            timesToLoop = self.codonComp[codon]
            while timesToLoop > 0:
                # loop once per every time that codon shows up
                self.aaComp[self._rna2aa_(codon)] += 1
                timesToLoop -= 1
        return self.aaComp

    def aaAsString(self):
        '''Outputs AAs as a string. Only use if self.saveAAString = True'''
        self.fullAAString = ""
        if self.saveAAString == False:
            return ""
        else:
            n = 0

            while n < len(self.cleanedNTstring) - 2:
                currentCodon = [self.cleanedNTstring[n], self.cleanedNTstring[n + 1], self.cleanedNTstring[n + 2]]
                self.fullAAString += self._rna2aa_(''.join(currentCodon))
                n += 3
            return self.fullAAString


    def _rna2aa_(self, rnaString):
        '''Takes in 3 letter RNA codon, outputs - (if stop codon) or resulting AA'''
        if rnaString in self.rnaCodonTable:
            return self.rnaCodonTable[rnaString]
        elif rnaString in self.dnaCodonTable:
            return self.dnaCodonTable[rnaString]
        else:
            print("Warning: Invalid RNA codon")

    def nucComposition(self):
        '''Returns dictionary of nucleotide composition'''
        return self.nucComp

    def codonComposition(self):
        '''Returns dictionary of codon compositon'''
        return self.codonComp

    def nucCount(self):
        '''Returns number of valid nucleotides'''
        if sum(self.nucComp.values()) == 0:
                # No valid nucleotides found
                self.nucNum = 0
                return self.nucNum
        else:
            self.nucNum = sum(self.nucComp.values())
        return self.nucNum

class ProteinParam:
    '''
    Calculates various properties of a protein given its AA structure:
    * Number of (valid) amino acids
    * Molecular weight
    * Molar extinction coefficient
    * Mass extinction coefficient
    * Theoretical pI
    * Amino acid composition

    By setting self.cytosine = False, the user can get molar and mass
    extinction coefficient under reducing conditions. This program will
    assume oxidizing conditions by default. Note that as we were asked
    to not change main(), the user can't edit the value of self.cytosine 
    at runtime; the user must edit the value in the code itself.

    Binary search is used to determine theoretical pI. It is capped to a
    hundred iterations and uses rounding in order to prevent excessively
    long computation time and/or infinite iterations.

    The user could smash their face across the keyboard or enter hiragana
    characters if they wanted, but only valid single-letter AA codes will
    be evaulated. Lowercase input will be handled as if it were uppercase.

    This one line in __init__()...

    self.aaDictIncludeZeros = dict.fromkeys(self.aa2mw,0)

    ...was based on Daniel Roseman's answer on SO.
    https://stackoverflow.com/questions/13712229/
    simultaneously-replacing-all-values-of-a-dictionary-to-zero-python


    Input: AA components of protein as a string
    Output: Number of AAs, molecular weight, breakdown of AA composition, 
            molar extinction coeff, mass exitinction coeff, & theoretical
            pI, all printed to the screen
    '''
    # These tables are for calculating:
    #     molecular weight (aa2mw), along with the mol. weight of H2O (mwH2O)
    #     absorbance at 280 nm (aa2abs280)
    #     pKa of positively charged Amino Acids (aa2chargePos)
    #     pKa of negatively charged Amino acids (aa2chargeNeg)
    #     and the constants aaNterm and aaCterm for pKa of the respective termini

    aa2mw = {
        'A': 89.093, 'G': 75.067, 'M': 149.211, 'S': 105.093, 'C': 121.158,
        'H': 155.155, 'N': 132.118, 'T': 119.119, 'D': 133.103, 'I': 131.173,
        'P': 115.131, 'V': 117.146, 'E': 147.129, 'K': 146.188, 'Q': 146.145,
        'W': 204.225, 'F': 165.189, 'L': 131.173, 'R': 174.201, 'Y': 181.189
    }

    mwH2O = 18.015
    aa2abs280 = {'Y': 1490, 'W': 5500, 'C': 125}

    aa2chargePos = {'K': 10.5, 'R': 12.4, 'H': 6}
    aa2chargeNeg = {'D': 3.86, 'E': 4.25, 'C': 8.33, 'Y': 10}
    aaNterm = 9.69
    aaCterm = 2.34

    def __init__(self, protein, cytosine):
        '''Initialize ProteinParam class and clean input to only have valid AAs'''
        self.cytosine = cytosine  # if True, cytosine is under oxidizing conditions

        # 1. Clean input, so it only contains valid AAs
        self.inputString = ""
        self.cleanString = ""
        for char in protein:
            if char.isalpha():
                char = char.upper()
                if char in self.aa2mw:
                    self.inputString = self.inputString + char
        self.cleanString = self.inputString
        # 2. Build a dictionary
        self.aaDictTemp = dict()  # only has AAs in user input
        self.aaDictIncludeZeros = dict.fromkeys(self.aa2mw, 0)  # includes all AAs
        for aa in self.cleanString:
            self.aaDictTemp[aa] = self.cleanString.count(aa)
            self.aaDictIncludeZeros[aa] = self.cleanString.count(aa)
        self.aaDict = self.aaDictTemp

    def aaCount(self):
        '''Count and return number of AAs from cleaned string (which only has valid AAs)'''
        return len(self.cleanString)

    def pI(self):
        '''Use binary search and _charge_() method to find theoretical pI.
        Capped at 100 iterations to prevent an infinite loop.'''
        self.searchRange = [0, 14.01]
        self.endFlag = False
        self.loopCounter = 0
        while self.endFlag == False:
            self.loopCounter += 1
            self.middleValue = ((self.searchRange[0] + self.searchRange[1]) / 2)
            self.chargeAtTopPH = self._charge_(self.searchRange[1])
            self.chargeAtBottomPH = self._charge_(self.searchRange[0])
            self.chargeAtTestPH = (self._charge_(self.middleValue))
            if round(self.chargeAtTestPH, 5) > 0:
                # move the floor up
                self.newFloor = self.middleValue
                self.newCeiling = self.searchRange[1]  # stays the same
                self.searchRange.clear()
                self.searchRange = [self.newFloor, self.newCeiling]
            elif round(self.chargeAtTestPH, 5) < 0:
                # move the ceiling down
                self.newFloor = self.searchRange[0]  # stays the same
                self.newCeiling = self.middleValue
                self.searchRange.clear()
                self.searchRange = [self.newFloor, self.newCeiling]
            else:
                self.endFlag = True
                return self.middleValue
            if self.loopCounter >= 100:
                # stop it from iterating forever
                print("Warning: Reached maximum number of iterations", end='')
                print(" for pI; reporting closest value after 100 iterations.")
                return self.middleValue

    def aaComposition(self):
        '''Return the dictionary of valid AAs in protein (ie, composition)'''
        return self.aaDictIncludeZeros
    

    def _charge_(self, pH):
        '''Return net charge of AA string at a given environmental pH'''
        self.pH = pH
        # First, calculate R, K, H, and N-terminus component
        self.iterPosList = ['R', 'K', 'H', 'N-terminus']
        self.posCharge = 0
        for aa in self.iterPosList:
            try:
                self.nAA = self.aaDictIncludeZeros[aa]
                self.pKA = self.aa2chargePos[aa]
            except KeyError:  # N-terminus
                self.nAA = 1
                self.pKA = self.aaNterm
            self.thisCharge = (self.nAA * ((10 ** self.pKA) / (10 ** self.pH + 10 ** self.pKA)))
            self.posCharge += self.thisCharge
        # Now, calculate D, E, C, Y, and C-terminus component
        self.iterNegList = ['D', 'E', 'C', 'Y', 'C-terminus']
        self.negCharge = 0
        for aa in self.iterNegList:
            try:
                self.nAA = self.aaDictIncludeZeros[aa]
                self.pKA = self.aa2chargeNeg[aa]
            except KeyError:  # C-terminus
                self.nAA = 1
                self.pKA = self.aaCterm
            self.thisCharge = (self.nAA * ((10 ** self.pH) / (10 ** self.pH + 10 ** self.pKA)))
            self.negCharge += self.thisCharge
        return self.posCharge - self.negCharge


    def molarExtinction(self):
        '''Calculate molar extinction coeff: E=(Y*E_{y})+(W*E_{W})+(C*E_{C})
        Can be calculated under oxidizing (default) or reducing conditions
        by changing the value of self.cytosine'''
        self.tyrosinePart = (self.aaDictIncludeZeros['Y']) * self.aa2abs280['Y']
        self.tryptophanPart = (self.aaDictIncludeZeros['W']) * self.aa2abs280['W']
        self.cytosinePart = int()
        if self.cytosine == True:
            self.cytosinePart = (self.aaDictIncludeZeros['C']) * self.aa2abs280['C']
        else:
            self.cytosinePart = 0
        self.molarExtinctionCoeff = self.tyrosinePart + self.tryptophanPart + self.cytosinePart
        return self.molarExtinctionCoeff


    def massExtinction(self):
        '''Calculate mass extinction using molar weight'''
        myMW = self.molecularWeight()
        return self.molarExtinction() / myMW if myMW else 0.0


    def molecularWeight(self):
        '''Calculate the molecular weight as sum MW of every AA in the dictionary 
        and subtracting water. Correctly handles instances where only one or no
        valid AAs exist.'''
        self.runningTotal = 0
        for key in self.aaDict:  # more efficient than self.aaDictIncludeZeros
            self.runningTotal += (self.aa2mw[key] - self.mwH2O) * self.aaDict[key]
        if self.aaCount() == 0:  # no valid AAs, so don't count MW of water
            return self.runningTotal
        else:
            return self.mwH2O + self.runningTotal
