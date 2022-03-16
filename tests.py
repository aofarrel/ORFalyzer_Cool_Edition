import os
import unittest

from orfalyzer import sequenceanalysis
import orfalyzer

class TestWeirdFastas(unittest.TestCase):

	# lateHeader: header row is not the first line
	def test_fasta_lateHeader(self):
		lateHeader = "asdf\nlalalalala\n>now I'm a header\nATGCCCCCCCCCTAA"
		with open('lateHeader.fasta', 'w') as f: f.write(lateHeader)
		lateReader = orfalyzer.sequenceanalysis.FastaReader(fname='lateHeader.fasta')
		for head, seq in lateReader.readFasta():
			assert head == "now I'm a header"
			assert seq == "ATGCCCCCCCCCTAA"
		os.remove('lateHeader.fasta')

	# multiHeader: has multiple headers
	def test_fasta_multiHeader(self):
		n = 0
		multiHeader = ">first header\nATGCCCCCCCCCTAA\n>second header\nATGGGGGGGGGTAA"
		with open('multiHeader.fasta', 'w') as f: f.write(multiHeader)
		multiReader = orfalyzer.sequenceanalysis.FastaReader(fname='multiHeader.fasta')
		for head, seq in multiReader.readFasta():
			assert head == ("first header" if n ==0 else "second header")
			assert seq == ("ATGCCCCCCCCCTAA" if n ==0 else "ATGGGGGGGGGTAA")
			n += 1
		os.remove('multiHeader.fasta')

	# chevronPlus: has extra chevrons all over the place
	def test_fasta_chevronPlus(self):
		n = 0
		chevronPlus = ">chevron with >techron<\nATGCCCCCCC>CCTAA>\n>wow this file is a mess\nATGCCCCCCCCCTAA"
		with open('chevronPlus.fasta', 'w') as f: f.write(chevronPlus)
		chevronReader = orfalyzer.sequenceanalysis.FastaReader(fname='chevronPlus.fasta')
		for head, seq in chevronReader.readFasta():
			assert head == ("chevron with >techron<" if n ==0 else "wow this file is a mess")
			assert seq == ("ATGCCCCCCC>CCTAA>" if n ==0 else "ATGCCCCCCCCCTAA")
			n += 1
		os.remove('chevronPlus.fasta')

	# emptyHeader: header is empty
	def test_fasta_emptyHeader(self):
		emptyHeader = ">\nATGCCCCCCCCCTAA"
		with open('emptyHeader.fasta', 'w') as f: f.write(emptyHeader)
		emptyReader = orfalyzer.sequenceanalysis.FastaReader(fname='emptyHeader.fasta')
		for head, seq in emptyReader.readFasta():
			assert head == ""
			assert seq == "ATGCCCCCCCCCTAA"
		os.remove('emptyHeader.fasta')

class TestWeirdReadingFrames(unittest.TestCase):
	global testID
	testID = 0
	settings = [-1,True,"temp%d.txt" % testID,['ATG'],['TAG', 'TGA', 'TAA'],True]

	def setupLyzerObject(self, orfToTest, header):
		lyzer = orfalyzer.sequenceanalysis.OrfFinder()
		global testID
		lyzer._processArgs_(*self.settings)
		lyzer._processHeader_(header, 0)
		lyzer._addsequence_(orfToTest)
		testID += testID
		return lyzer

	#no ATG -- should still be considered an ORF as the ATG may have come earlier
	def test_orf_noATG(self):
		noATG = "TAATCTTTTAAAGGGCCCTTTTAAAATC"
		noATGLyzer = self.setupLyzerObject(noATG, "noATG")
		global testID
		with open('temp%d.txt' % testID, 'r') as f: test = f.readlines()
		with open('testTruthFiles/test%d.txt' % testID, 'r') as g: truth = g.readlines()
		assert (test==truth)
		os.remove('temp%d.txt' % testID)
		testID += 1

	# these tests were never fully implemented as they have no truth files
	#
	#
	#short ORF
	# def test_orf_shortORF():
	# 	shortORF = "ATGTAA"
	# 	shortLyzer = setupLyzerObject(shortORF, "shortORF")
	# 	global testID
	# 	with open('temp%d.txt' % testID, 'r') as f: test = f.readlines()
	# 	with open('testTruthFiles/test%d.txt' % testID, 'r') as g: truth = g.readlines()
	# 	assert (test==truth)
	# 	os.remove('temp%d.txt' % testID)
	# 	testID += 1

	# # buncha newlines
	# def test_orf_newlines():
	# 	newlines = ">newlines\nATGATGATGATG\nTAANNNNNNNNNNNNATG\nCCCGGGCCCGGGTAA"

	# # ORF only contains unknown nucleotides
	# def test_orf_unknownORF():
	# 	unknownORF = ">unknownORF\nATGNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNTAA"

	# #no ORF
	# def test_orf_noORF():
	# 	noORF = "NNNNNNNNNNNN"
	# 	noLyzer = setupLyzerObject(noORF, "noORF")

if __name__ == '__main__':
    unittest.main()


