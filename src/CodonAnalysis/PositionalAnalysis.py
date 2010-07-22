'''
Created on 22 Jul 2010

@author: siddhu
'''
from typeutils.TypeChecker import require
from typeutils.Exceptions import FileError

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import single_letter_alphabet

from mappings import NUCLEOTIDE_AMBIGUITY_MAPPING

import types
import logging

##Based on BioPython code
#Port of Fasta Iterator from BioPython, removing bits I don't need and cleaning
#up the code
#This function is a generator
def fastaIterator(handle):
    #logging
    logging.getLogger().setLevel(logging.DEBUG)
    
    #Skip any text before the first record (e.g. blank lines, comments)
    while True:
        line = handle.readline()
        if line == "" :
            raise ValueError("Invalid file. Fasta File should contain at least one record.") 
            return #Premature end of file, or just empty?
        if line[0] == ">":
            break
    
    while True:
        if line[0]!=">":
            raise ValueError("Records in Fasta files should start with '>' character")
        else:
            descr = line[1:].rstrip()
            id   = descr.split()[0]
            name = id

        
        #a single fasta sequence might be split across multiple lines
        #read and create a list of them
        lines = []
        line = handle.readline() #jump blithely past the description line, on to
        #the next, and start looping
        while True:
            #if nothing read, break
            if not line : 
                break
            #if the next sequence is encountered, break
            if line[0] == ">": 
                break 
            
            #Remove trailing whitespace, and any internal spaces
            #(and any embedded \r which are possible in mangled files
            #when not opened in universal read lines mode)
            lines.append(line.rstrip().replace(" ","").replace("\r",""))
            line = handle.readline()
    
        #Return the record and then continue...
        record =SeqRecord(Seq("".join(lines), single_letter_alphabet),
                        id = id, name = name, description = descr)
        
        yield SeqRecord(Seq("".join(lines), single_letter_alphabet),id = id, name = name, 
                        description = descr)
        

        if not line : 
            return #StopIteration

##Class to perform fasta analysis    
class FastaAnalysis(object):
    logging.getLogger().setLevel(logging.DEBUG)
    @require(validKwargs = [], requiredCodonPositions = types.ListType)
    def __init__(self, handle, requiredCodonPositions):
        self.handle = handle
        #save the records in file obtained from the generator
        
        self.recordsInFile = fastaIterator(self.handle)
        
        self.requiredCodonPositions = requiredCodonPositions
            
    ##create a protein sequences out of the nucleotide sequence
    def createProteinSequences(self):
        codonsPerSequence = {}
        
        for nucRecord in self.recordsInFile:
            codonsPerSequence[nucRecord.id] =  self._makeProteinRecord(nucRecord)
        
        return codonsPerSequence
    
    ##protected function to convert a nucleotide sequence into a protein 
    #sequence
    def _makeProteinRecord(self, nucRecord):
        codons = []
        """Returns a list of protein codons with the translated sequence 
        (default table)."""
        #get the protein sequence
        proteinSequence = nucRecord.seq.translate(1,to_stop = False, 
                                                  stop_symbol = "*") 
        
        posCount = 0 #position in protein sequence count
        for char in proteinSequence:
            if posCount in self.requiredCodonPositions:
                #biopython converts ambiguous protein sequences into X
                # we shall do the mapping ourselves
                if char == 'X':
                    char = self._createAmbiguousProteinList(nucRecord, 
                                                           posCount *3)
                codons.append(char)
            
            posCount += 1
        
        return codons
    
    ##Really nasty function written in middle of night. But it works.
    #Replace with something better later.
    def _createAmbiguousProteinList(self, nucRecord, firstNucPosition):        
        nucleotideSequence = Seq("%s"%(nucRecord[firstNucPosition]+nucRecord[firstNucPosition + 1]
                                        +nucRecord[firstNucPosition + 2]))
        
        altNucleotideSeqs = []
        altNucleotideString = ""
        for nucleotide in nucleotideSequence:    
            if nucleotide.lower() in NUCLEOTIDE_AMBIGUITY_MAPPING.keys():
                #if there are no alternate nucleotide sequences yet, make them
                if altNucleotideSeqs == []:
                    for ambiguousNucleotide in NUCLEOTIDE_AMBIGUITY_MAPPING.get(nucleotide.lower()):
                        altNucleotideString = nucleotideSequence.__str__()
                        altNucleotideString = altNucleotideString.replace(nucleotide, ambiguousNucleotide.upper())
                        altNucleotideSeqs.append(Seq(altNucleotideString))
                else: #this is not hte first ambiguous nucleotide
                        tempAltNucleotideSeqs = []
                        for altNucleotideSeq in altNucleotideSeqs:
                            for ambiguousNucleotide in NUCLEOTIDE_AMBIGUITY_MAPPING.get(nucleotide.lower()):
                                altNucleotideString = altNucleotideSeq.__str__()
                                altNucleotideString = altNucleotideString.replace(nucleotide, ambiguousNucleotide.upper())
                                tempAltNucleotideSeqs.append(Seq(altNucleotideString))
                                
                            #remove the original altNucleotideSeq which has been replaced now
                        
                        #replace the existing list with the new list
                        altNucleotideSeqs = tempAltNucleotideSeqs
        ambiguousAminoAcidSeq = ""
        for altNucleotideSeq in altNucleotideSeqs:
            ambiguousAminoAcidSeq += altNucleotideSeq.translate().__str__() + "/"
        
        ambiguousAminoAcidSeq = ambiguousAminoAcidSeq[:-1] #remove the trailing slash
        
        ambiguousAminoAcidSeq = "(" + ambiguousAminoAcidSeq + ")"
        
        #self.logger.debug("Ambiguous amino acid seq resolved to " + ambiguousAminoAcidSeq)
        return ambiguousAminoAcidSeq
        
        
        
    