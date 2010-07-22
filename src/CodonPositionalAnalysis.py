from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

import string
import logging.config
import types
import os
import errno
from ConfigParser import ConfigParser
from ConfigParser import NoSectionError
import sys

from typeutils.TypeChecker import require
from typeutils.Exceptions import FileError
from mappings import NUCLEOTIDE_AMBIGUITY_MAPPING

LOGGER_CONFIG_FILE = os.path.join(os.path.expanduser('~'), '.config', 'HIVMutationAnalysis', 
                                  'logging.conf')

## @brief Class for performing positional analysis of codons.
#
# This class reads codons from FASTA files provided it as an argument, and
# produces a list of proteins, with ambiguities written in parantheses. 
# @ingroup codonposanalysis
# @author Siddhu Warrier (siddhuwarrier@gmail.com)
# @date 22/02/2010.
class CodonPositionalAnalysis:
    @require(validKwargs = [], fastaData = types.DictionaryType)
    def __init__(self, fastaData):
        #Check Logger file
        #if the file does not exist
        if not os.path.exists(LOGGER_CONFIG_FILE):
            raise FileError(errno.ENOENT, "Logger File %s does not exist"\
                            %LOGGER_CONFIG_FILE)
        #if the path exists, but the path does not refer to a file
        if not os.path.isfile(LOGGER_CONFIG_FILE):
            raise FileError(errno.EINVAL, "Invalid path to logger file: %s"\
                            %LOGGER_CONFIG_FILE)
        
        #configure the logger
        try:
            logging.config.fileConfig(LOGGER_CONFIG_FILE)
        except NoSectionError as noSectionError:
            print noSectionError.args
            raise NoSectionError(LOGGER_CONFIG_FILE)
        
        self.logger = logging.getLogger('CodonPositionalAnalysis')
        self.logger.debug("Logger set up...")
        
        #error check input data
        for fastaFile in fastaData.keys():
            fastaFilename = string.translate(fastaFile, None, "\"") #remove ""s if present
            self.logger.debug("Fasta filename" +  fastaFilename)
            if not os.path.isfile(os.path.abspath(fastaFilename)):
                raise FileError(errno.ENOENT, "Fasta file " + fastaFile + " does not exist.")
            
            if type(fastaData[fastaFile])!= types.ListType:
                raise TypeError("Expected codon positions to be of type <list>, but was of type %s instead."%type(self.fastaData[fastaFile]))
        
        self.fastaData = fastaData
    
    def createProteinSequence(self):
        for fastaFile in self.fastaData.keys():
            fastaFilename = string.translate(fastaFile, None, "\"") #remove ""s if present
            seqIterator = SeqIO.parse(open(fastaFilename, "rU"), "fasta")
        
        proteins = []
        for record in seqIterator:
            proteins.append(self.makeProteinRecord(record))
            
        proteins = tuple(proteins)
        self.logger.debug("Proteins produced:")
        self.logger.debug(proteins)
        #SeqIO.write(proteins, open("../out/proteins.fasta", "w"), "fasta")
        
        seqIterator.close() 

    def makeProteinRecord(self, nucRecord):
        """Returns a new SeqRecord with the translated sequence (default table)."""
        proteinSequence = nucRecord.seq.translate(1,to_stop = False, 
                                                  stop_symbol = "*") #get the protein sequence
        
        posCount = 0 #keep track of the position in the protein sequence
    
        newProteinString = proteinSequence.__str__()
        
        for char in proteinSequence:
            #look for ambiguous amino acids.
            if char == "X":
                #self.logger.debug("Sequence ID: %s, Ambiguity in position %d"%(nucRecord.id, posCount))
                #self.logger.debug("The nucleotide sequence that caused the problem is: %s"
                 #                 %(nucRecord[posCount *3]+nucRecord[posCount *3 + 1]+
                 #                   nucRecord[posCount *3  + 2]))
                #resolve the ambiguity using the mapping table.
                ambiguousAminoAcidSequence = self.createAmbiguousProteinList(nucRecord, posCount *3)
                newProteinString = newProteinString.replace("X", ambiguousAminoAcidSequence)
            posCount += 1
            
        proteinSequence = Seq(newProteinString)
        
        return SeqRecord(seq = proteinSequence, id= "trans_" + nucRecord.id)

    def createAmbiguousProteinList(self, nucRecord, firstNucPosition):        
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
    

def main():
    config = ConfigParser()
    openedFile = config.read(os.path.join(os.path.expanduser("~"),".config",
                                          "HIVMutationAnalysis", "HIVMutationAnalysis-config"))
    if openedFile == []:
        sys.stderr.write("File %s not found"%os.path.join(os.path.expanduser("~"), ".config", 
                                                          "HIVMutationAnalysis", "HIVMutationAnalysis-config"))
        print "Please create file."
        sys.exit(0)
    
    #get the fasta data into a dictionary
    fastaData = {}
    for section in config.sections():
        if section == "fasta-file":
            fastaData[config.get("fasta-file", "filename")] = eval(config.get("fasta-file", "codons")) 
            
    
    CodonPositionalAnalysis(fastaData).createProteinSequence()
    
if __name__ == "__main__":
    main()