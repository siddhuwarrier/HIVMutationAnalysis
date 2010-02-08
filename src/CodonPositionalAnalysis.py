from Bio import SeqIO
from Bio.Data import CodonTable
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from mappings import NUCLEOTIDE_AMBIGUITY_MAPPING
import sys

def createAmbiguousProteinList(nucRecord, firstNucPosition):
    
    print "%s"%(nucRecord[firstNucPosition]+nucRecord[firstNucPosition + 1]
                                        +nucRecord[firstNucPosition + 2])
    nucleotideSequence = Seq("%s"%(nucRecord[firstNucPosition]+nucRecord[firstNucPosition + 1]
                                    +nucRecord[firstNucPosition + 2]))
    
    nucleotideSequence = Seq("AAK")
    altNucleotideSeqs = []
    altNucleotideString = ""
    for nucleotide in nucleotideSequence:    
        if nucleotide.lower() in NUCLEOTIDE_AMBIGUITY_MAPPING.keys():
            print "ambiguous nucleotide:", nucleotide
            #if there are no alternate nucleotide sequences yet, make them
            if altNucleotideSeqs == []:
                for ambiguousNucleotide in NUCLEOTIDE_AMBIGUITY_MAPPING.get(nucleotide.lower()):
                    altNucleotideString = nucleotideSequence.__str__()
                    altNucleotideString = altNucleotideString.replace(nucleotide, ambiguousNucleotide.upper())
                    altNucleotideSeqs.append(Seq(altNucleotideString))
                    print altNucleotideString
            else: #this is not hte first ambiguous nucleotide
                    tempAltNucleotideSeqs = []
                    for altNucleotideSeq in altNucleotideSeqs:
                        print "Existing seq:" + altNucleotideSeq.__str__()
                        for ambiguousNucleotide in NUCLEOTIDE_AMBIGUITY_MAPPING.get(nucleotide.lower()):
                            altNucleotideString = altNucleotideSeq.__str__()
                            altNucleotideString = altNucleotideString.replace(nucleotide, ambiguousNucleotide.upper())
                            tempAltNucleotideSeqs.append(Seq(altNucleotideString))
                            print altNucleotideString
                        #remove the original altNucleotideSeq which has been replaced now
                    
                    #replace the existing list with the new list
                    altNucleotideSeqs = tempAltNucleotideSeqs
                            
    for altNucleotideSeq in altNucleotideSeqs:
        print altNucleotideSeq.translate() + "/",
    
       

def makeProteinRecord(nucRecord):
    """Returns a new SeqRecord with the translated sequence (default table)."""
    proteinSequence = nucRecord.seq.translate(1,to_stop = False, 
                                              stop_symbol = "*") #get the protein sequence
    
    posCount = 0 #keep track of the position in the protein sequence
    for char in proteinSequence:
        #look for ambiguous amino acids.
        if char == "X":
            print "Sequence ID: %s, Ambiguity in position %d"%(nucRecord.id, posCount)
            print "The nucleotide sequence that caused the problem is:",
            print "%s"%(nucRecord[posCount *3]+nucRecord[posCount *3 + 1]+nucRecord[posCount *3  + 2])
            createAmbiguousProteinList(nucRecord, posCount *3)
            sys.exit(0)
        posCount += 1
    
    
#    return SeqRecord(seq = nucRecord.seq.translate(1
#,to_stop = False, stop_symbol = "*"), \
#                     id = "trans_" + nucRecord.id, \
#                     description = "translation of CDS, using default table")

def main():
    seqIterator = SeqIO.parse(open("../input/all-blinded.fasta", "rU"), "fasta")
    
    for record in seqIterator:
        makeProteinRecord(record)
    
    #SeqIO.write(proteins, open("../out/proteins.fasta", "w"), "fasta")
    
    seqIterator.close()
if __name__ == "__main__":
    main()