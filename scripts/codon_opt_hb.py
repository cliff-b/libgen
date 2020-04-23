#!/programs/anaconda/bin/python
from Bio import SeqIO
import random
from Bio.Alphabet import generic_dna
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import Restriction
import copy
import re
import sys
random.seed(4594395793653479834)
xspacer = 'AgaagagcTGCCTCAATCTAGTCTTCATTCGCAGACCTAGCCGGATCGACCGTGCTACCGTTTCCGTAAGCCGAAACCGGAGAGCAAGAGGGCGTCATACC'
x_primers = {1: ['AGCGAAACCGTGCGT', 'TCTGGGTGCGCATCC'], 2: ['TGTCCCAGGTCGCAG', 'CCGTGCCACAACTCG'], 3: ['TCGCGGAGTTGAGGT', 'CTACGGCTCGGCAAC']}
yspacer = 'TGACTCCACAGATAGAACCGATCTCGCACTCGGAGGCCCAGAGTGTTTGACTCGTCCGGTACGTACTTTGCCTGGTAGAAATCGGTAACCGAAgcagt'
y_primers = {1: ['ACCTCCTGCGGCATT', 'TGGCTGAATCGCGCT'], 2: ['GACGTGCGGAGGTGA', 'TGTGCTCGTGGTGGA'], 3: ['ACGCCTGTGTCTGGT', 'GGCTGCGACGGTCTT']}


def readcodonusage(filename):
    # reads in usage file with probabilities
    handle = open(filename)
    strings = []
    for line in handle:
        splitline = line.split()
        for i in [0, 5, 10, 15]:
            strings.append(splitline[i:i + 5])

    codons = {}
    for string in strings:
        if string[1] in codons:
            dnacodon = string[0].replace('U', 'T')
            codons[string[1]].append([dnacodon, float(string[3])])
        else:
            dnacodon = string[0].replace('U', 'T')
            codons[string[1]] = [[dnacodon, float(string[3])]]

    # make a fractional accounting
    for aa in codons:
        runningaacount = 0
        for i in codons[aa]:
            runningaacount += i[1]
        for i in codons[aa]:
            fraction = i[1] / runningaacount
            i.append(fraction)
    # drop codons less than 15%
    for aa in codons:
        newlist = []
        for i in codons[aa]:
            if i[2] > 0.15:
                newlist.append(i)
        codons[aa] = newlist

    # do running total accounting
    for aa in codons:
        runningaacount = 0
        runningtotal = 0
        for i in codons[aa]:
            runningaacount += i[1]
        for i in codons[aa]:
            fraction = i[1] / runningaacount
            i[2] = fraction
            runningtotal += fraction
            i.append(runningtotal)

    # take last codon from each and make 1.0
    for aa in codons:
        codons[aa][-1][-1] = 1.0

    # make reverse codon dictionary
    revcodons = {}
    for aa in codons:
        for i in codons[aa]:
            revcodons[i[0]] = aa
    return codons, revcodons


def obtainfastasequences(filename):
    """gets fasta sequences from filename and turns them into a list of seq objects"""
    handle = open(filename)
    records = []
    for seqrecord in SeqIO.parse(handle, "fasta"):
        records.append(seqrecord)
    return records


def codonoptimize(codons, aarecords, usage_no):
    """creates a codon optimized DNA sequence from a protein sequence"""
    newrecords = []
    for record in aarecords:
        sequence = ""
        for letter in record.seq:
            rn = random.random()
            for i in codons[letter]:
                if rn < i[-1]:
                    sequence = sequence + i[0]
                    break
        newseq = Seq(sequence, generic_dna)
        newrecord = SeqRecord(newseq, id=record.id + '-C' + str(usage_no))
        newrecords.append(newrecord)
    return newrecords


def removesite(seqrecord, site, codons, revcodons):
    """remove restriction sites and repeated nucleotides by putting in a new codon"""
    index = seqrecord.seq.find(site)
    recognitionlength = len(site)
    # calculate how many codons can be changed and then change the first of them
    codonstochange = int((recognitionlength - 1) / 3) + 1
    for i in range(codonstochange):
        codonstart = i * 3 + index - (index % 3)
        codontochange = seqrecord.seq[codonstart:codonstart + 3]
        for possiblecodon in codons[revcodons[str(codontochange)]]:
            if possiblecodon[0] != str(codontochange):
                seqrecord.seq = seqrecord.seq[:codonstart] + possiblecodon[0] + seqrecord.seq[codonstart + 3:]
                return
    print("Could not change restriction site")


def changerestrictionsites(seqrecords, codons, revcodons):
    """uses remove site function to change restriction enzyme sites depending on location"""
    for seqrecord in seqrecords:
        rb = Restriction.RestrictionBatch([Restriction.AscI, Restriction.BspQI, Restriction.PstI, Restriction.EcoRI, Restriction.NotI, Restriction.BtsI, Restriction.BsaI])
        reanalysis = rb.search(seqrecord.seq)
        for key in reanalysis:
            for _ in reanalysis[key]:
                seqkey = Seq(key.site, generic_dna)
                removesite(seqrecord, seqkey, codons, revcodons)
    return seqrecords


def checkrepeats(seqrecords, codons, revcodons):
    """ check for repeated nucleotides, and if so, remove them and add new codon"""
    for seqrecord in seqrecords:
        patt = re.compile(r"([ATGCU])\1{3,}")
        repeats = patt.finditer(str(seqrecord.seq))
        for i in repeats:
            seqkey = Seq(i.group(), generic_dna)
            removesite(seqrecord, seqkey, codons, revcodons)


def addbuffers(seqrecords, f_buff, r_buff):
    """adds buffer regions to the ends of the sequence"""
    for seqrecord in seqrecords:
        seqrecord.seq = f_buff + seqrecord.seq + r_buff


def reversecomprecord(seqrecords):
    """reverses complements a sequence and changes its name to reflect that"""
    for seqrecord in seqrecords:
        seqrecord.seq = seqrecord.seq.reverse_complement()
        seqrecord.id = seqrecord.id[:-3] + "_rev" + seqrecord.id[-3:]


def addspacer(seqrecords, spacer, size, direction):
    """adds spacer DNA of arbitrary length. Can be added in either direction"""
    for seqrecord in seqrecords:
        if direction == -1:
            seqrecord.seq = spacer[direction * (size - len(seqrecord.seq)):] + seqrecord.seq
        else:
            seqrecord.seq = seqrecord.seq + spacer[:direction*(size - len(seqrecord.seq))]


# ============================================================================================================
if __name__ == '__main__':
    # read in codon optimization table and proteins to convert to DNA
    codondic, revcodondic = readcodonusage(sys.argv[2])
    aarecordseqs = obtainfastasequences(sys.argv[1])
    for usage in range(1, 4):
        # optimize codons and remove RE sites it makes
        seqoptimized = codonoptimize(codondic, aarecordseqs, usage)
        updatedseqrecords = changerestrictionsites(seqoptimized, codondic, revcodondic)
        # remove repeated nucleotides, must do several times because only matches first repeat seq per line
        for _ in range(3):
            checkrepeats(seqoptimized, codondic, revcodondic)

        # Recheck that we didn't create restriction sites removing repeated nucleotides
        updatedseqrecords = changerestrictionsites(updatedseqrecords, codondic, revcodondic)
        # copy seq records for reverse complement
        seqrecordscopy = copy.deepcopy(updatedseqrecords)

        # add buffersites and finalize sequences
        addbuffers(updatedseqrecords, "gTG", "CCAGAgagacc")
        addspacer(updatedseqrecords, yspacer, 199, -1)
        addbuffers(updatedseqrecords, y_primers[usage][0], y_primers[usage][1])
        addbuffers(updatedseqrecords, "", 'g')

        # Get reverse complement of sequences for the "X" libraries
        reversecomprecord(seqrecordscopy)
        # #add buffers and primers for "X" seqs
        addbuffers(seqrecordscopy, "TTA", "agg")
        addspacer(seqrecordscopy, xspacer, 194, 1)
        addbuffers(seqrecordscopy, x_primers[usage][0], x_primers[usage][1])
        addbuffers(seqrecordscopy, "", "ttcttc")

        for seq in updatedseqrecords:
            print(">" + seq.id + "," + seq.seq)
        for seq in seqrecordscopy:
            print(">" + seq.id + "," + seq.seq)
