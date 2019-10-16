import Bio.Blast.NCBIWWW as BioWeb
import operator
from Bio.Seq import Seq
import pandas as pd
pd.set_option('display.max_rows', 500)
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 1000)
pd.set_option('display.max_colwidth',40)

# set options here
filePath = r"C:\Users\Sulley\Downloads\dna2.fasta"
repeatNumber = 12 # Set number of repeats to analyze here!!!
# --------

seqDF = pd.DataFrame(columns=["geneID", "geneInfo", "Sequence"])
mergedSeq = ""
with open(filePath) as seqFile:
    entryNumber = -1
    for line in seqFile:
        if line.startswith(">"):
            entryNumber += 1
            ID = line[1:].split()[0]
            info = line[1:]
            seqDF.loc[entryNumber] = [ID,info,""] #Add ID and info leaving sequence data emtpy for now
        else:
            seqDF.iloc[entryNumber,2] += line.rstrip()
            mergedSeq += line.rstrip() # Get one big sequence to make repeat calculations on all sequences (kinda pointless I think)

#How many sequences in file
print(f"Number of records in file: {entryNumber + 1}")

def sequenceLengths():
    # Prints longest and shortest Sequences of a given multifasta file
    seqDF["SeqLength"] = seqDF["Sequence"].str.len() #Add sequence length column

    longestSeqValue = seqDF["SeqLength"].max() # Find longest sequence
    longestSeqsList = seqDF.index[seqDF["SeqLength"]==longestSeqValue].tolist() # Make list of all sequences with longest sequences (if more than one)
    for longestSeqIndex in longestSeqsList: #iterate through indices and print out row values
        longestSeqName, longestSeqLength = seqDF.iloc[longestSeqIndex,0],seqDF.iloc[longestSeqIndex,3]
        print(f"Longest Sequence(s): {longestSeqName}, {longestSeqLength} bp")

    shortestSeqValue = seqDF["SeqLength"].min() # Find shortest sequence
    shortestSeqsList = seqDF.index[seqDF["SeqLength"]==shortestSeqValue].tolist() # Make list of all sequences with longest sequences (if more than one)
    for shortestSeqIndex in shortestSeqsList: #iterate through indices and print out row values
        shortestSeqName, shortestSeqLength = seqDF.iloc[shortestSeqIndex,0],seqDF.iloc[shortestSeqIndex,3]
        print(f"Shortest Sequence(s): {shortestSeqName}, {shortestSeqLength} bp")
sequenceLengths()

def longestORF(sequence):
    # Determines longest ORF of a sequence (must include M and stop) and returns length of ORF in nt and start position within sequence
    # Input can be string (or Biopython sequence(?))

    try: # make biopython sequence of input sequence if string
        inputSeq = Seq(sequence)
    except: # otherwise likely already biopython sequence
        inputSeq = sequence

    frameLength = 0
    longestFrame = 0
    longestFrameStartPos = 0
    for readingFrame in range(3):
        aaSeq = inputSeq[readingFrame:].translate()
        start = False
        for aaPos, aa in enumerate(aaSeq):
            if aa =="M" and not start:
                start = True
                startPos = 3*aaPos + readingFrame + 1
                frameLength = 3
            elif aa == "*":
                if start:
                    if longestFrame < frameLength:
                        longestFrame = frameLength
                        longestFrameStartPos = startPos
                    start = False
            else:
                frameLength += 3
    return longestFrameStartPos, longestFrame

def ORFs():
    # Adds longest ORF column to each sequence of sequence data frame of the multifasta file. Then prints overall longest
    # ORF (that includes Met and stop codon)

    seqDF["longestORFStartPos"], seqDF["longestORFLength"] = seqDF["Sequence"], seqDF["Sequence"] # Duplicate sequence to then apply function
    seqDF["longestORFStartPos"] = seqDF["longestORFStartPos"].apply(longestORF) #apply get longest ORF function to copied sequences
    seqDF["longestORFLength"] = seqDF["longestORFStartPos"].str[1] # put length and start pos in different columns as function returns tuples
    seqDF["longestORFStartPos"] = seqDF["longestORFStartPos"].str[0]

    absLongestORF = seqDF["longestORFLength"].idxmax()
    print(f"The longest ORF in the file is in gene {seqDF.iloc[absLongestORF,0]} and is {seqDF.iloc[absLongestORF,5]}bp long"
    f", starting at position {seqDF.iloc[absLongestORF,4]}")

def repeats(sequence, repeatLength):
    # Counts all repeats of a given length in given sequences. Returns dictionary with repeat info and prints out the
    # most frequent repeat and number of repeats

    repeatsDict = {}
    for nt in range(len(sequence)-repeatLength): # limit so last call's repeat will end at sequence end
        repeatSeq = sequence[nt:nt+repeatLength]
        if repeatSeq in repeatsDict:
            repeatsDict[repeatSeq] += 1
        else:
            repeatsDict[repeatSeq] = 1
    repeatsByVal = sorted(repeatsDict.items(), key=operator.itemgetter(1), reverse=True)
    return repeatsByVal

def addRepeatsInfo(repeatLength):
    # adds column to sequence data frame with most common repeat of x length

    seqDF[f"{repeatLength}bp repeats"] = seqDF["Sequence"] # Duplicate sequence to then apply function
    seqDF[f"{repeatLength}bp repeats"] = seqDF[f"{repeatLength}bp repeats"].apply(repeats, args=(repeatLength,)) # apply get repeats info for given repeat length

fileRepeats = repeats(mergedSeq, repeatNumber)
print(f"The most repeats of {repeatNumber}bp in the entire file is {fileRepeats[0][0]}:{fileRepeats[0][1]} repeats")

ORFs()
addRepeatsInfo(repeatNumber)
print("----------------------")
print(seqDF)


