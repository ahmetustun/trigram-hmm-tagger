'''
Created on Mar 9, 2013

@author: arenduchintala
'''
import sys,re

wordCounts = {}
emissionCounts = {}
unigramCounts = {}
bigramCounts = {}
trigramCounts = {}

def usage():
    print "python MLEstimator [file with emmision and n-gram counts] [output modified gene.train file]"

def parseCountFile(inputFile):
    lines = open(inputFile, 'r').readlines()
    for line in lines:
        items = line.strip().split(" ")
        if (items[1].strip() == 'WORDTAG'):
            emissionCounts[items[2] + ',' + items[3]] = int(items[0])
            if wordCounts.has_key(items[3]):
                wordCounts[items[3]] = wordCounts[items[3]] + int(items[0])
            else:
                wordCounts[items[3]] = int(items[0])
                
        elif (items[1].strip() == '1-GRAM'):
            unigramCounts[items[2]] = int(items[0])
        elif (items[1].strip() == '2-GRAM'):
            bigramCounts[items[2] + ',' + items[3]] = int(items[0])
        elif (items[1].strip() == '3-GRAM'):
            trigramCounts[items[2] + ',' + items[3] + ',' + items[4]] = int(items[0])


    
if __name__ == "__main__":

    if len(sys.argv) == 3 or len(sys.argv) == 2:  # Expect exactly one argument: the training data file
        parseCountFile(sys.argv[1])
        if len(sys.argv) == 3:
            originalTrainingData = open(sys.argv[2], 'r').read()
            for word in sorted(wordCounts.keys()):
                if wordCounts[word] < 5:
                    rareWord = re.compile(r"'\b'+word+'\b'", re.IGNORECASE)
                    originalTrainingData = rareWord.sub( 'RARE', originalTrainingData)
        print originalTrainingData
    else:
        usage()
