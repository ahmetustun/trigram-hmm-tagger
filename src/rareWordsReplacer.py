'''
Created on Mar 9, 2013

@author: arenduchintala
'''

import sys,os,re

wordCounts = {}
emissionCounts = {}

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
                
    
if __name__ == "__main__":
    if len(sys.argv) == 3:  # Expect exactly one argument: the training data file
        parseCountFile(sys.argv[1])
        originalFile = open(sys.argv[2], 'r')
        originalTrainingLines = originalFile.readlines()
        modLines  = originalTrainingLines
        originalFile.close()
        
        writer = open('outputs/gene-rare.train','w')
        for i in range ( len(originalTrainingLines)):
            
            if  originalTrainingLines[i] != '\n':
               
                originalTrainingLines[i] = originalTrainingLines[i].strip()
                currentToken = originalTrainingLines[i].split(" ")[0].strip()
                currentTokenTag  = originalTrainingLines[i].split(" ")[1].strip()
                if ( wordCounts[currentToken] < 5):
                    modLines[i] = "_RARE_ " + currentTokenTag + '\n'
                else:
                    modLines[i] = currentToken + ' ' + currentTokenTag+'\n'
            else:
                modLines[i] = '\n'
            writer.write(modLines[i])
    
    writer.flush();
    writer.close();
