'''
Created on Mar 9, 2013

@author: arenduchintala
'''

import sys,os,re

wordCounts = {}
emissionCounts = {}

def replaceRareWord(token):
    return '_RARE_'
    if re.match("[0-9]", token):
        return '_RARE_NUMERIC_'
    elif (token.upper() == token):
        return '_RARE_ALL_CAPS_'
    elif (token[len(token)-1].upper() == token[len(token)-1]):
        return '_RARE_LAST_CAP_'
    else:
        return '_RARE_'
    
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
    
    parseCountFile('outputs/gene.count')
    originalFile = open('data/gene.train', 'r')
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
                modLines[i] = replaceRareWord(currentToken)+" " + currentTokenTag + '\n'
            else:
                modLines[i] = currentToken + ' ' + currentTokenTag+'\n'
        else:
            modLines[i] = '\n'
        writer.write(modLines[i])
    writer.flush();
    writer.close();
    os.system("python src/count_freqs.py outputs/gene-rare.train > outputs/gene-rare.count")