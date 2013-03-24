'''
Created on Mar 9, 2013

@author: arenduchintala
'''

import sys, os, re

wordCounts = {}
emissionCounts = {}
def oldRareWords(token):
    if len(re.findall(r'[0-9]', token)) > 0:
        return '_RARE_NUMERIC_'
    elif  len(re.findall(r'[A-Z]', token)) == len(token):
            return '_RARE_ALL_CAPS_'
    elif len(re.findall(r'[A-Z]', token[len(token) - 1])) == 1:
            return '_RARE_LAST_CAP_'
    else: 
        return '_RARE_'
        
def replaceRareWord(token):
    return oldRareWords(token)
    # try:
    #    float(token)
    #    return '_RARE_NUMERIC_'
    # except:
    #    
    #    if (token.upper() == token):
    #        return '_RARE_ALL_CAPS_'
    #    elif (token[len(token)-1].upper() == token[len(token)-1]):
    #        return '_RARE_LAST_CAP_'
    #    elif (token[0].upper() == token[0]):
    #        return '_RARE_FIRST_CAP_'
    #    elif (token.endswith('s')):
    #        return '_RARE_END_S_'
    #    else:
    #        return '_RARE_'
        
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
    
    parseCountFile('../data/gene.count')
    originalFile = open('../data/gene.train', 'r')
    originalTrainingLines = originalFile.readlines()
    modLines = originalTrainingLines
    originalFile.close()
    
    writer = open('../outputs/gene-rare-cat.train', 'w')
    for i in range (len(originalTrainingLines)):
        
        if  originalTrainingLines[i] != '\n':
           
            originalTrainingLines[i] = originalTrainingLines[i].strip()
            currentToken = originalTrainingLines[i].split(" ")[0].strip()
            currentTokenTag = originalTrainingLines[i].split(" ")[1].strip()
            if (wordCounts[currentToken] < 5):
                modLines[i] = replaceRareWord(currentToken) + " " + currentTokenTag + '\n'
            else:
                modLines[i] = currentToken + ' ' + currentTokenTag + '\n'
        else:
            modLines[i] = '\n'
        writer.write(modLines[i])
    writer.flush();
    writer.close();
    
    os.system("python ../src/count_freqs.py ../outputs/gene-rare-cat.train > ../outputs/gene-rare-cat.count")
