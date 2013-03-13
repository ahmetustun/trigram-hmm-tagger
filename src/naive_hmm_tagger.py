'''
Created on Mar 9, 2013

@author: arenduchintala
'''
import math
import sys

def MLEstimates(countFilePath):
    emissionCounts = {}
    seenObservations =  {}
    unigramCounts = {}
    bigramCounts = {}
    trigramCounts = {}


    lines = open(countFilePath, 'r').readlines()
    for line in lines:
        items = line.strip().split(" ")
        if (items[1].strip() == 'WORDTAG'):
            # 8 WORDTAG O localizes
            emissionCounts[items[2] + '->' + items[3]] = int(items[0])
        elif (items[1].strip() == '1-GRAM'):
            # 749 1-GRAM I-GENE
            unigramCounts[items[2]] = int(items[0])
        elif (items[1].strip() == '2-GRAM'):
            # 749 2-GRAM * I-GENE
            bigramCounts[items[2] + '->' + items[3]] = int(items[0])
        elif (items[1].strip() == '3-GRAM'):
            # 1813 3-GRAM I-GENE O STOP
            trigramCounts[items[2] + '->' + items[3] + '->' + items[4]] = int(items[0])
    
    emissionLogProbabilities = {}
    transitionLogProbabilities = {}
    
    #estimate emission probabilities
    for k, v in emissionCounts.items():
        state = k.split('->')[0]
        obs = k.split('->')[1]
        key = obs + '|' + state
        seenObservations[obs] = 1
        value = math.log(v / float(unigramCounts[state]))
        emissionLogProbabilities[key] = value
    
    #estimate transition probabilities
    for k,v in trigramCounts.items():
        #k = D,N,V 
        state_old = k.split('->')[0] + '->' + k.split('->')[1]
        #count of D,N,V over count of D,N
        value = math.log( v/ float(bigramCounts[state_old]))
        #probability_key = V/D,N
        probability_key = k.split('->')[2] + '|' + state_old
        transitionLogProbabilities[probability_key] = value
    
    return (emissionLogProbabilities, transitionLogProbabilities, seenObservations)

def getMaxEmission(emissions, observation):
    states = ['O', 'I-GENE']
    max_log_prob = -100
    max_state =''
    for state in states:
        if (emissions.has_key(observation+'|'+state) and emissions[observation+'|'+state] > max_log_prob):
            #print observation , state, emissions[observation+'|'+state]
            max_log_prob = emissions[observation+'|'+state]
            max_state = state
    #print observation, max_state, max_log_prob
    return max_state
    
if __name__ == "__main__":
    if len(sys.argv) == 4:
        countFile = sys.argv[1]
        testFile = sys.argv[2]
        outputFile = sys.argv[3]
        (ep, tp, so) = MLEstimates(countFile)
        writer = open(outputFile,'w')
        for line in open(testFile,'r').readlines():
            observation = line.strip()
            
            max_state = ''
            if observation.strip() != '':
                if (so.has_key(observation)):
                    max_state = getMaxEmission(ep, observation)
                else:
                    max_state = getMaxEmission(ep, '_RARE_')
            writer.write(observation +' ' + max_state +'\n')
        writer.flush()
        writer.close()
            

    
    
