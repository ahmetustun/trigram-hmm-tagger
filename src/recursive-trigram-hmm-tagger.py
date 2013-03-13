'''
Created on Mar 11, 2013

@author: arenduchintala
'''

import math
import sys

unigramCounts = {}
bigramCounts = {}
trigramCounts = {}
emissionCounts = {}
seenObservations = {}
emissions = {}
transitions = {}
stack = []

def MLEstimates(countFilePath):
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
            bigramCounts[items[2] + ',' + items[3]] = int(items[0])
        elif (items[1].strip() == '3-GRAM'):
            # 1813 3-GRAM I-GENE O STOP
            trigramCounts[items[2] + ',' + items[3] + ',' + items[4]] = int(items[0])

    # estimate emission probabilities
    for k, v in emissionCounts.items():
        state = k.split('->')[0]
        obs = k.split('->')[1]
        key = obs + '|' + state
        seenObservations[obs] = 1
        value = math.log(v / float(unigramCounts[state]))
        emissions[key] = value
    
    # estimate transition probabilities
    for k, v in trigramCounts.items():
        # k = D,N,V 
        state_old = k.split(',')[0] + ',' + k.split(',')[1]
        # count of D,N,V over count of D,N
        value = math.log(v / float(bigramCounts[state_old]))
        # probability_key = V/D,N
        probability_key = k.split(',')[2] + '|' + state_old
        transitions[probability_key] = value
    


def getWordPossibleTags(word):    
    states = []
    for state in unigramCounts.keys():
        if (emissions.has_key(word + '|' + state)):
            states.append(state)
    return states

def getRareWordCategory(word):
    return '_RARE_'



def pie(k, u, v, possible_states_per_word, words, maxargs):
   
    if  (u == '*' and v == '*'):
        return (0, maxargs)
    
    max = float("-Inf")
   
    for w in possible_states_per_word[k - 2]:
        if (transitions.has_key(v + '|' + w + ',' + u)):
            transition_lp = transitions[v + '|' + w + ',' + u]
        else:
            transition_lp = float("-Inf")
            
        (p, maxargs) = pie(k - 1, w, u, possible_states_per_word, words, maxargs)
        if (p + transition_lp >= max):
            max = p + transition_lp
             
            if (k - 2 > 0):
                maxargs[k - 2] = w
    
    if (emissions.has_key(words[k] + '|' + v)):
        emission_lp = emissions[words[k] + '|' + v]
    else:
        emission_lp = float("-Inf")
        
    return (max + emission_lp, maxargs)
        
def getViterbiSequences(possible_states_per_word, words):
    maxargs = {}
    for k in range(1,len(words)+1):
        uvs_lp = {}
        for v in possible_states_per_word[k]:
            for u in possible_states_per_word[k - 1]:
                (p, maxargs) = pie(k, u, v, possible_states_per_word, words, maxargs)
                print k, maxargs,p
                uvs_lp[u + ',' + v] = p
                

    max = float("-Inf")
    max_uv = ''
    for key, val in uvs_lp.items():
       if (transitions.has_key('STOP|' + key)):
           if (val + transitions['STOP|' + key] >= max):
               max = val + transitions['STOP|' + key]
               max_uv = key
       else:
           max = float("-Inf")
           max_uv = key

    maxargs[k - 1] = max_uv.split(',')[0]
    maxargs[k] = max_uv.split(',')[1]
    args = maxargs.values()
    
    print args, max, math.exp(max), '\n'
    return args


def getPossibleSates(sentence):
    words = {}
    possible_states_per_word = {}
    possible_states_per_word[-1] = ['*']
    possible_states_per_word[0] = ['*']
    i = 1
    for word in sentence.split("\n"):
        
        if (seenObservations.has_key(word)):
            words[i] = word
            possible_states_per_word[i] = getWordPossibleTags(word)
        else:
            rare_word = getRareWordCategory(word)
            words[i] = rare_word
            possible_states_per_word[i] = getWordPossibleTags(rare_word)
        # print word, possible_states_per_word[i]
        i += 1
    
    possible_states_per_word[i] = 'STOP'
    return (possible_states_per_word, words)

def getMaxEmission(emissions, observation):
    max_log_prob = float("-Inf")
    max_state = ''
    for state in unigramCounts.keys():
        if (emissions.has_key(observation + '|' + state) and emissions[observation + '|' + state] >= max_log_prob):
            # print observation , state, emissions[observation+'|'+state]
            max_log_prob = emissions[observation + '|' + state]
            max_state = state
    # print observation, max_state, max_log_prob
    return max_state
    
if __name__ == "__main__":
    if len(sys.argv) == 4:
        countFile = sys.argv[1]
        testFile = sys.argv[2]
        outputFile = sys.argv[3]
        MLEstimates(countFile)
        testSentences = open(testFile, 'r').read().split("\n\n")
        for a_sentence in testSentences:
            (possible_states_per_word, words) = getPossibleSates(a_sentence.strip())
            print a_sentence.strip().split("\n")
            max_porb_sequence = getViterbiSequences(possible_states_per_word, words)
            print '\n'
       
