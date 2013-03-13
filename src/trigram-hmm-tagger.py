'''
Created on Mar 10, 2013

@author: arenduchintala
'''
import math
import sys
from ViterbiState import ViterbiState
from BackTraceNode import BackTraceNode

unigramCounts = {}
bigramCounts = {}
trigramCounts = {}
emissionCounts = {}
seenObservations = {}
emissions = {}
transitions = {}
pi = {}

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

def unwrap(pi,pi_pos):
    
    return tags
        
def getViterbiProbability(possible_states_per_word, words):
    pi = {}
    pi['0,*,*'] = [0]
    pi_pos={}
    
    for k in range(1,len(words)+2):
        uvs_lp = {}
        for u in possible_states_per_word[k-1]:
            for v in possible_states_per_word[k]:
                max_p = float('-inf')
                new_bt_nodes = []
                for w in possible_states_per_word[k-2]:
                    print 'k=',k,'w=',w,'u=',u,'v=',v
                    if (transitions.has_key(v + '|' + w + ',' + u)):
                        q = transitions[v + '|' + w + ',' + u]
                    else:
                        q = float("-inf")
                    if(v == 'STOP'):
                        e = 0.0
                    elif (emissions.has_key(words[k] + '|' + v)):
                        e = emissions[words[k] + '|' + v]
                    else:
                        e = float("-inf")
                   
                    pi_key = str(k-1)+','+w+','+u
                    print 'searching pi_key',pi_key
                    p = max(pi[str(int(k-1))+','+w+','+u]) + q + e
                    #if(p > max_p):
                    #    max_p = p
                    #    arg_pi[k-2]= w
                        
                    new_pi_key = str(k)+','+u+','+v
                    print 'adding pi_key',new_pi_key, 'p=',p, 'pp=', math.exp(p)
                    if (pi.has_key(new_pi_key)):
                        pi[new_pi_key].append(p)
                    else:
                        pi[new_pi_key] = [p]
                    if (pi_pos.has_key(k)):
                        pi_pos[k].append(new_pi_key)
                    else:
                        pi_pos[k] = [new_pi_key]
            
    


    print pi
    print pi_pos
    unwrap(pi,pi_pos)
    return words.values()


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
    
    possible_states_per_word[i] = ['STOP']
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
        writer = open(outputFile,'w')
        MLEstimates(countFile)
        testSentences = open(testFile, 'r').read().split("\n\n")
        for a_sentence in testSentences:
            (possible_states_per_word, words) = getPossibleSates(a_sentence.strip())
            #print a_sentence.strip().split("\n")
            tags = getViterbiProbability(possible_states_per_word, words)
            original_words = a_sentence.strip().split('\n')
            for word,tag in zip(original_words,tags):
               writer.write(word+' '+tag+'\n')
            writer.write('\n')    
        writer.flush()
        writer.close()
       
