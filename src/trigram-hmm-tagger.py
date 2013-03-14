'''
Created on Mar 10, 2013

@author: arenduchintala
'''
import math
import sys,re
from ViterbiState import ViterbiState
from pprint import pprint
from BackTraceNode import BackTraceNode

unigramCounts = {}
bigramCounts = {}
trigramCounts = {}
emissionCounts = {}
seenObservations = {}
emissions = {}
transitions = {}
pi = {}
s2i = {}
i2s ={}

def populateStateMap(_state):
    s2i[_state] = _state
    return _state
    #if (_state == '*' or _state == 'STOP'):
    #    s2i[_state] = _state
    #elif (not s2i.has_key(_state)):
    #    s2i[_state] = str(len(s2i.keys()))
    #return s2i[_state]
                
def MLEstimates(countFilePath):
    populateStateMap('*')
    populateStateMap('STOP')
    
    lines = open(countFilePath, 'r').readlines()
    for line in lines:
        items = line.strip().split(" ")
        if (items[1].strip() == 'WORDTAG'):
            # 8 WORDTAG O localizes
            emissionCounts[populateStateMap(items[2]) + '->' + items[3]] = int(items[0])
        elif (items[1].strip() == '1-GRAM'):
            # 749 1-GRAM I-GENE
            
            unigramCounts[s2i[items[2]]] = int(items[0])
            
                
        elif (items[1].strip() == '2-GRAM'):
            # 749 2-GRAM * I-GENE
            bigramCounts[s2i[items[2]] + ',' + s2i[items[3]]] = int(items[0])
        elif (items[1].strip() == '3-GRAM'):
            # 1813 3-GRAM I-GENE O STOP
            trigramCounts[s2i[items[2]] + ',' + s2i[items[3]] + ',' + s2i[items[4]]] = int(items[0])

    

 

    # estimate emission probabilities
    for k, v in emissionCounts.items():
        state = k.split('->')[0].strip()
        obs = k.split('->')[1].strip()
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

    
def unwrap(pi, pi_pos):
    max_pos_end = 'STOP'
    max_pos_mid = []
    for pk in sorted(pi_pos.keys(), reverse=True):
       print pk, pi_pos[pk]
       max_at_pk = float("-inf")
       max_pos = ''
       for pos in pi_pos[pk]:
           print pos, pi[pos]
           if (pi[pos] > max_at_pk and pos.split(",").pop() == max_pos_end):
               max_pos = pos
       print max_pos
       max_pos_end = max_pos.split(",")[1] 
       max_pos_mid.append(max_pos_end)
    max_pos_mid.pop()
    max_pos_mid.reverse()
    return max_pos_mid
        
def getViterbiProbability(possible_states_per_word, words):
    pi = {}
    arg_pi = {}
    pi['0,*,*'] = 0 
    arg_pi['0,*,*'] = []
    
    for k in range(1, len(words) + 2):
        for u in possible_states_per_word[k - 1]:
          
            for v in possible_states_per_word[k]:
            
               

                ps = []
                pw = []
                for w in possible_states_per_word[k - 2]:
                   
                    print 'k=', k, 'w=', w, 'u=', u, 'v=', v
                    if (transitions.has_key(v + '|' + w + ',' + u)):
                        q = transitions[v + '|' + w + ',' + u]
                    else:
                        q = float("-inf")
                    if(v == 'STOP'):
                        e = 0.0
                    elif (emissions.has_key(words[k] + '|' + v)):
                        # print 'e = ',words[k],'|',v
                        e = emissions[words[k] + '|' + v]
                    else:
                        e = float("-inf")
                   
                    pi_key = str(k - 1) + ',' + w + ',' + u
                    print 'searching pi_key', pi_key
                    p = pi[str(k - 1) + ',' + w + ',' + u] + q + e
                    bt = list(arg_pi[str(k - 1) + ',' + w + ',' + u])
                    
                    bt.append(w)
                    if (v == 'STOP'):
                        bt.append(u)
                    ps.append(p)
                    pw.append(bt)
                max_p = max(ps)
                max_bt = pw[ps.index(max_p)]
                new_pi_key = str(k) + ',' + u + ',' + v
                # print 'pi[',new_pi_key, ']= pi[',pi_key,']=',math.exp(pi[str(int(k-1))+','+w+','+u]),' q=',math.exp(q),' * e=',math.exp(e), ' = ' ,math.exp(max_p)
                pi[new_pi_key] = max_p
                arg_pi[new_pi_key] = max_bt

            
    
    max_bt.pop(0)
    max_bt.pop(0)
    max_bt = map(lambda bt: i2s[bt], max_bt)
    print "back trace:", max_bt
    print "probability", max_p
    return max_bt



def getPossibleSates(sentence):
    words = {}
    possible_states_per_word = {}
    possible_states_per_word[-1] = ['*']
    possible_states_per_word[0] = ['*']
    i = 1
    for word in sentence.split("\n"):
        
        word = word.strip()
       
        if (seenObservations.has_key(word)):
            
            words[i] = word
            possible_states_per_word[i] = getWordPossibleTags(word)
        else:
            rare_word = replaceRareWord(word)
            print rare_word
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

    countFile = 'outputs/gene-rare.count'
    testFile = 'data/gene.test'
    outputFile = 'outputs/gene_test.p2.out'
    writer = open(outputFile, 'w')
    MLEstimates(countFile)
    i2s = {v:k for k, v in s2i.items()}  # create a number to state mapping for displaying the state sequence in the end
    testSentences = open(testFile, 'r').read().split("\n\n")
    for a_sentence in testSentences:
        (possible_states_per_word, words) = getPossibleSates(a_sentence.strip())
        print a_sentence.strip().split("\n")
        tags = getViterbiProbability(possible_states_per_word, words)
        original_words = a_sentence.strip().split('\n')
        for word, tag in zip(original_words, tags):
            print word, tag
            writer.write(word + ' ' + tag + '\n')
        writer.write('\n')    
    writer.flush()
    writer.close()
       
